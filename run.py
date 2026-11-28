import argparse
from copy import deepcopy
import glob
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time
from src.common.util.MetaUtil import MetadataUtils
from src.common.util.LogUtil import setup_logger
from src.common.util.CmdUtil import _run_cmd, _run_cmds_parallel
from src.common.util.ProcessUtil import SignalExit, terminate_all_children
from src.common.util.SchemaValidatorUtil import SchemaValidator, smart_cast
from src.common.util.EnvUtil import is_path_like
from src.common.util.CallbackUtil import ExecutionCallback, DEFAULT_CALLBACK_URL
from node import runCoCulture, runMERIP, runRNAseq, runncRNAseq, runCLIP, runMutation, runPacVar, runKARRseq, runPeakCalling, runRIPseq, runQuantMS, runtRNAseq, runscRNAseq, runFiberseq, runLRtranscriptome, runDR, runPopulation_genomics, runspatial_transcriptomics
import logging
from typing import Dict, Any, List, Optional, Tuple, Set, Union
logger = setup_logger(__name__, level=logging.DEBUG)

def dict_set_by_path(d: Dict[str, Any], keys: List[str], value: Any) -> None:
    """Set a value in a nested dict by a sequence of keys.

    Creates intermediate dicts as needed. The value is passed through
    smart_cast for automatic type conversion.

    Args:
        d: Target dict (modified in-place).
        keys: Sequence of keys forming the path (e.g. ["genome", "fasta"]).
        value: Value to set at the nested location.
    """
    for k in keys[:-1]:
        if k not in d or not isinstance(d[k], dict):
            d[k] = {}
        d = d[k]
    d[keys[-1]] = smart_cast(value)

def parse_dot_args(extra_args: Dict[str, Any]) -> Dict[Tuple[str, ...], Any]:
    """Extract dot-notation keys from extra_args into nested key tuples.

    Converts {"genome.fasta": "/path"} → {("genome", "fasta"): "/path"}
    for use with dict_set_by_path.

    Args:
        extra_args: Dict of extra CLI arguments.

    Returns:
        Dict mapping key tuples to values for dot-notation arguments.
    """
    dot_args = {}
    for k, v in list(extra_args.items()):
        if '.' in k:
            dot_args[tuple(k.split('.'))] = v
    return dot_args


def _load_model_json(model_json_file: str) -> Dict[str, Any]:
    """Load model JSON template from disk."""
    with open(model_json_file, 'r', encoding='utf-8') as f:
        return json.load(f)



def _detect_total_mem_mb(default_mb: int = 470000) -> int:
    """Total RAM in MB for the snakemake --resources budget.

    Reads CGroup V2 limit first (containers/slurm set it), then /proc/meminfo.
    Budget is 90% of total (reserves ~10% for the OS). Falls back to
    default_mb when neither is readable.
    """
    try:
        with open("/sys/fs/cgroup/memory.max") as fh:
            v = fh.read().strip()
        if v != "max":
            total_mb = int(v) // (1024 * 1024)
            if total_mb > 0:
                return int(total_mb * 0.90)
    except (OSError, ValueError):
        pass
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    total_mb = int(line.split()[1]) // 1024
                    return int(total_mb * 0.90)
    except (OSError, ValueError, IndexError):
        pass
    return default_mb


MIN_MEM_SHARE_MB = 8000

# I/O 带宽预算(设计 §8.3):盘类型 → 可持续混合读写带宽(MB/s)。全部为待校准
# 暂定值(HDD 取现场 sda 饱和吞吐 76~100 MB/s 留余量;SSD/NVMe 为数量级估值),
# 由校准扩展(设计 §10)实测重定。
DISK_SUSTAINABLE_MB_S = {"hdd": 150, "sata_ssd": 450, "nvme": 1200}
IO_SAFETY_FACTOR = 0.8
MIN_DISK_SHARE_MB_S = 10   # 多工作流切分时每工作流带宽份额下限(经验值,待校准)
TMP_FREE_WARN_MB = 50000  # tmp 剩余空间告警阈值(经验值:STAR 单 job tmp ~17 GB,待校准)


def _detect_disk_budget_mb_s(outdir: str) -> int:
    """I/O 带宽预算(MB/s) = io_safety_factor × outdir 所在盘的可持续带宽(设计 §8.3/§8.4)。

    盘类型探测:st_dev → /sys/dev/block/<maj>:<min> 的 DEVNAME(含 nvme → NVMe)
    与 queue/rotational(分区缺 queue 时向父盘回退一层)。任何探测失败回退 HDD
    暂定常量并告警,绝不拒绝启动(设计 §8.6)。
    """
    kind = "hdd"
    try:
        st_dev = os.stat(outdir).st_dev
        dev_dir = os.path.realpath(
            f"/sys/dev/block/{os.major(st_dev)}:{os.minor(st_dev)}")
        name = os.path.basename(dev_dir)
        try:
            with open(os.path.join(dev_dir, "uevent")) as fh:
                for line in fh:
                    if line.startswith("DEVNAME="):
                        name = line.strip().split("=", 1)[1]
                        break
        except OSError:
            pass
        rotational = None
        for d in (dev_dir, os.path.dirname(dev_dir)):
            try:
                with open(os.path.join(d, "queue", "rotational")) as fh:
                    rotational = fh.read().strip()
                break
            except OSError:
                continue
        if "nvme" in name:
            kind = "nvme"
        elif rotational == "0":
            kind = "sata_ssd"
    except OSError:
        logger.warning("disk budget: cannot probe outdir device, "
                       "falling back to HDD constant")
    return int(DISK_SUSTAINABLE_MB_S[kind] * IO_SAFETY_FACTOR)


def _check_tmp_free_space(tmp_root: str) -> None:
    """落 tmp 前的可用空间检查(设计 §9 边界):不足只告警不拦截,失败由重试兜底。"""
    try:
        st = os.statvfs(tmp_root)
        free_mb = (st.f_bavail * st.f_frsize) // (1024 * 1024)
        if free_mb < TMP_FREE_WARN_MB:
            logger.warning(
                "tmp dir %s has only %d MB free (< %d MB): concurrent jobs may "
                "fill it; job failures fall back to --restart-times retries",
                tmp_root, free_mb, TMP_FREE_WARN_MB)
    except OSError:
        pass


def detect_host_profile() -> Dict[str, Any]:
    """Detect {hostname, cpus, mem_mb} of the current host (resource-balancing-proposal §3 第 1 点).

    cpus: sched affinity first (containers/taskset bind below os.cpu_count),
    then cgroup v2 cpu.max quota, then os.cpu_count(). mem_mb reuses
    _detect_total_mem_mb() (cgroup memory.max -> /proc/meminfo, 90% of total).
    """
    hostname = socket.gethostname()
    try:
        cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        cpus = os.cpu_count() or 4
    try:
        with open("/sys/fs/cgroup/cpu.max") as fh:
            quota_s, period_s = fh.read().split()
        if quota_s != "max":
            cpus = min(cpus, max(1, int(quota_s) // int(period_s)))
    except (OSError, ValueError):
        pass
    return {"hostname": hostname, "cpus": cpus, "mem_mb": _detect_total_mem_mb()}


def _workflow_input_bytes(input_json: str) -> Optional[int]:
    """Total bytes of existing files referenced by a workflow config JSON.

    Workload weight for the multi-workflow memory split (resource-balancing-proposal §3 第 3 点).
    Returns None when no referenced file exists yet so callers fall back to
    an equal split.
    """
    try:
        with open(input_json, "r", encoding="utf-8") as f:
            cfg = json.load(f)
    except (OSError, ValueError):
        return None
    total = 0
    seen: Set[str] = set()

    def _walk(obj: Any) -> None:
        nonlocal total
        if isinstance(obj, str):
            if is_path_like(obj) and obj not in seen:
                seen.add(obj)
                try:
                    total += os.path.getsize(obj)
                except OSError:
                    pass
        elif isinstance(obj, dict):
            for v in obj.values():
                _walk(v)
        elif isinstance(obj, list):
            for v in obj:
                _walk(v)

    _walk(cfg)
    return total if total > 0 else None


def _split_mem_budget(total_mb: int, weights: List[Optional[int]]) -> List[int]:
    """Split the memory budget across workflows (resource-balancing-proposal §3 第 3 点,已定:加权).

    Weighted by input bytes; equal split when weights are unavailable.
    Every share keeps at least min(total//n, MIN_MEM_SHARE_MB) so a tiny
    workflow cannot starve below usability.
    """
    n = len(weights)
    if n <= 1:
        return [total_mb]
    w = [x or 0 for x in weights]
    if sum(w) <= 0:
        return [total_mb // n] * n
    shares = [int(total_mb * x / sum(w)) for x in w]
    floor = min(total_mb // n, MIN_MEM_SHARE_MB)
    return [max(s, floor) for s in shares]


def _split_disk_budget(total_mb_s: int, weights: List[Optional[int]]) -> List[int]:
    """Split the disk bandwidth budget (MB/s) across workflows sharing one device.

    Same rules as _split_mem_budget (input-bytes weighted, equal split when
    weights unavailable, floor share so tiny workflows keep a usable slice) —
    设计 §8.4 第 3 点.
    """
    n = len(weights)
    if n <= 1:
        return [total_mb_s]
    w = [x or 0 for x in weights]
    if sum(w) <= 0:
        return [total_mb_s // n] * n
    shares = [int(total_mb_s * x / sum(w)) for x in w]
    floor = min(total_mb_s // n, MIN_DISK_SHARE_MB_S)
    return [max(s, floor) for s in shares]


def _inject_resource_caps(input_json: str, mem_mb: int, cpus: int,
                          disk_budget_mb_s: Optional[int] = None) -> None:
    """Add config["resource_caps"] to a generated workflow config JSON.

    common.smk's rule_res clamps declared mem_mb/threads/io_bandwidth_mb_s
    against these caps so no job can declare more than the process budget.
    The values MUST equal the --resources values passed on the same snakemake
    command line. disk_budget_mb_s (设计 §8.5) clamps io_bandwidth_mb_s;
    omitted/0 -> no io clamping (run_module.py standalone behaviour).
    """
    with open(input_json, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    caps = {"mem_mb": int(mem_mb), "cpus": int(cpus)}
    if disk_budget_mb_s:
        caps["disk_bandwidth_budget_mb_s"] = int(disk_budget_mb_s)
    cfg["resource_caps"] = caps
    with open(input_json, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)


def _write_run_summary(outdir: str, payload: Dict[str, Any]) -> None:
    """Write run_summary.json (效果评判数据来源, proposal §7.4). Best effort."""
    try:
        path = os.path.join(outdir, "run_summary.json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        logger.info("run summary written: %s", path)
    except OSError as e:
        logger.warning("run summary not written: %s", e)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments and collect extra dot-notation overrides.

    Handles three forms of extra arguments:
      --key=value, --key value, --key v1 v2 v3

    Dot-notation keys (e.g. --genome.fasta /path) are stored in args.extra_args
    for later injection into the workflow config via dict_set_by_path.

    Returns:
        argparse.Namespace with all parsed arguments plus extra_args dict.
    """
    parser = argparse.ArgumentParser(description="workflow")
    parser.add_argument('-m','--meta', type=str, default=None, help='meta input file or data dir which condatain fastq file')
    parser.add_argument('--design', '--design-file', dest='design', type=str, default=None,
        help='experiment design JSON (sample_params.design shape: per-tissue '
             'groups + comparisons with groups/ref/coding). The single entry for '
             'experiment comparison design; loaded into sample_params.design and '
             'derived into group_pairs / designPairs')
    parser.add_argument('-w','--workflow_name', type=str, nargs='+',
        choices=["all", "CoCulture", "MERIP", "RNAseq", "ncRNAseq", "scRNAseq", "tRNAseq", "CLIP", "Mutation", "PacVar", "KARRseq", "PeakCalling", "QuantMS", "Fiberseq", "LRtranscriptome", "DR", "RIPseq", "Population_genomics", "spatial_transcriptomics"],
        default=['CoCulture'], help='workflow name(s), multiple for parallel execution; "all" = every registered workflow')
    parser.add_argument('-o','--output_dir','--output-dir', dest='output_dir', type=str, default=None, help='output dir')
    parser.add_argument('-t','--threads', type=int, default=None,
        help='total threads passed to snakemake --cores (split across workflows); '
             'default: auto = 90%% of host CPUs, at least 1 core kept free')
    parser.add_argument('--dry-run', action='store_true',
        help='dry run: build the DAG and exit without executing run-block code')
    parser.add_argument('--touch', action='store_true', help='touch output files to update timestamps without re-running (cannot combine with --dry-run)')
    parser.add_argument('--clean-temp', action='store_true',
                        help='clean the run outdir(s) temp space (tmp/ plus container/tool cache dirs) and exit')
    parser.add_argument('--test', type=str, nargs='?', const='dag', choices=['dag', 'shell'], metavar='MODE',
        help='test mode with built-in fixtures; workflow(s) via -w. '
             'dag = only build the DAG; shell = execute run-block code and write '
             "each rule's shell receipt without running tools (see modules.md). "
             'Default MODE: dag')
    parser.add_argument('--log', type=str, default='workflow.log', help='log file')
    parser.add_argument('--callback-url', type=str, default=None,
        help='execution callback endpoint (default: env BIPLATFORM_CALLBACK_URL '
             'or ' + DEFAULT_CALLBACK_URL + ')')
    parser.add_argument('--callback-execution-id', type=str, default=None,
        help='execution id reported to the callback endpoint (default: env '
             'BIPLATFORM_EXEC_ID)')
    parser.add_argument('--callback-token', type=str, default=None,
        help='per-execution callback token (default: env '
             'BIPLATFORM_CALLBACK_TOKEN)')
    parser.add_argument('--no-callback', action='store_true',
        help='disable the execution callback even when id/token are configured')
    parser.add_argument('--conda-prefix', type=str, default=None, help='conda prefix for snakemake (required when NOT using --sdm apptainer)')
    parser.add_argument('--sdm', action='store_const', const='apptainer', default=None,
        help='use apptainer container backend (SIF images). When set, --use-conda is omitted')
    parser.add_argument('--mem-mb', type=int, default=None, dest='mem_mb',
        help='global mem_mb budget passed to snakemake --resources; jobs declaring '
             'resources.mem_mb are scheduled against it. Default: auto-detect '
             '(90%% of total RAM)')
    parser.add_argument('--io-budget-mb-s', type=int, default=None, dest='io_budget_mb_s',
        help='total disk bandwidth budget (MB/s) passed to snakemake --resources as '
             'disk_bandwidth_budget_mb_s; jobs declaring resources.io_bandwidth_mb_s '
             'are packed against it (design §8). Default: auto-detect from each '
             "output disk's type; 0 disables I/O-aware scheduling")
    parser.add_argument('--singularity-args', type=str, default=None,
        help='arguments for singularity/apptainer, e.g. --singularity-args \'--bind /path1,/path2\' (use with --sdm)')
    parser.add_argument(
        '--rerun-trigger', '--rerun-triggers',
        dest='rerun_trigger',
        nargs='+',
        default=["mtime"],
        choices=["code", "input", "mtime", "params", "software-env"],
        help='snakemake rerun-triggers, e.g. code input mtime params software-env'
    )
    parser.add_argument('--conda-frontend', type=str, choices=["conda", "mamba"], default="mamba", help='conda frontend for snakemake')
    parser.add_argument('--forcerun', type=str, nargs='+', default=None,
        help='force re-run specific rules or files; a bare RULE name runs ALL jobs of '
             'that rule (wildcards filled automatically), format: RULE, RULE:wildcard=value, '
             'or /path/to/file, e.g. --forcerun scanpy_qc')
    parser.add_argument('--forcerun-downstream', type=str, nargs='+', default=None,
        help='force re-run these rules/files AND continue through their downstream '
             'dependents (the whole workflow runs; unlike --forcerun which stops at '
             'the named rules), format: RULE, RULE:wildcard=value, or /path/to/file, '
             'e.g. --forcerun-downstream somatic_spectrum')
    parser.add_argument('--unlock', action='store_true', help='unlock snakemake working directory')
    parser.add_argument('--restart-times', dest='restart_times', type=int, default=2,
        help='snakemake --restart-times: retry failed jobs up to N times (resource-balancing-proposal §5; '
             'resources mem_mb attempt-escalation x2**(attempt-1) in rule_res)')
    parser.add_argument('--no-schema-validate', dest='schema_validate', action='store_false',
        default=True, help='disable schema-based type casting for extra args')
    parser.add_argument('--list-target-rules', dest='list_target_rules', action='store_true', help='list all target rules in the workflow and exit, no wildcard expansion')
    parser.add_argument('--list-rules', dest='list_rules', action='store_true', help='list all rules in the workflow and exit')
    parser.add_argument(
        '--snakemake-args',
        nargs=argparse.REMAINDER,
        default=[],
        help='additional arguments forwarded to snakemake; place them after this flag'
    )

    args, unknown = parser.parse_known_args()
    extra_args: Dict[str, Any] = {}
    i = 0
    while i < len(unknown):
        arg = unknown[i]
        if arg.startswith('--'):
            key = arg[2:]
            if '=' in key:
                k, v = key.split('=', 1)
                extra_args[k] = v
            else:
                values = []
                while i + 1 < len(unknown) and not unknown[i + 1].startswith('--'):
                    values.append(unknown[i + 1])
                    i += 1
                if not values:
                    extra_args[key] = True
                elif len(values) == 1:
                    extra_args[key] = values[0]
                else:
                    extra_args[key] = values
        i += 1
    args.extra_args = extra_args
    return args


def _detect_singularity(snakemake_args: Optional[List[str]]) -> bool:
    """Check if --sdm apptainer/singularity is present in raw snakemake_args (backward compat)."""
    snakemake_args = snakemake_args or []
    return any(
        arg in ("apptainer", "singularity")
        for i, arg in enumerate(snakemake_args)
        if i > 0 and snakemake_args[i - 1] == "--sdm"
    )


def _collect_bind_paths(config_path: str) -> List[str]:
    """Scan a config JSON file and collect directory paths to bind-mount.

    Recursively walks all string values; strings containing '/' are treated as
    paths. Files yield their parent directory; directories are kept as-is.
    Subdirectories are collapsed into their parents to minimise the bind list.
    """
    import json as _json
    import os as _os

    with open(config_path, "r", encoding="utf-8") as f:
        cfg = _json.load(f)

    dirs: Set[str] = set()

    def _walk(obj: Any) -> None:
        if is_path_like(obj) and os.path.exists(obj):
            p = _os.path.abspath(obj.strip())
            if _os.path.isfile(p):
                dirs.add(_os.path.dirname(p))
            else:
                dirs.add(p)
        elif isinstance(obj, dict):
            for v in obj.values():
                _walk(v)
        elif isinstance(obj, list):
            for v in obj:
                _walk(v)

    _walk(cfg)

    result = set(dirs)
    for d in list(dirs):
        for other in dirs:
            if d != other and d.startswith(other.rstrip("/") + "/"):
                result.discard(d)
                break
    return sorted(result)

def _merge_singularity_args(
    json_bind_paths: List[str],
    singularity_args: Optional[str] = None,
    tmp_dir: Optional[str] = None,
) -> Optional[str]:
    """
    Merge bind paths collected from the JSON configuration with
    user-provided Singularity arguments.

    JSON-derived bind paths are always preserved. User-provided
    bind paths are added to them. Duplicate paths are removed while
    preserving their original order.

    Other Singularity arguments, such as --cleanenv, are preserved.
    """
    bind_paths: List[str] = []

    def add_bind_paths(paths: List[str]) -> None:
        for path in paths:
            path = path.strip()
            if path and path not in bind_paths:
                bind_paths.append(path)

    add_bind_paths(json_bind_paths)

    if tmp_dir:
        # 规则以 tmp_dir 绝对路径放临时文件(设计 §9),容器内需同名可见,
        # 故自绑定;/tmp 映射保留给默认写 /tmp 的工具。
        add_bind_paths([tmp_dir, f"{tmp_dir}:/tmp"])
    else:
        add_bind_paths(["/tmp"])

    if not singularity_args:
        return (
            "--bind " + ",".join(bind_paths)
            if bind_paths
            else None
        )

    match = re.search(
        r"(?:^|\s)--bind(?:=|\s+)([^\s]+)",
        singularity_args,
    )

    if match:
        user_bind = match.group(1)

        add_bind_paths(user_bind.split(","))

        merged_bind = "--bind " + ",".join(bind_paths)

        start, end = match.span()
        singularity_args = (
            singularity_args[:start]
            + merged_bind
            + singularity_args[end:]
        )
    else:
        singularity_args = (
            singularity_args.rstrip()
            + " --bind "
            + ",".join(bind_paths)
        )

    return singularity_args.strip()

def build_snakemake_cmd(root_dir: str, smk: str, input_json: str, threads: int,
                        conda_prefix: Optional[str], rerun_trigger: List[str],
                        dry_run: bool, conda_frontend: str, snakemake_args: Optional[List[str]],
                        sdm: Optional[str] = None, singularity_args: Optional[str] = None,
                        forcerun: Optional[List[str]] = None, unlock: bool = False,
                        touch: bool = False, list_target_rules: bool = False,
                        list_rules: bool = False,
                        forcerun_downstream: Optional[List[str]] = None,
                        mem_mb: Optional[int] = None,
                        disk_budget_mb_s: Optional[int] = None,
                        restart_times: Optional[int] = None) -> List[str]:
    """Build the snakemake CLI command list for a given subworkflow.

    Configures conda or apptainer container backend, collects bind paths
    from the config JSON, merges user singularity args, and appends
    forcerun targets if specified. forcerun targets scope the DAG to the
    named rules and stop there; forcerun_downstream targets additionally
    continue through their dependents (whole-workflow scope).

    Args:
        root_dir: Project root directory (contains subworkflow/ and config/).
        smk: Subworkflow snakefile name (e.g. "RNAseq.smk").
        input_json: Path to the workflow config JSON file.
        threads: Number of threads for snakemake.
        conda_prefix: Conda prefix directory for snakemake --use-conda.
        rerun_trigger: List of snakemake rerun-trigger reasons.
        dry_run: If True, append --dry-run.
        conda_frontend: Conda frontend ("conda" or "mamba").
        snakemake_args: Additional raw arguments forwarded to snakemake.
        sdm: Snakemake deployment method (e.g. "apptainer"). None for conda mode.
        singularity_args: Extra singularity/apptainer arguments (e.g. "--bind /path1,/path2").
        forcerun: List of rule names to force re-run (with optional wildcards).
        forcerun_downstream: Like forcerun, but also re-run the rules'
            downstream dependents (whole-workflow scope).
        unlock: If True, append --unlock.
        touch: If True, append --touch.
        list_target_rules: If True, append --list-target-rules.
        list_rules: If True, append --list-rules.
    Returns:
        List[str]: Complete snakemake command as a list of arguments.
    """
    snakemake_args = snakemake_args or []
    use_singularity = sdm is not None or _detect_singularity(snakemake_args)

    cmd = [
        "snakemake",
        "-s",
        f"{root_dir}/subworkflow/{smk}",
        "--configfile",
        input_json,
        "--cores",
        str(threads),
        "--rerun-triggers",
        *rerun_trigger,
    ]
    resources = []
    if mem_mb:
        resources.append(f"mem_mb={mem_mb}")
    if disk_budget_mb_s:
        resources.append(f"disk_bandwidth_budget_mb_s={disk_budget_mb_s}")
    if resources:
        cmd += ["--resources", *resources]
    if disk_budget_mb_s:
        # 未声明 io_bandwidth_mb_s 的规则的地板权重(设计 §8.6/§8.7);
        # 限流关闭(0/未设)时不加,行为与现状逐比特一致。
        cmd += ["--default-resources", "io_bandwidth_mb_s=1"]
    if restart_times:
        cmd += ["--restart-times", str(restart_times)]
    # 临时目录(设计 §9):tmp_dir 未配置时默认 /tmp(自动分流,2026-11-27 已决)。
    # 容器内同名可见:_merge_singularity_args 绑定 tmp_dir 自身及其 /tmp 映射。
    with open(input_json, "r", encoding="utf-8") as _f:
        _cfg = json.load(_f)
    tmp_dir = _cfg.get("tmp_dir")
    if tmp_dir and not os.path.isabs(tmp_dir):
        tmp_dir = os.path.join(os.path.dirname(input_json), tmp_dir)
    if tmp_dir:
        os.makedirs(tmp_dir, exist_ok=True)
    _check_tmp_free_space(tmp_dir or "/tmp")

    if os.environ.get("OMICS_SCRIPT_ONLY", "") == "1":
        pass
    elif not use_singularity:
        cmd += [
            "--conda-prefix",
            conda_prefix,
            "--use-conda",
            "--conda-frontend",
            conda_frontend,
        ]
    else:
        if sdm is not None:
            cmd += ["--sdm", sdm]

        json_bind_paths = _collect_bind_paths(input_json)

        singularity_args = _merge_singularity_args(
            json_bind_paths,
            singularity_args,
            tmp_dir=tmp_dir,
        )

        if singularity_args:
            cmd += [
                "--singularity-args",
                singularity_args,
            ]
    if dry_run:
        cmd.append("--dry-run")
    if touch:
        cmd.append("--touch")
    if unlock:
        cmd.append("--unlock")
    if forcerun or forcerun_downstream:
        forced = list(forcerun or []) + list(forcerun_downstream or [])
        plain_targets = [t for t in forced if not (":" in t and "=" in t)]
        wildcard_targets = [t for t in forced if ":" in t and "=" in t]
        if forcerun_downstream:
            cmd[1:1] = ["all"]
        elif plain_targets:
            cmd[1:1] = plain_targets
        if plain_targets:
            cmd.append("--forcerun")
            cmd.extend(plain_targets)
        if wildcard_targets:
            cmd.append("--target-jobs")
            cmd.extend(wildcard_targets)
            cmd.append("--force")
    if list_target_rules:
        cmd.append("--list-target-rules")
    if list_rules:
        cmd.append("--list-rules")
    if snakemake_args:
        cmd.extend(snakemake_args)
    return cmd

def _prefix_forcerun_targets(targets: List[str], wf_name: str) -> List[str]:
    """Auto-prefix forcerun-style targets with the workflow name.

    Supports three formats:
      - File paths (/path/to/file) — passed through as-is
      - Rule:wildcard syntax (rule:wildcard=value) — prefix rule name only
      - Plain rule names — prefix with workflow name ("all" and already
        prefixed names pass through untouched).

    Args:
        targets: Raw target strings from --forcerun / --forcerun-downstream.
        wf_name: Workflow name used as the rule prefix.

    Returns:
        Prefixed target list, order preserved.
    """
    prefixed: List[str] = []
    for t in targets:
        if "/" in t:
            prefixed.append(t)
            continue
        if ":" in t and "=" in t:
            rule_part, wildcards = t.split(":", 1)
            if not rule_part.startswith(f"{wf_name}_"):
                rule_part = f"{wf_name}_{rule_part}"
            prefixed.append(f"{rule_part}:{wildcards}")
            continue
        if ":" in t:
            rule_part, rest = t.split(":", 1)
            rest = ":" + rest
        else:
            rule_part, rest = t, ""
        if rule_part == "all" or rule_part.startswith(f"{wf_name}_"):
            prefixed.append(t)
        else:
            prefixed.append(f"{wf_name}_{rule_part}{rest}")
    return prefixed

def _enumerate_rule_instances(
    root_dir: str, smk: str, input_json: str, rule: str, logger,
) -> List[Dict[str, str]]:
    """Enumerate every DAG job of a rule via a dry-run probe.

    Runs ``snakemake --dry-run --until R --forcerun R`` (a single-rule
    --until set cannot hit snakemake's #823 deadlock) and parses each job
    block's ``wildcards:`` line. Forced jobs are listed even when their
    outputs are up to date, so all wildcard instances of the rule appear.
    Returns one wildcard dict per job (empty dict when the rule has no
    wildcards).

    Raises:
        RuntimeError: snakemake missing from PATH, or the probe failed —
            its output is the DAG-build error the real run would hit.
        ValueError: the rule has no jobs in the current DAG.
    """
    if shutil.which("snakemake") is None:
        raise RuntimeError("snakemake not found on PATH (needed to expand rule targets)")
    probe = [
        "snakemake", "-s", f"{root_dir}/subworkflow/{smk}",
        "--configfile", input_json, "--cores", "1",
        "--rerun-triggers=mtime", "--dry-run",
        "--until", rule, "--forcerun", rule,
    ]
    probe_cwd = tempfile.mkdtemp(prefix="omics_target_probe_")
    try:
        result = subprocess.run(probe, capture_output=True, text=True, cwd=probe_cwd)
    finally:
        shutil.rmtree(probe_cwd, ignore_errors=True)
    out = (result.stdout or "") + "\n" + (result.stderr or "")
    if result.returncode != 0:
        raise RuntimeError(
            f"target expansion probe failed for rule '{rule}' "
            f"(rc={result.returncode}):\n{out[-3000:]}")
    instances: List[Dict[str, str]] = []
    for block in re.findall(rf"(?ms)^rule {re.escape(rule)}:\n(.*?)(?=^\S|\Z)", out):
        inst: Dict[str, str] = {}
        m = re.search(r"(?m)^\s+wildcards: (.*)$", block)
        if m and m.group(1).strip():
            for pair in m.group(1).strip().split(", "):
                key, sep, value = pair.partition("=")
                if not sep or not re.fullmatch(r"\w+", key):
                    raise RuntimeError(
                        f"cannot parse probe wildcards for rule '{rule}': {m.group(1)!r}")
                inst[key] = value
        instances.append(inst)
    if not instances:
        raise ValueError(
            f"rule '{rule}' has no jobs in the current DAG — check the rule "
            f"name (run.py --list-rules) and the workflow config")
    return instances


def _expand_wildcard_rule_targets(
    targets: List[str],
    wf_name: str,
    root_dir: str,
    smk: str,
    input_json: str,
    logger,
) -> List[str]:
    """Expand bare rule names into rule:wildcard=value targets.

    Snakemake rejects a wildcarded rule as a CLI target ("Target rules may
    not contain wildcards") and --target-jobs cannot enumerate instances,
    so a bare wildcarded rule name is expanded into concrete
    ``rule:wildcard=value`` specs (one per DAG job of the rule, wildcards
    filled automatically) before it becomes a target — a bare rule name
    means "run all jobs of this rule". File paths, explicit
    ``rule:wildcard=value`` specs, ``all`` and rules without wildcards pass
    through untouched. Order is preserved, duplicates removed.

    Args:
        targets: Prefixed forcerun targets (bare names already carry the
            ``{wf_name}_`` prefix).
        wf_name: Workflow name (used for logging only).
        root_dir: Workflow root (contains subworkflow/).
        smk: Subworkflow snakefile name.
        input_json: Path to the raw.json.
        logger: Logger for diagnostics.

    Returns:
        Target list with bare wildcarded rules replaced by their
        ``rule:wildcard=value`` expansions (order preserved, duplicates
        removed).

    Raises:
        RuntimeError / ValueError: see _enumerate_rule_instances.
    """
    if not targets:
        return targets

    expanded: List[str] = []
    for t in targets:
        if (":" in t or "=" in t or "/" in t or t == "all"
                or os.path.exists(t)):
            expanded.append(t)
            continue
        instances = _enumerate_rule_instances(root_dir, smk, input_json, t, logger)
        wildcarded = [inst for inst in instances if inst]
        if not wildcarded:
            expanded.append(t)
            continue
        for inst in wildcarded:
            expanded.append(
                t + ":" + ",".join(f"{k}={v}" for k, v in sorted(inst.items())))
        logger.info(
            "[%s] Expanded bare rule '%s' into %d wildcard targets: %s",
            wf_name, t, len(wildcarded), expanded[-len(wildcarded):],
        )

    seen = set()
    deduped = []
    for t in expanded:
        if t not in seen:
            seen.add(t)
            deduped.append(t)
    return deduped


WORKFLOW_DISPATCH = {
    "CoCulture":  lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("CoCulture.smk", runCoCulture(cfg, sid, indir, outdir, rf)),
    "MERIP":      lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("MERIP.smk",     runMERIP(cfg, sid, indir, outdir, rf)),
    "RNAseq":     lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("RNAseq.smk",    runRNAseq(cfg, sid, gp, indir, outdir, rf)),
    "ncRNAseq":   lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("ncRNAseq.smk",  runncRNAseq(cfg, sid, sp, indir, outdir, rf)),
    "CLIP":       lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("CLIP.smk",      runCLIP(cfg, sid, indir, outdir, rf)),
    "Mutation":   lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("Mutation.smk",  runMutation(cfg, sid, sp, indir, outdir, rf)),
    "PacVar":     lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("PacVar.smk",    runPacVar(cfg, sid, gp, indir, outdir, rf)),
    "KARRseq":    lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("KARRseq.smk",   runKARRseq(cfg, sid, indir, outdir, rf)),
    "PeakCalling":lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("PeakCalling.smk",runPeakCalling(cfg, sid, sp, indir, outdir, rf)),
    "RIPseq":     lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("RIPseq.smk",   runRIPseq(cfg, sid, sp, indir, outdir, rf)),
    "QuantMS":    lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("QuantMS.smk",   runQuantMS(cfg, sid, indir, outdir, rf)),
    "tRNAseq":    lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("tRNAseq.smk",   runtRNAseq(cfg, sid, indir, outdir, meta, rf)),
    "scRNAseq":   lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("scRNAseq.smk",  runscRNAseq(cfg, sid, indir, outdir, rf, ci)),
    "Fiberseq":   lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("Fiberseq.smk",  runFiberseq(cfg, sid, gp, indir, outdir, rf)),
    "LRtranscriptome": lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("LRtranscriptome.smk", runLRtranscriptome(cfg, sid, indir, outdir, rf)),
    "DR":         lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("DR.smk",         runDR(cfg, sid, sp, gp, indir, outdir, rf)),
    "Population_genomics": lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("Population_genomics.smk", runPopulation_genomics(cfg, sid, indir, outdir, rf)),
    "spatial_transcriptomics": lambda cfg, sid, sp, gp, indir, outdir, meta, rf, ci: ("spatial_transcriptomics.smk", runspatial_transcriptomics(cfg, sid, indir, outdir, rf)),
}




def touch_test_inputs(meta_path: str) -> int:
    """测试模式:为 meta 引用的输入路径创建占位空文件/目录。

    只为通过存在性检查(MetaUtil 会跳过缺文件的样本),不放任何内容;
    路径口径与 MetaUtil.load_meta 一致——相对路径以 meta 文件所在目录为基准。
    占位文件是运行时产物,不进 git(见 .gitignore)。
    """
    import csv as _csv
    base = os.path.dirname(os.path.abspath(meta_path))
    created = 0
    with open(meta_path, encoding="utf-8") as fh:
        head = fh.read(2048)
        fh.seek(0)
        sep = "\t" if head.count("\t") >= head.count(",") else ","
        for row in _csv.DictReader(fh, delimiter=sep):
            for v in row.values():
                if not isinstance(v, str) or "/" not in v:
                    continue
                p = v if os.path.isabs(v) else os.path.join(base, v)
                if os.path.exists(p):
                    continue
                if v.endswith("/") or "." not in os.path.basename(v.rstrip("/")):
                    os.makedirs(p, exist_ok=True)
                else:
                    os.makedirs(os.path.dirname(p), exist_ok=True)
                    open(p, "a").close()
                created += 1
    return created


def _cleanup_script_only_placeholders(outdir: str) -> None:
    """删除测试/--script-only 模式留下的产物占位(清单由 common.smk 的
    script_only_touch 维护于 log/cmd/.script_only_placeholders)。
    只删清单内的缺失占位,真实输出目录零残留;正常运行无清单,直接返回。"""
    manifest = os.path.join(outdir, "log", "cmd", ".script_only_placeholders")
    if not os.path.isfile(manifest):
        return
    removed = 0
    with open(manifest) as fh:
        paths = [line.strip() for line in fh if line.strip()]
    for p in paths:
        try:
            if os.path.isdir(p):
                os.rmdir(p)
            elif os.path.isfile(p):
                os.remove(p)
            removed += 1
        except OSError:
            pass
    os.remove(manifest)
    if removed:
        logger.info(f"script-only placeholders cleaned: {removed} ({outdir})")


def setup_test_args(args: argparse.Namespace, root_dir: str) -> argparse.Namespace:
    """Configure args for --test (dag|shell) mode.

    Workflow selection comes from -w (args.workflow_name); MODE picks the
    depth: dag = build the DAG only, shell = execute run-block code and write
    shell receipts without running tools.
    Modifies args in-place and returns it.
    """
    import shutil

    TEST_DIR = os.path.join(root_dir, "assests", "test")

    if args.workflow_name and "all" in args.workflow_name:
        args.workflow_name = list(WORKFLOW_DISPATCH)

    if args.test == "shell":
        args.dry_run = False
        os.environ["OMICS_SCRIPT_ONLY"] = "1"
        args.snakemake_args = list(args.snakemake_args or []) + ["--forceall", "--keep-going"]
    else:
        args.dry_run = True

    if args.output_dir:
        args.output_dir = os.path.join(args.output_dir, "test")
    else:
        args.output_dir = tempfile.mkdtemp(prefix="omics_test_")
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    logger.info(f"Output: {args.output_dir}")

    args.log = os.path.join(args.output_dir, "test.log")
    os.makedirs(args.output_dir, exist_ok=True)

    args._test_meta_map = {}
    for wf in args.workflow_name:
        meta_path = os.path.join(TEST_DIR, f"meta_{wf}.tsv")
        if os.path.isfile(meta_path):
            args._test_meta_map[wf] = meta_path
            continue
        variants = sorted(
            glob.glob(os.path.join(TEST_DIR, f"meta_{wf}_*.tsv")))
        if not variants:
            logger.info(f"Warning: no meta for {wf} at {meta_path}")
            continue
        if len(variants) == 1:
            args._test_meta_map[wf] = variants[0]
            continue
        args._test_meta_map[wf] = variants[0]
        for i, variant in enumerate(variants, start=1):
            args._test_meta_map[f"{wf}:{variant}"] = variant
        logger.info(
            f"{len(variants)} meta variants for {wf}: "
            + ", ".join(os.path.basename(v) for v in variants)
            + f" (default: {os.path.basename(variants[0])})")

    args.meta = None

    total_placeholders = 0
    for meta_path in set(args._test_meta_map.values()):
        if meta_path and os.path.isfile(meta_path):
            total_placeholders += touch_test_inputs(meta_path)
    if total_placeholders:
        logger.info(f"test input placeholders created: {total_placeholders}")

    test_data = os.path.join(TEST_DIR, "data")
    GENOME = "GRCm39"
    schema_validator = SchemaValidator()
    schema_validator._schema_dir = os.path.join(root_dir, "config")
    args._test_base_paths = schema_validator.generate_test_paths(test_data, GENOME)
    args._test_genome = GENOME

    args.conda_prefix = os.path.join(args.output_dir, ".conda")
    os.makedirs(args.conda_prefix, exist_ok=True)

    return args


def setup_normal_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate args for normal (non-test) mode."""
    if args.workflow_name and "all" in args.workflow_name:
        args.workflow_name = list(WORKFLOW_DISPATCH)
    if not args.meta:
        logger.info("Error: -m/--meta is required (unless --test is used)")
        exit(1)
    if args.design and not os.path.isfile(args.design):
        logger.info(f"Error: --design file not found: {args.design}")
        exit(1)
    if not args.output_dir:
        logger.info("Error: -o/--output_dir is required (unless --test is used)")
        exit(1)
    use_singularity = args.sdm is not None or _detect_singularity(args.snakemake_args)
    if not use_singularity and not args.conda_prefix:
        logger.info("Error: --conda-prefix is required when not using --sdm apptainer")
        exit(1)
    args._test_meta_map = None
    return args


def print_test_summary(test_results: Dict[str, Tuple[bool, str]]) -> None:
    """Print summary of test workflow results."""
    logger.info(f"[Results ({len(test_results)} workflows)")
    passed = [k for k, (ok, _) in test_results.items() if ok]
    failed = [k for k, (ok, _) in test_results.items() if not ok]
    for wf in sorted(test_results.keys()):
        ok, err = test_results[wf]
        status = "PASS" if ok else "FAIL"
        logger.info(f"  [{status}] {wf}")
        if not ok and err:
            logger.info(f"         {err.splitlines()[0] if err else 'unknown error'}")
    logger.info(f"\n  Passed: {len(passed)}, Failed: {len(failed)}, Total: {len(test_results)}")
    if failed:
        logger.error(f"  Failed workflows: {', '.join(failed)}")


_TMP_CACHE_DIRS = (".cache", ".config", ".keras", ".parallel", ".RepeatMaskerCache")
_TIDK_KEEP = os.path.join(".local", "share", "tidk", "tidk_database.csv")


_TMP_WORKFLOW_ENTRIES = ("star_align", "star_index", "bam_sort", "gatk")


def _run_tmp_dir(abs_outdir: str) -> Optional[str]:
    """The run's explicit tmp dir from its raw.json; None = 默认 /tmp 分流(设计 §9)。"""
    raw_path = os.path.join(abs_outdir, "raw.json")
    tmp_dir = None
    if os.path.isfile(raw_path):
        try:
            with open(raw_path, encoding="utf-8") as f:
                tmp_dir = json.load(f).get("tmp_dir")
        except (OSError, ValueError):
            tmp_dir = None
    if not tmp_dir:
        return None
    return tmp_dir if os.path.isabs(tmp_dir) else os.path.join(abs_outdir, tmp_dir)


def _run_tmp_paths(abs_outdir: str) -> List[str]:
    """Directories that may hold this run's temp files.

    Explicit tmp_dir → the run's own dir(整体清空);默认 /tmp 分流时只认
    工作流子目录(star_align、star_index*、bam_sort、gatk),绝不整扫 /tmp;
    兼容旧布局 <outdir>/tmp(存在才收)。
    """
    root = _run_tmp_dir(abs_outdir)
    if root and os.path.realpath(root) != "/tmp":
        return [root]
    paths = []
    legacy = os.path.join(abs_outdir, "tmp")
    if os.path.isdir(legacy):
        paths.append(legacy)
    for name in _TMP_WORKFLOW_ENTRIES:
        p = os.path.join("/tmp", name)
        if os.path.isdir(p):
            paths.append(p)
    try:
        for name in sorted(os.listdir("/tmp")):
            if name.startswith("star_index_") and os.path.isdir(os.path.join("/tmp", name)):
                paths.append(os.path.join("/tmp", name))
    except OSError:
        pass
    return paths


def _tree_bytes(path: str) -> int:
    total = 0
    if os.path.isfile(path):
        return 0 if os.path.islink(path) else os.path.getsize(path)
    for dp, _dns, fns in os.walk(path):
        for fn in fns:
            fp = os.path.join(dp, fn)
            if not os.path.islink(fp):
                try:
                    total += os.path.getsize(fp)
                except OSError:
                    pass
    return total


def report_temp_residue(abs_outdir: str, logger: logging.Logger) -> None:
    paths = [p for p in _run_tmp_paths(abs_outdir) if os.path.isdir(p)]
    n = sum(_tree_bytes(p) for p in paths)
    if n:
        logger.info("tmp residue in %s: %.1f MB (auto-cleaned after success; "
                    "use --clean-temp to clean now)", ", ".join(paths), n / 2 ** 20)


def clean_run_temp(abs_outdir: str, logger: logging.Logger,
                   include_caches: bool = False) -> None:
    """Clear a finished run's temp space.

    Always empties the run's temp dirs (explicit tmp_dir 整体清空;默认 /tmp
    分流时只清工作流子目录,目录本身保留)。With include_caches, also
    removes container/tool cache dirs at the output root, keeping the
    downloaded tidk database under .local.

    Args:
        abs_outdir: absolute workflow output dir (holds raw.json).
        logger: run logger.
        include_caches: also remove cache dirs (.cache/.config/.keras/...).
    """
    freed = 0
    for tmp_dir in _run_tmp_paths(abs_outdir):
        if not os.path.isdir(tmp_dir):
            continue
        freed += _tree_bytes(tmp_dir)
        for entry in os.listdir(tmp_dir):
            p = os.path.join(tmp_dir, entry)
            if os.path.islink(p) or os.path.isfile(p):
                os.remove(p)
            else:
                shutil.rmtree(p, ignore_errors=True)
        logger.info("cleaned tmp dir %s", tmp_dir)
    if freed:
        logger.info("cleaned tmp total: %.1f MB", freed / 2 ** 20)
    if include_caches:
        cfreed = 0
        for name in _TMP_CACHE_DIRS:
            p = os.path.join(abs_outdir, name)
            if os.path.isdir(p):
                cfreed += _tree_bytes(p)
                shutil.rmtree(p, ignore_errors=True)
        local = os.path.join(abs_outdir, ".local")
        keep = os.path.realpath(os.path.join(abs_outdir, _TIDK_KEEP))
        if os.path.isdir(local):
            for dp, _dns, fns in os.walk(local, topdown=False):
                for fn in fns:
                    fp = os.path.join(dp, fn)
                    if os.path.realpath(fp) != keep:
                        try:
                            cfreed += os.path.getsize(fp)
                            os.remove(fp)
                        except OSError:
                            pass
                for dn in _dns:
                    try:
                        os.rmdir(os.path.join(dp, dn))
                    except OSError:
                        pass
        logger.info("cleaned cache dirs under %s (%.1f MB)",
                    abs_outdir, cfreed / 2 ** 20)


def execute_workflows(args: argparse.Namespace, root_dir: str, logger: logging.Logger) -> None:
    """Execute all configured workflows.

    In test mode (args._test_meta_map is set), runs each workflow and collects results.
    In normal mode, runs workflows directly.
    """
    workflow_names: List[str] = args.workflow_name
    n_workflows = len(workflow_names)

    def _get_meta(wf_name: str) -> Optional[str]:
        if args._test_meta_map:
            return args._test_meta_map.get(wf_name, args.meta)
        return args.meta

    def _resolve_test_meta(meta_path: str, test_data_dir: str) -> str:
        """Resolve relative FASTQ paths in test meta file to absolute paths."""
        import pandas as pd
        import tempfile
        df = pd.read_csv(meta_path, sep="\t")
        path_cols = [c for c in df.columns if c.startswith("fastq_") or c in ("bam", "pbi", "ms_file")]
        changed = False
        for col in path_cols:
            for i, val in df[col].items():
                if pd.isna(val) or not isinstance(val, str):
                    continue
                val_str = str(val)
                if os.path.isabs(val_str) and os.path.exists(val_str):
                    continue
                candidate = os.path.join(test_data_dir, val_str)
                if os.path.isabs(val_str) and not os.path.exists(candidate):
                    for marker in ("fastq/", ".fq.gz", ".bam"):
                        idx = val_str.find(marker)
                        if idx > 0:
                            candidate = os.path.join(test_data_dir, val_str[idx:])
                            break
                if os.path.exists(candidate):
                    df.at[i, col] = candidate
                    changed = True
                elif not os.path.isabs(val_str):
                    logger.warning(f"[test] Path not found: {candidate}")
        if changed:
            tmp = meta_path + ".resolved"
            df.to_csv(tmp, sep="\t", index=False)
            return tmp
        return meta_path

    ref_outdir = os.path.join(args.output_dir, workflow_names[0])
    abs_ref_outdir = os.path.abspath(ref_outdir)
    run_log_dir = os.path.join(abs_ref_outdir, "log")
    os.makedirs(run_log_dir, exist_ok=True)
    run_log_path = os.path.join(
        run_log_dir,
        f"run_{time.strftime('%Y%m%d_%H%M%S')}_pid{os.getpid()}.log",
    )
    _run_log_fh = logging.FileHandler(run_log_path, encoding="utf-8")
    _run_log_fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"))
    logging.getLogger().addHandler(_run_log_fh)
    logger.info(f"Run log file: {run_log_path}")
    design = None
    if args.design:
        with open(args.design, "r", encoding="utf-8") as fh:
            design = json.load(fh)
        logger.info(f"design model loaded from {args.design} "
                    f"({sum(len((b or {}).get('comparisons') or {}) for b in (design or {}).values())} comparisons)")
    first_meta = _get_meta(workflow_names[0])
    if first_meta and os.path.isfile(first_meta):
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, meta=first_meta, workflow=workflow_names[0], design=design)
    else:
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, fastq_dir=first_meta, workflow=workflow_names[0], design=design)
    samples_info_dict, sample_pairs, group_pairs, raw_fastq_dir, raw_files, cellranger_input_dict = metadataUtil.run()

    threads_per_workflow = max(1, args.threads // n_workflows)
    if n_workflows > 1:
        logger.info(f"Parallel mode: {n_workflows} workflows, "
                    f"{args.threads} total threads -> {threads_per_workflow} per workflow")

    test_results: Dict[str, Tuple[bool, str]] = {}
    smk_cmds: list[tuple[list[str], str]] = []
    pending: List[Dict[str, Any]] = []
    for wf_name in workflow_names:
        abs_outdir = os.path.abspath(os.path.join(args.output_dir, wf_name))
        os.makedirs(abs_outdir, exist_ok=True)
        report_temp_residue(abs_outdir, logger)
        wf_design = design

        try:
            if args._test_meta_map:
                wf_meta = _get_meta(wf_name)
                if wf_meta and os.path.isfile(wf_meta):
                    test_data = os.path.join(root_dir, "assests", "test", "data")
                    wf_meta = _resolve_test_meta(wf_meta, test_data)
                    metadataUtil = MetadataUtils(outdir=abs_outdir, meta=wf_meta, workflow=wf_name, design=wf_design)
                else:
                    metadataUtil = MetadataUtils(outdir=abs_outdir, fastq_dir=wf_meta, workflow=wf_name, design=wf_design)
                samples_info_dict, sample_pairs, group_pairs, raw_fastq_dir, raw_files, cellranger_input_dict = metadataUtil.run()

            model_json = os.path.join(root_dir, f"config/{wf_name}.json")
            workflow_config = _load_model_json(model_json)

            if args.schema_validate and args.extra_args:
                sv = SchemaValidator()
                sv._schema_dir = os.path.join(root_dir, "config")
                try:
                    sv.load_workflow(wf_name)
                    casted, errors = sv.cast_extra_args(args.extra_args)
                    if errors:
                        for err in errors:
                            logger.error(f"[{wf_name}] 参数校验失败: {err}")
                        raise ValueError(
                            f"extra_args 校验失败，共 {len(errors)} 个错误"
                        )
                    args.extra_args = casted
                except FileNotFoundError:
                    pass

            if hasattr(args, '_test_base_paths'):
                test_data = os.path.join(os.path.join(root_dir, "assests", "test"), "data")
                base_paths = args._test_base_paths

                def _is_path(val: Any) -> bool:
                    if val is None: return True
                    if not isinstance(val, str): return False
                    if "/" in val: return True
                    return False

                def _make_test_path(key: str, test_data: str, genome: str) -> str:
                    from pathlib import Path
                    ref = Path(test_data) / "ref"
                    ref.mkdir(parents=True, exist_ok=True)
                    for name in ("smallrna", "rRNA", "access", "repeat", "decoy"):
                        if name in key:
                            f = ref / name
                            f.touch(exist_ok=True)
                            return str(f)
                    if key == "assembly_dir":
                        asm_root = Path(test_data) / "telomere" / "asm"
                        for _s in ("TestSample1", "TestSample2"):
                            _ad = asm_root / _s / "assembly"
                            _ad.mkdir(parents=True, exist_ok=True)
                            _fa = _ad / "asm.bp.p_ctg.fa"
                            if not _fa.exists():
                                _fa.write_text(
                                    f">{_s}_ctg1\n" + "TTAGGG" * 300
                                    + "ACGT" * 5000 + "CCCTAA" * 300 + "\n")
                        return str(asm_root)
                    if key.endswith(("_dir", "_index")):
                        d = Path(test_data) / "index" / key.replace("_dir", "").replace("_index", "")
                        d.mkdir(parents=True, exist_ok=True)
                        return str(d)
                    if key.endswith("_prefix"):
                        d = Path(test_data) / "index" / key.replace("_index_prefix", "").replace("_prefix", "")
                        d.mkdir(parents=True, exist_ok=True)
                        return str(d / genome)
                    g = ref / genome
                    g.touch(exist_ok=True)
                    return str(g)

                wf_schema: Dict[str, Any] = {}
                try:
                    sv_inject = SchemaValidator()
                    sv_inject._schema_dir = os.path.join(root_dir, "config")
                    sv_inject.load_workflow(wf_name)
                    wf_schema = sv_inject._schema
                except FileNotFoundError:
                    wf_schema = {}

                def _schema_node(dotted: str) -> Dict[str, Any]:
                    node: Any = wf_schema
                    for seg in dotted.split("."):
                        if not isinstance(node, dict):
                            return {}
                        props = node.get("properties")
                        if isinstance(props, dict) and seg in props:
                            node = props[seg]
                        elif seg in node:
                            node = node[seg]
                        else:
                            return {}
                    return node if isinstance(node, dict) else {}

                def _inject(cfg: Dict[str, Any], prefix: str, wf_extra: Dict[str, Any]) -> None:
                    for field, val in cfg.items():
                        dotted = f"{prefix}.{field}" if prefix else field
                        if isinstance(val, dict):
                            _inject(val, dotted, wf_extra)
                        elif _is_path(val):
                            declared = _schema_node(dotted).get("type")
                            decl = declared if isinstance(declared, list) else (
                                [declared] if declared is not None else [])
                            if declared is not None and not any(
                                    t in ("str", "null") for t in decl):
                                continue
                            wf_extra[dotted] = base_paths.get(dotted, _make_test_path(field, test_data, args._test_genome))

                wf_extra: Dict[str, Any] = dict(args.extra_args) if hasattr(args, 'extra_args') else {}
                _inject(workflow_config, "", wf_extra)

                flat_args = {k: v for k, v in wf_extra.items() if '.' not in k}
                workflow_config.update(flat_args)
                dot_args = parse_dot_args(wf_extra)
                for key_tuple, v in dot_args.items():
                    dict_set_by_path(workflow_config, list(key_tuple), v)
            else:
                flat_args = {k: v for k, v in args.extra_args.items() if '.' not in k}
                workflow_config.update(flat_args)
                dot_args = parse_dot_args(args.extra_args)
                for key_tuple, v in dot_args.items():
                    dict_set_by_path(workflow_config, list(key_tuple), v)

            if wf_design is not None:
                workflow_config.setdefault("sample_params", {})["design"] = wf_design
            elif getattr(metadataUtil, "column_design", None) and not (
                    workflow_config.get("sample_params", {}).get("design")):
                workflow_config.setdefault("sample_params", {})["design"] = \
                    metadataUtil.column_design

            if getattr(args, '_test_meta_map', None) and wf_name == "scRNAseq":
                _wf = workflow_config.setdefault("Params", {}).setdefault("workflow", {})
                _wf.setdefault("counters", ["scTE", "cellranger"])
                _wf.setdefault("aligner", "cellranger")

            if wf_name not in WORKFLOW_DISPATCH:
                raise ValueError(f"Unknown workflow name: {wf_name}")

            smk, input_json = WORKFLOW_DISPATCH[wf_name](
                deepcopy(workflow_config), samples_info_dict, sample_pairs, group_pairs,
                raw_fastq_dir, abs_outdir, _get_meta(wf_name), raw_files, cellranger_input_dict
            )

            forcerun_targets: Optional[List[str]] = None
            if args.forcerun:
                forcerun_targets = _expand_wildcard_rule_targets(
                    _prefix_forcerun_targets(args.forcerun, wf_name),
                    wf_name, root_dir, smk, input_json, logger,
                )
            forcerun_downstream_targets: Optional[List[str]] = None
            if args.forcerun_downstream:
                forcerun_downstream_targets = _prefix_forcerun_targets(
                    args.forcerun_downstream, wf_name)
            pending.append({
                "wf": wf_name,
                "smk": smk,
                "input_json": input_json,
                "outdir": abs_outdir,
                "forcerun": forcerun_targets,
                "forcerun_downstream": forcerun_downstream_targets,
                "input_bytes": _workflow_input_bytes(input_json),
            })

        except Exception as e:
            if args._test_meta_map:
                test_results[wf_name] = (False, str(e))
                logger.error(f"[{wf_name}] Config/build failed: {e}")
            else:
                raise

    run_started = time.time()
    weights = [p["input_bytes"] for p in pending]

    # I/O 带宽预算(设计 §8.4):同一物理盘上的工作流按同一权重口径切分该盘预算,
    # 解决 P4(多工作流各自无上限压同一块盘);--io-budget-mb-s 0 显式关闭限流。
    io_budget_arg = getattr(args, "io_budget_mb_s", None)
    if io_budget_arg is not None and io_budget_arg < 0:
        raise ValueError("--io-budget-mb-s must be >= 0")
    disk_shares: Dict[int, int] = {}
    _dev_groups: Dict[int, List[int]] = {}
    for _i, _p in enumerate(pending):
        try:
            _dev = os.stat(_p["outdir"]).st_dev
        except OSError:
            _dev = 0
        _dev_groups.setdefault(_dev, []).append(_i)
    for _dev, _idxs in _dev_groups.items():
        if io_budget_arg == 0:
            for _i in _idxs:
                disk_shares[_i] = 0
            continue
        _total = io_budget_arg if io_budget_arg else _detect_disk_budget_mb_s(
            pending[_idxs[0]]["outdir"])
        _g_weights = [pending[_i]["input_bytes"] for _i in _idxs]
        for _i, _s in zip(_idxs, _split_disk_budget(_total, _g_weights)):
            disk_shares[_i] = _s

    mem_shares = _split_mem_budget(args.mem_mb, weights)
    for idx, (p, share) in enumerate(zip(pending, mem_shares)):
        disk_share = disk_shares.get(idx, 0)
        _inject_resource_caps(p["input_json"], mem_mb=share, cpus=threads_per_workflow,
                              disk_budget_mb_s=disk_share or None)
        if args.schema_validate:
            try:
                _sv = SchemaValidator()
                _sv._schema_dir = os.path.join(root_dir, "config")
                _sv.load_workflow(p["wf"])
                with open(p["input_json"], "r", encoding="utf-8") as _fh:
                    _errs = _sv.validate(json.load(_fh))
                if _errs:
                    for _e in _errs:
                        logger.error(f"[{p['wf']}] [schema] {_e}")
                    raise ValueError(
                        f"{p['wf']}: merged config failed schema validation "
                        f"({len(_errs)} error(s); fix the config or pass "
                        f"--no-schema-validate to skip)")
            except FileNotFoundError:
                pass
        cmd = build_snakemake_cmd(
            root_dir, p["smk"], p["input_json"], threads_per_workflow,
            args.conda_prefix, args.rerun_trigger, args.dry_run,
            args.conda_frontend, args.snakemake_args,
            sdm=args.sdm, singularity_args=args.singularity_args,
            forcerun=p["forcerun"], unlock=args.unlock,
            touch=args.touch, list_target_rules=args.list_target_rules,
            list_rules=args.list_rules,
            forcerun_downstream=p["forcerun_downstream"],
            mem_mb=share,
            disk_budget_mb_s=disk_share or None,
            restart_times=args.restart_times,
        )
        logger.info("[%s] mem_mb share=%d, disk budget=%s MB/s (input_bytes=%s)",
                    p["wf"], share, disk_share or "off", p["input_bytes"])
        logger.info(f"[{p['wf']}] {cmd}")
        smk_cmds.append((cmd, p["outdir"]))

    def _summary_payload(status: Dict[str, str]) -> Dict[str, Any]:
        """run_summary.json payload(效果评判 §7.4)。"""
        return {
            "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "makespan_s": round(time.time() - run_started, 1),
            "host": getattr(args, "host_profile", {}) or {},
            "dry_run": bool(args.dry_run),
            "threads_total": args.threads,
            "threads_per_workflow": threads_per_workflow,
            "mem_budget_mb": args.mem_mb,
            "io_budget_mb_s": getattr(args, "io_budget_mb_s", None),
            "workflows": [
                {"name": p["wf"], "input_bytes": p["input_bytes"],
                 "mem_mb_share": s,
                 "disk_budget_mb_s_share": disk_shares.get(i, 0),
                 "status": status.get(p["wf"], "unknown")}
                for i, (p, s) in enumerate(zip(pending, mem_shares))
            ],
        }

    if args._test_meta_map:
        import re as _re

        def _missing_outputs_only(exc: Exception) -> bool:
            """True = snakemake 输出中的每个 Error in rule 块都是
            MissingOutputException(测试模式跳过执行的预期结果)。"""
            cause = getattr(exc, "__cause__", None)
            out = ((getattr(cause, "stdout", "") or "")
                   + (getattr(cause, "stderr", "") or ""))
            if "MissingOutputException" not in out:
                logger.info("[test-judge] no MissingOutputException in snakemake output")
                return False
            blocks = out.split("Error in rule")[1:] or []
            if not blocks:
                logger.info("[test-judge] no 'Error in rule' blocks in snakemake output")
                return False
            for blk in blocks:
                names = _re.findall(r"\b([A-Za-z_]+(?:Exception|Error))\b", blk)
                names = [n for n in names if n != "WorkflowError"]
                logger.info(f"[test-judge] block exceptions: {sorted(set(names))}")
                if any(n != "MissingOutputException" for n in names):
                    return False
            return True

        exec_status: Dict[str, str] = {}
        for cmd, cwd in smk_cmds:
            wf = os.path.basename(cwd)
            try:
                _run_cmd(cmd, cwd=cwd)
                test_results[wf] = (True, "")
                exec_status[wf] = "ok"
            except Exception as e:
                if _missing_outputs_only(e):
                    test_results[wf] = (True, "回执已生成;产物按测试模式约定不落盘(MissingOutputException 已豁免)")
                    exec_status[wf] = "ok"
                else:
                    test_results[wf] = (False, str(e)[:200])
                    exec_status[wf] = "failed"
            finally:
                _cleanup_script_only_placeholders(cwd)

        _write_run_summary(args.output_dir, _summary_payload(exec_status))
        print_test_summary(test_results)
        wf_failed = [k for k, (ok, _) in test_results.items() if not ok]
        if wf_failed:
            logger.error("[test] FAILED: workflows: %s", ", ".join(wf_failed))
            raise SystemExit(1)
    else:
        exec_status = {p["wf"]: "not-run" for p in pending}
        try:
            if n_workflows == 1:
                _run_cmd(smk_cmds[0][0], cwd=smk_cmds[0][1])
            else:
                logger.info(f"Launching {n_workflows} snakemake processes in parallel...")
                _run_cmds_parallel(smk_cmds)
                logger.info("All workflows completed.")
            for wf in exec_status:
                exec_status[wf] = "ok"
        except BaseException as exc:
            interrupted = isinstance(exc, (KeyboardInterrupt, SignalExit))
            for wf in exec_status:
                if exec_status[wf] == "not-run":
                    exec_status[wf] = "interrupted" if interrupted else "failed"
            if interrupted:
                logger.error("run interrupted (signal): writing run summary "
                             "with status 'interrupted' — NOT a normal completion")
            _write_run_summary(args.output_dir, _summary_payload(exec_status))
            raise
        finally:
            for _, cwd in smk_cmds:
                _cleanup_script_only_placeholders(cwd)
        _write_run_summary(args.output_dir, _summary_payload(exec_status))
        if all(s == "ok" for s in exec_status.values()):
            for _, cwd in smk_cmds:
                clean_run_temp(cwd, logger)


def _build_execution_callback(args: argparse.Namespace,
                              logger: logging.Logger) -> Optional[ExecutionCallback]:
    """Resolve the execution callback from CLI flags and environment.

    CLI flags win over environment variables; the callback arms only when
    BOTH the execution id and the token resolve (both required by the
    platform API). A partially configured setup logs a warning and stays
    disarmed.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments (callback-url / callback-execution-id /
        callback-token / no-callback).
    logger : logging.Logger

    Returns
    -------
    ExecutionCallback or None
        None when disabled (--no-callback) or not fully configured.
    """
    if getattr(args, "no_callback", False):
        logger.info("execution callback disabled by --no-callback")
        return None
    url = (args.callback_url or os.environ.get("BIPLATFORM_CALLBACK_URL")
           or DEFAULT_CALLBACK_URL)
    exec_id = (args.callback_execution_id
               or os.environ.get("BIPLATFORM_EXEC_ID"))
    token = args.callback_token or os.environ.get("BIPLATFORM_CALLBACK_TOKEN")
    if not exec_id or not token:
        if exec_id or token:
            logger.warning(
                "execution callback not armed: BOTH execution id and token "
                "are required (id: --callback-execution-id or env "
                "BIPLATFORM_EXEC_ID; token: "
                "--callback-token or env BIPLATFORM_CALLBACK_TOKEN)")
        return None
    cb = ExecutionCallback(url, exec_id, token, logger=logger)
    logger.info(f"execution callback armed: url={url} executionId={cb.execution_id}")
    return cb


if __name__ == "__main__":
    args = parse_args()

    callback = _build_execution_callback(args, logger)

    _signal_state = {"received": False}

    def _force_exit(signum):
        terminate_all_children(force=True)
        os._exit(128 + signum)

    def _on_terminate(signum, frame):
        if _signal_state["received"]:
            _force_exit(signum)
        _signal_state["received"] = True
        raise SignalExit(signum)

    def _on_sigint(signum, frame):
        if _signal_state["received"]:
            _force_exit(signum)
        _signal_state["received"] = True
        raise KeyboardInterrupt(f"interrupted by signal {signum}")

    signal.signal(signal.SIGTERM, _on_terminate)
    signal.signal(signal.SIGINT, _on_sigint)

    try:
        host_profile = detect_host_profile()
        args.host_profile = host_profile
        if args.threads is None:
            cpus = host_profile["cpus"]
            args.threads = max(1, min(int(cpus * 0.9), cpus - 1))
            logger.info("cores auto: %d (host %s, cpus=%d)",
                        args.threads, host_profile["hostname"], host_profile["cpus"])
        if args.mem_mb is None:
            args.mem_mb = host_profile["mem_mb"]
        ROOT_DIR = os.path.dirname(__file__)

        if args.touch and args.dry_run:
            logger.warning("--touch 和 --dry-run 同时指定，--dry-run 会阻止 --touch 生效，将忽略 --touch")
            args.touch = False

        if args.test is not None:
            setup_test_args(args, ROOT_DIR)
        else:
            setup_normal_args(args)

        logger = setup_logger("root", level=logging.DEBUG, log_file=args.log)
        if args.clean_temp:
            for wf in args.workflow_name:
                clean_run_temp(os.path.abspath(os.path.join(args.output_dir, wf)),
                               logger, include_caches=True)
            raise SystemExit(0)
        execute_workflows(args, ROOT_DIR, logger)
    except SystemExit as e:
        code = e.code
        if isinstance(code, str):
            if callback is not None:
                callback.send("FAILED", message=code, exit_code=1)
            raise
        code = 0 if code is None else int(code)
        if callback is not None:
            callback.send("SUCCESS" if code == 0 else "FAILED",
                          message="" if code == 0 else f"run.py exited with code {code}",
                          exit_code=code)
        raise
    except (KeyboardInterrupt, SignalExit) as exc:
        signum = getattr(exc, "signum", signal.SIGINT)
        sig_name = signal.Signals(signum).name
        exit_code = 128 + signum
        logger.error("=" * 70)
        logger.error("run.py INTERRUPTED by %s — this is NOT a normal "
                     "completion; exit code %d", sig_name, exit_code)
        logger.error("=" * 70)
        terminate_all_children(force=True)
        if callback is not None:
            callback.send("CANCELLED", message=f"interrupted ({sig_name})",
                          exit_code=exit_code)
        sys.exit(exit_code)
    except Exception as e:
        if callback is not None:
            callback.send("FAILED", message=f"{type(e).__name__}: {e}", exit_code=1)
        raise
    if callback is not None:
        callback.send("SUCCESS", exit_code=0)