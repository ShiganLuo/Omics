import argparse
from copy import deepcopy
import json
import os
import re
from src.common.util.MetaUtil import MetadataUtils
from src.common.util.LogUtil import setup_logger
from src.common.util.CmdUtil import _run_cmd, _run_cmds_parallel
from src.common.util.SchemaValidatorUtil import SchemaValidator
from src.common.util.EnvUtil import is_path_like
from node import runCoCulture, runMERIP, runRNAseq, runncRNAseq, runCLIP, runMutation, runPacVar, runKARRseq, runPeakCalling, runQuantMS, runtRNAseq
import logging
from typing import Dict, Any
logger = setup_logger(__name__, level=logging.DEBUG)

def smart_cast(val):
    """尝试将字符串转换为 int/float/bool，否则原样返回；list 逐元素转换"""
    if isinstance(val, list):
        return [smart_cast(v) for v in val]
    if isinstance(val, str):
        if val.lower() in {"true", "false"}:
            return val.lower() == "true"
        try:
            if val.startswith("0") and len(val) > 1 and not val.startswith("0."):
                return val  # 避免八进制等
            return int(val)
        except Exception:
            pass
        try:
            return float(val)
        except Exception:
            pass
    return val

def dict_set_by_path(d, keys, value):
    """递归设置嵌套字典的值，keys为key列表，自动类型转换"""
    for k in keys[:-1]:
        if k not in d or not isinstance(d[k], dict):
            d[k] = {}
        d = d[k]
    d[keys[-1]] = smart_cast(value)

def parse_dot_args(extra_args):
    """从extra_args中提取点号语法参数，返回{(k1,k2,...):v}"""
    dot_args = {}
    for k, v in list(extra_args.items()):
        if '.' in k:
            dot_args[tuple(k.split('.'))] = v
    return dot_args


def _load_model_json(model_json_file: str) -> Dict[str, Any]:
    """Load model JSON template from disk."""
    with open(model_json_file, 'r', encoding='utf-8') as f:
        return json.load(f)
    


def parse_args():
    parser = argparse.ArgumentParser(description="workflow")
    parser.add_argument('-m','--meta', type=str, default=None, help='meta input file or data dir which condatain fastq file')
    parser.add_argument('-w','--workflow_name', type=str, nargs='+',
        choices=["CoCulture", "MERIP", "RNAseq", "ncRNAseq", "CLIP", "Mutation", "PacVar", "KARRseq", "PeakCalling", "QuantMS", "tRNAseq"],
        default=['CoCulture'], help='workflow name(s), multiple for parallel execution')
    parser.add_argument('-o','--output_dir', type=str, default=None, help='output dir')
    parser.add_argument('-t','--threads', type=int, default=10, help='threads')
    parser.add_argument('--dry-run', action='store_true', help='dry run')
    parser.add_argument('--test', type=str, nargs='?', const='all', metavar='WORKFLOW',
        help='run dry-run test for a workflow (or "all"). Auto-sets meta/output/dry-run from test/ directory')
    parser.add_argument('--log', type=str, default='workflow.log', help='log file')
    parser.add_argument('--conda-prefix', type=str, default=None, help='conda prefix for snakemake (required when NOT using --sdm apptainer)')
    parser.add_argument('--sdm', action='store_const', const='apptainer', default=None,
        help='use apptainer container backend (SIF images). When set, --use-conda is omitted')
    parser.add_argument('--singularity-args', type=str, default=None,
        help='arguments for singularity/apptainer, e.g. --singularity-args \'--bind /path1,/path2\' (use with --sdm)')
    parser.add_argument(
        '--rerun-trigger', '--rerun-triggers',
        dest='rerun_trigger',
        nargs='+',
        default=["code", "input", "mtime", "params", "software-env"],
        choices=["code", "input", "mtime", "params", "software-env"],
        help='snakemake rerun-triggers, e.g. code input mtime params software-env'
    )
    parser.add_argument('--conda-frontend', type=str, choices=["conda", "mamba"], default="mamba", help='conda frontend for snakemake')
    parser.add_argument('--forcerun', type=str, nargs='+', default=None,
        help='force re-run specific jobs without downstream, e.g. --forcerun trimming_Paired:sample_id=S1')
    parser.add_argument(
        '--snakemake-args',
        nargs=argparse.REMAINDER,
        default=[],
        help='additional arguments forwarded to snakemake; place them after this flag'
    )
    
    # 支持 --key=value、--key value、--key v1 v2 v3 三种形式的额外参数
    # 多值参数在碰到下一个 --key 或到达末尾时停止收集
    args, unknown = parser.parse_known_args()
    extra_args = {}
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


def _detect_singularity(snakemake_args):
    """Check if --sdm apptainer/singularity is present in raw snakemake_args (backward compat)."""
    snakemake_args = snakemake_args or []
    return any(
        arg in ("apptainer", "singularity")
        for i, arg in enumerate(snakemake_args)
        if i > 0 and snakemake_args[i - 1] == "--sdm"
    )


def _collect_bind_paths(config_path):
    """Scan a config JSON file and collect directory paths to bind-mount.

    Recursively walks all string values; strings containing '/' are treated as
    paths. Files yield their parent directory; directories are kept as-is.
    Subdirectories are collapsed into their parents to minimise the bind list.
    """
    import json as _json
    import os as _os

    with open(config_path, "r", encoding="utf-8") as f:
        cfg = _json.load(f)

    dirs = set()

    def _walk(obj):
        if is_path_like(obj):
            # Resolve to absolute path; use as-is if not exists (don't resolve symlinks)
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

    # Collapse: remove any directory that is a child of another in the set
    result = set(dirs)
    for d in list(dirs):
        for other in dirs:
            if d != other and d.startswith(other.rstrip("/") + "/"):
                result.discard(d)
                break
    return sorted(result)

def _merge_singularity_args(
    json_bind_paths,
    singularity_args=None,
):
    """
    Merge bind paths collected from the JSON configuration with
    user-provided Singularity arguments.

    JSON-derived bind paths are always preserved. User-provided
    bind paths are added to them. Duplicate paths are removed while
    preserving their original order.

    Other Singularity arguments, such as --cleanenv, are preserved.
    """
    bind_paths = []

    def add_bind_paths(paths):
        for path in paths:
            path = path.strip()
            if path and path not in bind_paths:
                bind_paths.append(path)

    # Paths collected from config.
    add_bind_paths(json_bind_paths)

    # Always make /tmp available inside the container.
    add_bind_paths(["/tmp"])

    if not singularity_args:
        return (
            "--bind " + ",".join(bind_paths)
            if bind_paths
            else None
        )

    # Find an existing --bind argument.
    match = re.search(
        r"(?:^|\s)--bind(?:=|\s+)([^\s]+)",
        singularity_args,
    )

    if match:
        user_bind = match.group(1)

        # Add user-defined paths.
        add_bind_paths(user_bind.split(","))

        # Replace only the existing --bind argument.
        merged_bind = "--bind " + ",".join(bind_paths)

        start, end = match.span()
        singularity_args = (
            singularity_args[:start]
            + merged_bind
            + singularity_args[end:]
        )
    else:
        # User supplied other singularity arguments but no --bind.
        singularity_args = (
            singularity_args.rstrip()
            + " --bind "
            + ",".join(bind_paths)
        )

    return singularity_args.strip()

def build_snakemake_cmd(root_dir, smk, input_json, threads, conda_prefix, rerun_trigger,
                        dry_run, conda_frontend, snakemake_args, sdm=None, singularity_args=None,
                        forcerun=None):
    snakemake_args = snakemake_args or []
    # Determine container backend: explicit --sdm flag, or legacy --snakemake-args --sdm
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
    if not use_singularity:
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
        )

        if singularity_args:
            cmd += [
                "--singularity-args",
                singularity_args,
            ]
    if dry_run:
        cmd.append("--dry-run")
    if forcerun:
        cmd.append("--until")
        cmd.extend(forcerun)
        cmd.append("--forcerun")
        cmd.extend(forcerun)
    if snakemake_args:
        cmd.extend(snakemake_args)
    return cmd


WORKFLOW_DISPATCH = {
    "CoCulture":  lambda cfg, sid, sp, gp, indir, outdir, meta: ("CoCulture.smk", runCoCulture(cfg, sid, indir, outdir)),
    "MERIP":      lambda cfg, sid, sp, gp, indir, outdir, meta: ("MERIP.smk",     runMERIP(cfg, sid, indir, outdir)),
    "RNAseq":     lambda cfg, sid, sp, gp, indir, outdir, meta: ("RNAseq.smk",    runRNAseq(cfg, sid, gp, indir, outdir)),
    "ncRNAseq":   lambda cfg, sid, sp, gp, indir, outdir, meta: ("ncRNAseq.smk",  runncRNAseq(cfg, sid, indir, outdir)),
    "CLIP":       lambda cfg, sid, sp, gp, indir, outdir, meta: ("CLIP.smk",      runCLIP(cfg, sid, indir, outdir)),
    "Mutation":   lambda cfg, sid, sp, gp, indir, outdir, meta: ("Mutation.smk",  runMutation(cfg, sid, sp, indir, outdir)),
    "PacVar":     lambda cfg, sid, sp, gp, indir, outdir, meta: ("PacVar.smk",    runPacVar(cfg, sid, indir, outdir)),
    "KARRseq":    lambda cfg, sid, sp, gp, indir, outdir, meta: ("KARRseq.smk",   runKARRseq(cfg, sid, indir, outdir)),
    "PeakCalling":lambda cfg, sid, sp, gp, indir, outdir, meta: ("PeakCalling.smk",runPeakCalling(cfg, sid,sp, indir, outdir)),
    "QuantMS":    lambda cfg, sid, sp, gp, indir, outdir, meta: ("QuantMS.smk",   runQuantMS(cfg, sid, indir, outdir)),
    "tRNAseq":    lambda cfg, sid, sp, gp, indir, outdir, meta: ("tRNAseq.smk",   runtRNAseq(cfg, sid, indir, outdir, meta)),
}


# ============================================================
# Test path generation
# ============================================================


def setup_test_args(args, root_dir: str):
    """Configure args for --test mode.

    Resolves workflow names, output directory, meta files, and test paths.
    Modifies args in-place and returns it.
    """
    import shutil

    TEST_DIR = os.path.join(root_dir, "assests", "test")

    # All registered workflows
    ALL_WORKFLOWS = list(WORKFLOW_DISPATCH.keys())
    if args.test == "all":
        args.workflow_name = ALL_WORKFLOWS
    elif args.test in ALL_WORKFLOWS:
        args.workflow_name = [args.test]
    else:
        logger.info(f"Unknown workflow: {args.test}")
        logger.info(f" Available: {ALL_WORKFLOWS} or 'all'")
        exit(1)

    # Output to {cwd or --output-dir}/test
    base_out = args.output_dir if args.output_dir else os.getcwd()
    args.output_dir = os.path.join(base_out, "test")
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    logger.info(f"Output: {args.output_dir}")

    args.dry_run = True
    args.log = os.path.join(args.output_dir, "test.log")
    os.makedirs(args.output_dir, exist_ok=True)

    # Build per-workflow meta map dynamically
    args._test_meta_map = {}
    for wf in args.workflow_name:
        meta_path = os.path.join(TEST_DIR, f"meta_{wf}.tsv")
        if os.path.isfile(meta_path):
            args._test_meta_map[wf] = meta_path
        else:
            logger.info(f"Warning: no meta for {wf} at {meta_path}")

    args.meta = None  # will be resolved per-workflow via _get_meta()

    # Store schema validator and test data for per-workflow path injection
    test_data = os.path.join(TEST_DIR, "data")
    GENOME = "GRCm39"
    schema_validator = SchemaValidator()
    schema_validator._schema_dir = os.path.join(root_dir, "config")
    args._test_base_paths = schema_validator.generate_test_paths(test_data, GENOME)
    args._test_genome = GENOME

    # Local conda-prefix (avoid permission issues)
    args.conda_prefix = os.path.join(args.output_dir, ".conda")
    os.makedirs(args.conda_prefix, exist_ok=True)

    return args


def setup_normal_args(args):
    """Validate args for normal (non-test) mode."""
    if not args.meta:
        logger.info("Error: -m/--meta is required (unless --test is used)")
        exit(1)
    if not args.output_dir:
        logger.info("Error: -o/--output_dir is required (unless --test is used)")
        exit(1)
    # Determine container backend: explicit --sdm flag, or legacy --snakemake-args --sdm
    use_singularity = args.sdm is not None or _detect_singularity(args.snakemake_args)
    if not use_singularity and not args.conda_prefix:
        logger.info("Error: --conda-prefix is required when not using --sdm apptainer")
        exit(1)
    args._test_meta_map = None
    return args


def print_test_summary(test_results: Dict[str, tuple]):
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
        exit(1)


def execute_workflows(args, root_dir: str, logger):
    """Execute all configured workflows.

    In test mode (args._test_meta_map is set), runs each workflow and collects results.
    In normal mode, runs workflows directly.
    """
    workflow_names = args.workflow_name
    n_workflows = len(workflow_names)

    def _get_meta(wf_name):
        if args._test_meta_map:
            return args._test_meta_map.get(wf_name, args.meta)
        return args.meta

    # Use first workflow's output dir for metadata (or a shared parent if multi)
    ref_outdir = os.path.join(args.output_dir, workflow_names[0])
    abs_ref_outdir = os.path.abspath(ref_outdir)
    first_meta = _get_meta(workflow_names[0])
    if first_meta and os.path.isfile(first_meta):
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, meta=first_meta)
    else:
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, fastq_dir=first_meta)
    samples_info_dict, sample_pairs, group_pairs, raw_fastq_dir = metadataUtil.run()

    # Thread allocation: user-specified total threads split across workflows
    threads_per_workflow = max(1, args.threads // n_workflows)
    if n_workflows > 1:
        logger.info(f"Parallel mode: {n_workflows} workflows, "
                    f"{args.threads} total threads -> {threads_per_workflow} per workflow")

    # Prepare each workflow
    test_results = {}  # wf_name -> (passed: bool, error: str)
    smk_cmds: list[tuple[list[str], str]] = []  # (cmd, cwd) pairs
    for wf_name in workflow_names:
        abs_outdir = os.path.abspath(os.path.join(args.output_dir, wf_name))
        os.makedirs(abs_outdir, exist_ok=True)

        try:
            # In test mode, reload metadata per workflow (different meta files)
            if args._test_meta_map:
                wf_meta = _get_meta(wf_name)
                if wf_meta and os.path.isfile(wf_meta):
                    metadataUtil = MetadataUtils(outdir=abs_outdir, meta=wf_meta)
                else:
                    metadataUtil = MetadataUtils(outdir=abs_outdir, fastq_dir=wf_meta)
                samples_info_dict, sample_pairs, group_pairs, raw_fastq_dir = metadataUtil.run()

            model_json = os.path.join(root_dir, f"config/{wf_name}.json")
            workflow_config = _load_model_json(model_json)

            # In test mode, inject test paths for ALL path-like fields in config
            if hasattr(args, '_test_base_paths'):
                test_data = os.path.join(os.path.join(root_dir, "assests", "test"), "data")
                base_paths = args._test_base_paths

                def _is_path(val):
                    if val is None: return True
                    if not isinstance(val, str): return False
                    if "/" in val: return True
                    return False

                def _make_test_path(key, test_data, genome):
                    from pathlib import Path
                    ref = Path(test_data) / "ref"
                    ref.mkdir(parents=True, exist_ok=True)
                    for name in ("smallrna", "rRNA", "access", "repeat", "decoy"):
                        if name in key: return str(ref / name)
                    if key.endswith(("_dir", "_index")):
                        d = Path(test_data) / "index" / key.replace("_dir", "").replace("_index", "")
                        d.mkdir(parents=True, exist_ok=True)
                        return str(d)
                    if key.endswith("_prefix"):
                        d = Path(test_data) / "index" / key.replace("_index_prefix", "").replace("_prefix", "")
                        d.mkdir(parents=True, exist_ok=True)
                        return str(d / genome)
                    return str(ref / genome)

                def _inject(cfg, prefix, wf_extra):
                    for field, val in cfg.items():
                        dotted = f"{prefix}.{field}" if prefix else field
                        if isinstance(val, dict):
                            _inject(val, dotted, wf_extra)
                        elif _is_path(val):
                            wf_extra[dotted] = base_paths.get(dotted, _make_test_path(field, test_data, args._test_genome))

                wf_extra = dict(args.extra_args) if hasattr(args, 'extra_args') else {}
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

            if wf_name not in WORKFLOW_DISPATCH:
                raise ValueError(f"Unknown workflow name: {wf_name}")

            smk, input_json = WORKFLOW_DISPATCH[wf_name](
                deepcopy(workflow_config), samples_info_dict, sample_pairs, group_pairs,
                raw_fastq_dir, abs_outdir, _get_meta(wf_name)
            )

            # Auto-prefix forcerun targets with workflow name if not already prefixed.
            # In subworkflows, rules are renamed via "use rule ... as <wf>_...", so
            # the user can write "function_gsea" instead of "RNAseq_function_gsea".
            # For wildcards targets like "trimming_Paired:sample_id=S1", only the
            # rule name (before ":") gets prefixed.
            forcerun_targets = None
            if args.forcerun:
                forcerun_targets = []
                for t in args.forcerun:
                    if ":" in t:
                        rule_part, rest = t.split(":", 1)
                        rest = ":" + rest
                    else:
                        rule_part, rest = t, ""
                    if rule_part == "all" or rule_part.startswith(f"{wf_name}_"):
                        forcerun_targets.append(t)
                    else:
                        forcerun_targets.append(f"{wf_name}_{rule_part}{rest}")

            cmd = build_snakemake_cmd(
                root_dir, smk, input_json, threads_per_workflow,
                args.conda_prefix, args.rerun_trigger, args.dry_run,
                args.conda_frontend, args.snakemake_args,
                sdm=args.sdm, singularity_args=args.singularity_args,
                forcerun=forcerun_targets,
            )
            logger.info(f"[{wf_name}] {cmd}")
            smk_cmds.append((cmd, abs_outdir))

        except Exception as e:
            if args._test_meta_map:
                test_results[wf_name] = (False, str(e))
                logger.error(f"[{wf_name}] Config/build failed: {e}")
            else:
                raise

    # Execute snakemake commands
    if args._test_meta_map:
        # Test mode: run each workflow, catch errors
        for cmd, cwd in smk_cmds:
            wf = os.path.basename(cwd)
            try:
                _run_cmd(cmd, cwd=cwd)
                test_results[wf] = (True, "")
            except Exception as e:
                test_results[wf] = (False, str(e)[:200])

        print_test_summary(test_results)
    else:
        if n_workflows == 1:
            _run_cmd(smk_cmds[0][0], cwd=smk_cmds[0][1])
        else:
            logger.info(f"Launching {n_workflows} snakemake processes in parallel...")
            _run_cmds_parallel(smk_cmds)
            logger.info("All workflows completed.")


if __name__ == "__main__":
    args = parse_args()
    ROOT_DIR = os.path.dirname(__file__)

    # Configure args based on mode
    if args.test is not None:
        setup_test_args(args, ROOT_DIR)
    else:
        setup_normal_args(args)

    logger = setup_logger("root", level=logging.DEBUG, log_file=args.log)
    execute_workflows(args, ROOT_DIR, logger)