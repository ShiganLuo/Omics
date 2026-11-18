"""Common utilities for all modules.

This file is included by modules to provide shared utilities.
Usage:
    include: "../common/common.smk"

Then you can use:
    - setup_logger: Create and configure a logger
    - _require: Validate a required config input, naming the JSON key to fill
    - time, shutil, os, sys: Standard library modules
    - ROOT_DIR: Project root directory from config
    - sif: Resolve SIF container path for a conda env YAML
"""

import sys
import os
import time
import shutil
import tempfile
from snakemake.logging import logger
import subprocess
import shlex
import re
# Get ROOT_DIR from config (set by run.py)
ROOT_DIR = config.get("ROOT_DIR", ".")

# Output-root constant for command-receipt backups: subworkflows pass
# OUTDIR_ROOT (the workflow's total output directory, absolute). Receipt
# backups live flat under OUTDIR_ROOT/log/cmd/ — file names already carry
# rule + sample + timestamp, so no per-key subdirectories. The directory
# is created here once; rules only build their file path.
OUTDIR_ROOT = config.get("OUTDIR_ROOT") or "output"
RECEIPT_DIR = os.path.join(OUTDIR_ROOT, "log", "cmd")
os.makedirs(RECEIPT_DIR, exist_ok=True)

# 测试模式(run.py --test / --script-only 置环境变量 OMICS_SCRIPT_ONLY=1):
# 输出 shell 脚本但不执行——run 块末尾以 SCRIPT_ONLY 判断跳过 bash 执行,回执
# 脚本照常双写落盘,供跨版本 diff 对照参数传递语义。
SCRIPT_ONLY = os.environ.get("OMICS_SCRIPT_ONLY", "") == "1"


def script_only_touch(out):
    """测试模式:为产物写占位文件,保证 snakemake job 能"成功"走完整个 DAG。

    不 touch 占位的话,每个 job 都以 MissingOutputException"失败",下游被
    全部跳过,回执只能覆盖 DAG 根规则(实测 4/55)。约束:只创建缺失的产物
    (绝不 touch 已存在的真实产物,不污染其 mtime);创建的占位登记进
    <OUTDIR_ROOT>/log/cmd/.script_only_placeholders 清单,run.py 在 snakemake
    结束后统一删除,真实输出目录零残留。
    """
    if not SCRIPT_ONLY:
        return
    manifest = os.path.join(RECEIPT_DIR, ".script_only_placeholders")
    created = []
    items = [out] if isinstance(out, str) else list(out)
    for f in items:
        f = str(f)
        if f.endswith("/") or "." not in os.path.basename(f.rstrip("/")):
            if not os.path.exists(f):
                os.makedirs(f, exist_ok=True)
                created.append(f)
        elif not os.path.exists(f):
            d = os.path.dirname(f)
            if d:
                os.makedirs(d, exist_ok=True)
            open(f, "a").close()
            created.append(f)
    if created:
        with open(manifest, "a") as fh:
            fh.write("\n".join(created) + "\n")

# Ensure src directory is in sys.path for importing common modules
_src_dir = os.path.join(ROOT_DIR, "src")
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

# Import common utilities
try:
    from common.util.LogUtil import setup_logger
except ImportError as e:
    raise ImportError(
        f"Failed to import common.LogUtil.setup_logger. "
        f"Please ensure:\n"
        f"1. ROOT_DIR is set correctly in config (current: {ROOT_DIR})\n"
        f"2. Directory exists: {_src_dir}\n"
        f"3. File exists: {os.path.join(_src_dir, 'common', 'LogUtil.py')}\n"
        f"Original error: {e}"
    )

# ---------------------------------------------------------------------------
# Required-input helper
# ---------------------------------------------------------------------------

def _require(value, key, tool=""):
    """Return a configured input path or fail naming the JSON key to fill.

    ``key`` is the dotted config JSON key to set (e.g. ``genome.gtf``);
    ``tool`` optionally tags the error message with the module name.
    """
    if not value:
        who = f" for {tool}" if tool else ""
        raise ValueError(
            f"Missing reference file{who}: set '{key}' in the "
            f"workflow config JSON, then rerun."
        )
    return value


# ---------------------------------------------------------------------------
# Per-rule resource helper (config["resources"])
# ---------------------------------------------------------------------------

def _clamp_res_value(key, value, caps):
    """Clamp one resource value against the host budget (方案 A clamp).

    Only ``mem_mb`` and ``threads`` are clamped; every other key, and every
    value that fits the budget, passes through unchanged — rule semantics
    are preserved. ``caps`` comes from ``config["resource_caps"]`` injected
    by run.py (the per-process budget share); empty caps disable clamping
    (standalone run_module.py runs keep today's behaviour).
    """
    if not caps or isinstance(value, bool) or not isinstance(value, (int, float)):
        return value
    if key == "mem_mb" and caps.get("mem_mb"):
        cap = int(caps["mem_mb"] * 0.8)
    elif key == "threads" and caps.get("cpus"):
        cap = int(caps["cpus"])
    else:
        return value
    if value > cap:
        logger.warning(
            f"rule_res: {key}={value} exceeds host budget "
            f"(resource_caps={caps}), clamped to {cap}")
        return cap
    return value


def _resource_caps_from_env():
    """resource_caps fallback from the OMICS_RESOURCE_CAPS environment var.

    Snakemake module configs are hand-assembled per module and carry no
    resource_caps; run.py's CmdUtil exports OMICS_RESOURCE_CAPS from the
    workflow configfile so module rules clamp against the same per-process
    budget. Absent or invalid -> no clamping (run_module.py behaviour).
    """
    raw = os.environ.get("OMICS_RESOURCE_CAPS")
    if not raw:
        return {}
    try:
        import json
        caps = json.loads(raw)
        return caps if isinstance(caps, dict) else {}
    except ValueError:
        return {}


def input_size_mb(input_files):
    """Sum of input file sizes in MB (方案 B 规模度量,已定:字节数).

    ``input_files`` is a snakemake InputFiles instance (or any iterable of
    paths). Returns None when no file exists / total is 0 so callers fall
    back to the static value — never 0 (0 would crush mem_mb to min).
    """
    total = 0
    try:
        for f in input_files or []:
            try:
                total += os.path.getsize(str(f))
            except (OSError, TypeError):
                pass
    except TypeError:
        return None
    return total / 1e6 if total > 0 else None


def input_size_named(input_files, key):
    """Size in MB of the named input *key* when present, else of all inputs.

    Rules whose input block is positional (no named keys) would raise
    AttributeError on input.<key>; fall back to the whole input list so the
    byte-sum metric still applies.
    """
    try:
        sub = getattr(input_files, key)
    except AttributeError:
        sub = input_files
    return input_size_mb(sub)


_VALID_METRICS = ("fastq_mb", "bam_mb", "sample_count")


def sample_count():
    """Number of samples in the run (B.2 sample_count metric, for aggregate rules).

    Counts distinct sample ids across the metadata forms set by run.py
    (genome_paired_samples / genome_single_samples / samples). Returns None
    when no sample metadata is present so callers fall back to static values.
    """
    ids = set()
    for key in ("genome_paired_samples", "genome_single_samples"):
        for _, entries in (config.get(key) or {}).items():
            try:
                for pair in entries:
                    for s in (pair if isinstance(pair, (list, tuple)) else [pair]):
                        if s:
                            ids.add(str(s))
            except TypeError:
                pass
    for s in config.get("samples") or []:
        ids.add(str(s))
    return len(ids) if ids else None


def _validate_scale(name, scale):
    """Reject invalid scale blocks early (dry-run parse time, B.6)."""
    if not isinstance(scale, dict):
        raise ValueError(f"resources.{name}.scale must be a dict, got {type(scale).__name__}")
    missing = [k for k in ("size_metric", "ref_input_size", "fixed_frac",
                           "min_mem_mb", "max_mem_mb") if k not in scale]
    if missing:
        raise ValueError(f"resources.{name}.scale missing required keys: {missing}")
    if scale["size_metric"] not in _VALID_METRICS:
        raise ValueError(f"resources.{name}.scale.size_metric must be one of {_VALID_METRICS}")
    if not (0 <= float(scale["fixed_frac"]) <= 1):
        raise ValueError(f"resources.{name}.scale.fixed_frac must be in [0,1], got {scale['fixed_frac']}")
    if not (int(scale["min_mem_mb"]) <= int(scale["max_mem_mb"])):
        raise ValueError(f"resources.{name}.scale: min_mem_mb > max_mem_mb")
    if float(scale["ref_input_size"]) <= 0:
        raise ValueError(f"resources.{name}.scale.ref_input_size must be > 0")


def _scaled_mem_mb(name, merged, input_size):
    """方案 B 两段式缩放:固定开销 + 数据相关(B.3)。

    mem = base × (fixed_frac + (1 − fixed_frac) × input_size / ref_input_size)
    """
    scale = merged.get("scale")
    if not scale or input_size is None:
        return None  # caller keeps the static value
    _validate_scale(name, scale)
    base = float(merged.get("mem_mb", 0))
    if base <= 0:
        return None
    ref = float(scale["ref_input_size"])
    frac = float(scale["fixed_frac"])
    val = base * (frac + (1.0 - frac) * float(input_size) / ref)
    return int(min(max(val, int(scale["min_mem_mb"])), int(scale["max_mem_mb"])))


class _RuleRes(dict):
    """dict with host-clamped mem_mb/threads on access.

    Clamping happens in get()/__getitem__ rather than only at merge time so
    call-site fallbacks — rule_res(x).get("mem_mb", 32000) — cannot bypass
    the clamp and revive the P3 "job larger than budget" deadlock.
    """

    def __init__(self, caps):
        super().__init__()
        self._caps = caps

    def get(self, key, default=None):
        return _clamp_res_value(key, super().get(key, default), self._caps)

    def __getitem__(self, key):
        return _clamp_res_value(key, super().__getitem__(key), self._caps)


def rule_res(name, input_size=None, attempt=1):
    """Resolve resources for rule *name* from config["resources"].

    Four steps, order fixed (方案 B; clamp MUST be last so scaled/attempt
    values cannot bypass the host budget):
      1) static merge: resources.default < resources[name];
      2) scale: when a ``scale`` block is configured AND input_size is given,
         mem_mb follows the fixed-cost + data term formula (B.3);
      3) attempt escalation: attempt>1 doubles mem_mb (方案 C, 2**(attempt-1));
      4) host clamp via config["resource_caps"] / OMICS_RESOURCE_CAPS.

    Callers that pass no input_size skip step 2 and behave exactly like
    today — the 279 existing call sites need no change.

    Usage in a rule:
        threads: rule_res("fiberseq_track").get("threads", 4)
        resources:
            mem_mb = lambda wildcards, input, attempt: rule_res(
                "fiberseq_track", input_size=input_size_mb(input),
                attempt=attempt).get("mem_mb", 65000),
        # scaled pilot rules additionally pass input_size=input_size_mb(input)
    """
    res_cfg = config.get("resources") or {}
    merged = _RuleRes(config.get("resource_caps") or _resource_caps_from_env())
    merged.update(res_cfg.get("default") or {})
    merged.update(res_cfg.get(name) or {})

    scaled = _scaled_mem_mb(name, merged, input_size)
    if scaled is not None:
        merged["mem_mb"] = scaled
    if attempt and int(attempt) > 1 and merged.get("mem_mb"):
        merged["mem_mb"] = int(merged["mem_mb"]) * (2 ** (int(attempt) - 1))
    return merged


# ---------------------------------------------------------------------------
# Isobaric label sanity check
# ---------------------------------------------------------------------------
# TMT6/10/11/16/18 and iTRAQ 4/8 reporter m/z values (monoisotopic).
_REPORTER_MZ = [
    114.1276, 115.1247, 116.1217, 117.1189, 118.1345, 119.1316, 121.1344,
    126.1277, 127.1248, 127.1311, 128.1281, 128.1344, 129.1312, 129.1375,
    130.1345, 131.1379, 132.1412, 133.1382, 134.1416, 134.1446, 135.1416,
]


def check_reporter_ions(mzml_path, max_spectra=400):
    """Scan an mzML for isobaric reporter ions (sampled across the whole file).

    Returns ``(has_reporters, detail)``. A file counts as labelled when at
    least 3 sampled MS2 spectra carry >= 3 peaks within 20 ppm of known
    TMT/iTRAQ reporter m/z values. Spectra are reservoir-sampled over the
    ENTIRE file: early-gradient MS2 are often contaminant spectra without
    reporters even in labelled runs, so sampling only the head of the file
    yields false 'unlabelled' verdicts. Arrays that cannot be decoded with
    the standard library (e.g. numpress) are skipped rather than failing.

    Parameters
    ----------
    mzml_path : str
        mzML file to scan (indexed or plain).
    max_spectra : int
        Size of the reservoir sample of MS2 spectra.

    Returns
    -------
    tuple
        (bool, str) — whether reporter ions were found, plus a detail line.
    """
    import base64
    import random
    import struct
    import zlib
    import xml.etree.ElementTree as ET

    rng = random.Random(1)
    reservoir = []          # sampled (mz_raw, is_64bit) arrays
    total_ms2 = 0
    for _, el in ET.iterparse(mzml_path, events=("end",)):
        if el.tag.endswith("spectrum"):
            mslevel = None
            mz_raw = b""
            mz_is_64bit = True
            for node in el.iter():
                if node.tag.endswith("cvParam"):
                    acc = node.get("accession")
                    if acc == "MS:1000511":
                        mslevel = node.get("value")
            for bda in el.iter():
                if not bda.tag.endswith("binaryDataArray"):
                    continue
                is_mz = False
                is_64bit = True
                raw = b""
                for sub in bda:
                    if sub.tag.endswith("cvParam"):
                        acc = sub.get("accession")
                        if acc == "MS:1000514":
                            is_mz = True
                        elif acc == "MS:1000521":
                            is_64bit = False
                    elif sub.tag.endswith("binary"):
                        txt = sub.text or ""
                        try:
                            raw = zlib.decompress(base64.b64decode(txt))
                        except Exception:
                            try:
                                raw = base64.b64decode(txt)
                            except Exception:
                                raw = b""
                if is_mz:
                    mz_raw = raw
                    mz_is_64bit = is_64bit
            el.clear()
            try:
                is_ms2 = mslevel is not None and int(float(mslevel)) == 2
            except (TypeError, ValueError):
                is_ms2 = False
            if not is_ms2 or not mz_raw:
                continue
            total_ms2 += 1
            item = (mz_raw, mz_is_64bit)
            if len(reservoir) < max_spectra:
                reservoir.append(item)
            else:
                j = rng.randrange(total_ms2)
                if j < max_spectra:
                    reservoir[j] = item

    checked = 0
    hit_spectra = 0
    total_hits = 0
    for mz_raw, is_64bit in reservoir:
        checked += 1
        try:
            fmt = "<%dd" if is_64bit else "<%df"
            size = 8 if is_64bit else 4
            n = len(mz_raw) // size
            mzs = struct.unpack(fmt % n, mz_raw[:n * size])
        except Exception:
            continue
        hits = 0
        for rm in _REPORTER_MZ:
            tol = rm * 20e-6
            for mz in mzs:
                if abs(mz - rm) <= tol:
                    hits += 1
                    break
        total_hits += hits
        if hits >= 3:
            hit_spectra += 1
    detail = (
        f"{hit_spectra}/{checked} sampled MS2 spectra carry reporter-mass "
        f"peaks ({total_hits} reporter-mass peaks in total, {total_ms2} MS2 in file)"
    )
    # Labelled data shows reporter peaks in most peptide spectra; unlabelled
    # data produces only random m/z coincidences (~0.1 hits/spectrum, rarely
    # more than one >=3-hit spectrum per 400). Require at least 3 such
    # spectra to call the file labelled.
    return hit_spectra >= 3, detail


# ---------------------------------------------------------------------------
# Container (SIF) path resolution
# ---------------------------------------------------------------------------
_ENV_MAP = config.get("env", {})


def sif(yaml_filename: str, rule_name: str = "") -> str:
    """Resolve the SIF container path for a conda env YAML.

    Priority (specific over general, explicit config over filesystem):

    1. ``config["env"][<rule_name>]`` — explicit per-rule mapping
    2. ``config["env"][<yaml_stem>]`` — explicit module-level mapping,
       where ``yaml_stem`` is the YAML filename without extension
       (e.g. ``star.yaml`` -> ``star``)
    3. ``<env_dir>/<rule_name>/<yaml_stem>.sif`` — per-rule container
       directory, used only when the file exists on disk
    4. ``<env_dir>/<module_dir>/<yaml_stem>.sif`` — module-level
       fallback, where ``module_dir`` is ``workflow.basedir``'s basename
       and ``env_dir`` comes from ``config["env"]["env_dir"]``

    Parameters
    ----------
    yaml_filename : str
        The conda YAML filename (same value passed to ``conda:`` in a rule).
    rule_name : str, optional
        Name of the enclosing rule. Enables the per-rule container lookup
        (steps 1 and 3). Empty string disables per-rule resolution.

    Returns
    -------
    str
        Path to the ``.sif`` file.

    Notes
    -----
    Call sites pass the rule name explicitly, e.g.::

        container: sif("star.yaml", "star_align")
    """
    stem = os.path.splitext(os.path.basename(yaml_filename))[0]
    logger.debug(f"stem:{stem}, yaml_filename: {yaml_filename}, rule_name: {rule_name}")
    logger.debug(f"_ENV_MAP: {_ENV_MAP}")

    # 1. Explicit per-rule mapping in config["env"]
    if rule_name and rule_name in _ENV_MAP:
        return _ENV_MAP[rule_name]

    # 2. Explicit module-level mapping in config["env"]
    if stem in _ENV_MAP:
        return _ENV_MAP[stem]

    env_dir = _ENV_MAP.get("env_dir")
    if not env_dir:
        raise ValueError(
            f"SIF path for '{stem}' (rule '{rule_name}') not found in "
            f"config['env'] and no 'env_dir' fallback set. Please add "
            f"'\"{rule_name}\": \"/path/to/{stem}.sif\"' or "
            f"'\"{stem}\": \"/path/to/{stem}.sif\"' or "
            f"'\"env_dir\": \"/path/to/env\"' to the 'env' section of your config."
        )

    # 3. Per-rule container directory (rule-specific SIF wins over module SIF)
    if rule_name:
        rule_sif = os.path.join(env_dir, rule_name, f"{stem}.sif")
        if os.path.isfile(rule_sif):
            return rule_sif

    # 4. Fallback: env_dir / <module_dir> / <stem>.sif
    module_dir = os.path.basename(workflow.basedir)
    return os.path.join(env_dir, module_dir, f"{stem}.sif")
