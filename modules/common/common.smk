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

def rule_res(name):
    """Resolve resources for rule *name* from config["resources"].

    Merge order: resources.default < resources[name]. Missing keys fall
    back through default; rules with no entry inherit default entirely.
    Expected shape in the workflow config JSON:

        "resources": {
            "default": {"threads": 2, "mem_mb": 8000},
            "fiberseq_track": {"threads": 4, "mem_mb": 65000}
        }

    Usage in a rule:
        threads: rule_res("fiberseq_track").get("threads", 4)
        resources:
            mem_mb = rule_res("fiberseq_track").get("mem_mb", 65000)
    """
    res_cfg = config.get("resources") or {}
    merged = dict(res_cfg.get("default") or {})
    merged.update(res_cfg.get(name) or {})
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


def cmd_dir(logdir, key="group"):
    """Resolve the command-receipt directory <output_root>/log/cmd/<key>.

    The receipt location is a fixed convention: ALWAYS <output_root>/log/cmd,
    regardless of which logdir variant a module received (subworkflows pass
    <root>/log, <root>/log/sample or <root>/log/group). The log root is
    recovered from the given logdir by locating its trailing "log" segment;
    relative logdirs are anchored at ROOT_DIR first.

    key: sample/tissue/contrast/... wildcard value (str) or "group".
    Returns the receipt dir path (created on disk).
    """
    _ld = str(logdir)
    if not os.path.isabs(_ld):
        _ld = os.path.join(ROOT_DIR, _ld)
    _ld = os.path.normpath(_ld)
    # walk up to the .../log segment so .../log/sample and .../log/group
    # both land on <root>/log
    parts = _ld.split(os.sep)
    if "log" in parts:
        _log_root = os.sep.join(parts[: parts.index("log") + 1])
    else:
        _log_root = _ld  # unusual layout: keep the given dir as-is
    _d = os.path.join(_log_root, "cmd", str(key))
    os.makedirs(_d, exist_ok=True)
    return _d


def sif(yaml_filename: str) -> str:
    """Resolve the SIF container path for a conda env YAML.

    Looks up ``config["env"][<yaml_stem>]`` where ``yaml_stem`` is the
    YAML filename without extension (e.g. ``star.yaml`` -> ``star``).
    If not found, falls back to ``<env_dir>/<module_dir>/<yaml_stem>.sif``
    where ``env_dir`` is taken from ``config["env"]["env_dir"]``.

    Parameters
    ----------
    yaml_filename : str
        The conda YAML filename (same value passed to ``conda:`` in a rule).

    Returns
    -------
    str
        Path to the ``.sif`` file.
    """
    stem = os.path.splitext(os.path.basename(yaml_filename))[0]
    logger.debug(f"stem:{stem}, yaml_filename: {yaml_filename}")
    logger.debug(f"_ENV_MAP: {_ENV_MAP}")
    # 1. Explicit mapping in config["env"]
    if stem in _ENV_MAP:
        return _ENV_MAP[stem]

    # 2. Fallback: env_dir / <module_dir> / <stem>.sif
    env_dir = _ENV_MAP.get("env_dir")
    if not env_dir:
        raise ValueError(
            f"SIF path for '{stem}' not found in config['env'] "
            f"and no 'env_dir' fallback set. Please add "
            f"'\"{stem}\": \"/path/to/{stem}.sif\"' or "
            f"'\"env_dir\": \"/path/to/env\"' to the 'env' section of your config."
        )
    module_dir = os.path.basename(workflow.basedir)
    return os.path.join(env_dir, module_dir, f"{stem}.sif")
