"""Common utilities for all modules.

This file is included by modules to provide shared utilities.
Usage:
    include: "../common/common.smk"

Then you can use:
    - setup_logger: Create and configure a logger
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
# Container (SIF) path resolution
# ---------------------------------------------------------------------------
_ENV_MAP = config.get("env", {})


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
