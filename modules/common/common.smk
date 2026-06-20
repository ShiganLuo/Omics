"""Common utilities for all modules.

This file is included by modules to provide shared utilities.
Usage:
    include: "../common/common.smk"

Then you can use:
    - setup_logger: Create and configure a logger
    - time, shutil, os, sys: Standard library modules
    - ROOT_DIR: Project root directory from config
"""

import sys
import os
import time
import shutil

# Get ROOT_DIR from config (set by run.py)
ROOT_DIR = config.get("ROOT_DIR", ".")

# Ensure src directory is in sys.path for importing common modules
_src_dir = os.path.join(ROOT_DIR, "src")
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

# Import common utilities
try:
    from common.LogUtil import setup_logger
except ImportError as e:
    raise ImportError(
        f"Failed to import common.LogUtil.setup_logger. "
        f"Please ensure:\n"
        f"1. ROOT_DIR is set correctly in config (current: {ROOT_DIR})\n"
        f"2. Directory exists: {_src_dir}\n"
        f"3. File exists: {os.path.join(_src_dir, 'common', 'LogUtil.py')}\n"
        f"Original error: {e}"
    )
