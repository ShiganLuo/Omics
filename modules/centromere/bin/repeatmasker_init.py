#!/usr/bin/env python3
"""Initialize the RepeatMasker library cache once per environment.

Links a custom Dfam FamDB h5 into RepeatMasker's famdb library directory
(RepeatMasker cannot take a custom h5 path directly), then warms the library
cache with a tiny dummy annotation so concurrent per-sample runs do not race
on makeblastdb. The marker file is written last, so downstream rules start
only after the cache is ready.
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Optional

logger = logging.getLogger("repeatmasker_init")


def setup_logging() -> None:
    """Configure stderr logging for the engine entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def famdb_lib_dir(repeatmasker: str) -> Optional[str]:
    """Locate RepeatMasker's famdb library directory from its executable.

    Args:
        repeatmasker: RepeatMasker executable (name or path).

    Returns:
        The ``share/RepeatMasker/Libraries/famdb`` directory derived from the
        executable location, or None when the executable is not resolvable.
    """
    exe = shutil.which(repeatmasker)
    if exe is None:
        return None
    return os.path.join(os.path.dirname(exe), "../share/RepeatMasker/Libraries/famdb")


def link_dfam_h5(dfam_h5: str, repeatmasker: str) -> None:
    """Link the custom Dfam h5 into the famdb library directory.

    RepeatMasker picks up FamDB files by name inside its library directory.
    Linking is skipped (with a log message) when the target already exists or
    the library directory is not host-visible (container env: expect the
    FamDB to be provided inside the env or via bind mount).

    Args:
        dfam_h5: custom Dfam FamDB h5 file.
        repeatmasker: RepeatMasker executable (name or path).

    Raises:
        FileNotFoundError: the RepeatMasker executable cannot be resolved.
    """
    if not os.path.exists(dfam_h5):
        logger.warning(
            "No valid Dfam h5 provided for RepeatMasker (%s missing), using "
            "default libraries. Please ensure the default libraries include "
            "relevant satellite DNA for the target species.", dfam_h5)
        return
    logger.warning(
        "RepeatMasker does not support specifying a custom Dfam h5 directly; "
        "linking it into the library dir under a name RepeatMasker "
        "recognizes. Ensure the Dfam h5 matches your RepeatMasker version and "
        "includes relevant satellite DNA for the target species.")
    lib_dir = famdb_lib_dir(repeatmasker)
    if lib_dir is None:
        raise FileNotFoundError(
            f"RepeatMasker executable not found: {repeatmasker}; "
            f"set 'Procedure.RepeatMasker' in the workflow config JSON.")
    target = os.path.join(lib_dir, os.path.basename(dfam_h5))
    try:
        if os.path.exists(target):
            logger.info(
                "Custom Dfam h5 already exists in RepeatMasker library "
                "directory: %s, skipping linking.", target)
        elif not os.path.exists(lib_dir):
            logger.warning(
                "RepeatMasker library directory not visible from host "
                "(container env): %s. Skipping linking; expecting FamDB "
                "inside the env or via bind mount.", lib_dir)
        else:
            os.symlink(dfam_h5, target)
            logger.info(
                "Linked custom Dfam h5 to RepeatMasker library directory: %s",
                target)
    except OSError as exc:
        logger.warning(
            "Could not link custom Dfam h5 into %s: %s. Continuing; expecting "
            "FamDB inside the env or via bind mount.", lib_dir, exc)


def warm_cache(repeatmasker: str, species: str, threads: int) -> None:
    """Build the RepeatMasker library cache with a tiny dummy annotation.

    Args:
        repeatmasker: RepeatMasker executable.
        species: RepeatMasker -species term (e.g. 'Mus musculus').
        threads: parallel workers for RepeatMasker (-pa).

    Raises:
        RuntimeError: RepeatMasker exited non-zero.
    """
    tmpdir = tempfile.mkdtemp()
    try:
        dummy_fa = os.path.join(tmpdir, "dummy.fa")
        with open(dummy_fa, "w") as f:
            f.write(">dummy\nATCGATCGATCG\n")
        cmd = [
            repeatmasker,
            "-species", species,
            "-pa", str(threads),
            "-dir", tmpdir,
            "-gff",
            dummy_fa,
        ]
        logger.info("Executing: %s", " ".join(cmd))
        result = subprocess.run(cmd)
        if result.returncode != 0:
            raise RuntimeError(
                f"RepeatMasker cache warmup failed (exit {result.returncode}) "
                f"for species '{species}'; check the rule log for tool output.")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main() -> None:
    """Parse arguments, link the Dfam h5, warm the cache, write the marker."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeatmasker", default="RepeatMasker",
                        help="RepeatMasker executable (default: RepeatMasker)")
    parser.add_argument("--species", default="",
                        help="RepeatMasker -species term (e.g. 'Mus musculus'); "
                             "set via Params.RepeatMasker.species")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--dfam-h5", dest="dfam_h5", default="",
                        help="Custom Dfam FamDB h5 to link into the library dir")
    parser.add_argument("--marker", required=True,
                        help="Marker file to touch when the cache is ready")
    args = parser.parse_args()
    setup_logging()

    if not args.species:
        raise RuntimeError(
            "Missing RepeatMasker species: set 'Params.RepeatMasker.species' "
            "in the workflow config JSON (never guess a species).")

    link_dfam_h5(args.dfam_h5, args.repeatmasker)
    warm_cache(args.repeatmasker, args.species, args.threads)
    os.makedirs(os.path.dirname(args.marker) or ".", exist_ok=True)
    with open(args.marker, "w") as f:
        f.write("done\n")
    logger.info("RepeatMasker library cache initialized.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("%s", exc)
        sys.exit(1)
