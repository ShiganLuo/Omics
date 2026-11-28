"""Discover pseudobulk DEG enrichment units under the pseudobulk output tree.

Scans ``<root>/<tissue>/<counter>/<mode_dir>/deg_<tissue>/<contrast>/<cell_type>/``
for ``<contrast>_<cell_type>_DEG.csv`` tables and writes one
``counter/mode/tissue/contrast/cell_type`` unit line per table to the output
file (empty file when nothing is found).
"""
import argparse
import glob
import json
import logging
import os
import sys
from typing import Dict, List

logger = logging.getLogger("discover_units")


def parse_args(argv: List[str]) -> argparse.Namespace:
    """Parse CLI arguments; modes/mode-dirs arrive as JSON (config-native form)."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True,
                        help="pseudobulk output tree root")
    parser.add_argument("--modes", required=True,
                        help='JSON list of modes to scan, e.g. ["gene","te"]')
    parser.add_argument("--mode-dirs", required=True,
                        help='JSON dict mode -> DEG dir name, e.g. {"gene":"deg_pseudobulk_gene"}')
    parser.add_argument("--out", required=True,
                        help="units.txt output path (one unit per line)")
    return parser.parse_args(argv)


def scan_units(root: str, modes: List[str], mode_dirs: Dict[str, str]) -> List[str]:
    """List 'counter/mode/tissue/contrast/cell_type' units under root.

    A unit is valid only when the DEG table lives at the exact pseudobulk
    contract path (``deg_<tissue>`` mirrors ``<tissue>``, file name is
    ``<contrast>_<cell_type>_DEG.csv``) and its mode dir maps back to a
    requested mode. Returns sorted unit lines for deterministic output.
    """
    mode_of = {v: k for k, v in mode_dirs.items()}
    unknown = sorted(set(modes) - set(mode_dirs))
    if unknown:
        logger.warning(f"modes without a known mode dir are ignored: {unknown}")
    pattern = os.path.join(root, "*", "*", "deg_pseudobulk*",
                           "deg_*", "*", "*", "*_DEG.csv")
    units = set()
    for path in glob.glob(pattern):
        rel = os.path.relpath(path, root).split(os.sep)
        if len(rel) != 7:
            continue
        tissue, counter, mode_dir, deg_tissue, contrast, cell_type, fname = rel
        if deg_tissue != f"deg_{tissue}":
            continue
        if fname != f"{contrast}_{cell_type}_DEG.csv":
            continue
        mode = mode_of.get(mode_dir)
        if mode is None or mode not in modes:
            continue
        units.add(f"{counter}/{mode}/{tissue}/{contrast}/{cell_type}")
    return sorted(units)


def write_units(units: List[str], out_path: str) -> None:
    """Write sorted unit lines to out_path (empty file when no unit found)."""
    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("\n".join(units) + ("\n" if units else ""))


def setup_logging() -> None:
    """Configure logging to stderr with timestamped level format."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        stream=sys.stderr,
    )


def main() -> None:
    """Parse arguments, scan units, write the units file, log the count."""
    setup_logging()
    args = parse_args(sys.argv[1:])
    modes = json.loads(args.modes)
    mode_dirs = json.loads(args.mode_dirs)
    units = scan_units(args.root, modes, mode_dirs)
    write_units(units, args.out)
    logger.info(f"Discovered {len(units)} unit(s) under {args.root}")


if __name__ == "__main__":
    main()
