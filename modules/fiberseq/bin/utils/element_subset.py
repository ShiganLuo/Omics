#!/usr/bin/env python3
"""Extract region-overlapping subsets of big Fiber-seq element BEDs.

The per-fiber element BEDs (`ft extract` BED12, m6a/msp/nuc/cpg/fire) are
multi-GB bgzip files; region-focused runs (quant_matrix.py, censat_profile.py)
would stream each file whole. This script takes one compressed element BED and
a small regions BED, and writes the subset of records whose span overlaps a
region (expanded by --margin), preserving input order and de-duplicating
records that touch several query intervals. Feeding the subset to those
engines is exactly equivalent to whole-file streaming (their per-anchor
overlap filter does the rest).

Compression is delegated to `zcat` (C-speed); the file itself is only read.

Usage:
    python element_subset.py --elements sample.m6a.bed.gz --regions telocen.bed \
        --margin 15000 --outdir subset/ --sample S1
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_BIN_DIR = Path(__file__).resolve().parent
_ROOT_DIR = _BIN_DIR.parent.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("ElementSubset", log_file=str(log_file))
    else:
        lg = setup_logger("ElementSubset")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def merged_query_intervals(regions_path: Path, margin: int
                           ) -> Dict[str, List[Tuple[int, int]]]:
    """Per-chromosome sorted disjoint [start, end) intervals, margin-applied."""
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open(regions_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            start = max(0, int(f[1]) - margin)
            end = int(f[2]) + margin
            by_chrom[f[0]].append((start, end))
    out: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, spans in by_chrom.items():
        spans.sort()
        merged = [list(spans[0])]
        for s, e in spans[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[chrom] = [(s, e) for s, e in merged]
    return out


def overlaps(spans: Sequence[Tuple[int, int]], start: int, end: int) -> bool:
    """Any interval in the sorted disjoint list overlapping [start, end)."""
    lo, hi = 0, len(spans)
    while lo < hi:  # first span with span_end > start
        mid = (lo + hi) // 2
        if spans[mid][1] <= start:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(spans) and spans[lo][0] < end


def subset_one(element: Path, regions_path: Path, margin: int,
               outdir: Path) -> Path:
    queries = merged_query_intervals(regions_path, margin)
    n_query = sum(len(v) for v in queries.values())
    out = outdir / f"{element.name.replace('.bed.gz', '')}.subset.bed"
    outdir.mkdir(parents=True, exist_ok=True)

    opener = (lambda: subprocess.Popen(
        ["zcat", str(element)], stdout=subprocess.PIPE, text=True))
    seen = set()
    n_in = n_out = 0
    with opener() as proc, open(out, "w") as dst:
        assert proc.stdout is not None
        for line in proc.stdout:
            n_in += 1
            if line.startswith(("#", "track", "browser")):
                continue
            c1 = line.find("\t")
            c2 = line.find("\t", c1 + 1)
            c3 = line.find("\t", c2 + 1)
            chrom = line[:c1]
            start = int(line[c1 + 1:c2])
            end = int(line[c2 + 1:c3])
            spans = queries.get(chrom)
            if not spans or not overlaps(spans, start, end):
                continue
            c4 = line.find("\t", c3 + 1)
            name = line[c3 + 1:c4] if c4 > 0 else line[c3 + 1:].rstrip("\n")
            key = (chrom, start, end, name)
            if key in seen:
                continue
            seen.add(key)
            dst.write(line if line.endswith("\n") else line + "\n")
            n_out += 1
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"zcat failed on {element} (rc={rc})")
    logger.info(f"{element.name}: {n_out}/{n_in} records kept "
                f"({n_query} query intervals, margin={margin}) -> {out}")
    return out


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="element_subset.py",
        description="Region-overlap subset of a compressed element BED.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--elements", type=Path, nargs="+", required=True,
                   help="element BED(.gz) files")
    p.add_argument("--regions", type=Path, required=True,
                   help="regions BED (query)")
    p.add_argument("--margin", type=int, default=0,
                   help="expand each query region by this many bp")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sample", type=str, default="sample",
                   help="sample name (for the log only)")
    p.add_argument("-l", "--log", type=Path, default=None)
    p.add_argument("-n", "--dry-run", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in list(args.elements) + [args.regions]:
        if not Path(f).exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would subset {len(args.elements)} element files "
                    f"against {args.regions} -> {args.outdir}")
        return 0
    for element in args.elements:
        subset_one(Path(element), args.regions, args.margin, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
