#!/usr/bin/env python3
"""Shuffle fiber intervals to build a null distribution for FDR peak calling.

Reads fiber alignment intervals (BED3+, e.g. from `bedtools bamtobed` or an
awk conversion of `samtools view`), moves every interval to a random position
(uniform within its chromosome by default), and writes the shuffled BED for
`ft call-peaks --shuffled` [5]. Stdlib-only, deterministic with --seed.

Usage:
    python shuffle_bed.py -i fibers.bed -g genome.chrom.sizes -o shuffled.bed --seed 1
"""
from __future__ import annotations

import argparse
import logging
import random
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence

_SRC_DIR = Path(__file__).resolve().parent.parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("ShuffleBed", log_file=str(log_file))
    else:
        lg = setup_logger("ShuffleBed")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def read_chrom_sizes(path: str) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    return sizes


def run(args: argparse.Namespace) -> int:
    sizes = read_chrom_sizes(str(args.genome))
    rng = random.Random(args.seed)
    n_in, n_out, n_skip = 0, 0, 0
    out_lines = []
    with open(args.input) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            n_in += 1
            chrom = f[0]
            length = int(f[2]) - int(f[1])
            gsize = sizes.get(chrom)
            if gsize is None:
                n_skip += 1
                continue
            if length >= gsize:
                n_skip += 1
                continue
            if args.keep_chrom:
                new_start = rng.randint(0, gsize - length)
            else:
                chrom = rng.choice(list(sizes))
                gsize = sizes[chrom]
                if length >= gsize:
                    n_skip += 1
                    continue
                new_start = rng.randint(0, gsize - length)
            name = f[3] if len(f) > 3 else f"shuf_{n_out}"
            out_lines.append(f"{chrom}\t{new_start}\t{new_start + length}\t{name}\n")
            n_out += 1

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        fh.writelines(out_lines)
    tmp.replace(out_path)
    logger.info(f"Shuffled {n_out}/{n_in} intervals ({n_skip} skipped) -> {out_path}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="shuffle_bed.py",
        description="Shuffle fiber intervals to build a null for FDR peak calling.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", type=Path, required=True, help="fiber intervals BED")
    p.add_argument("-g", "--genome", type=Path, required=True,
                   help="chrom sizes (.chrom.sizes or .fai)")
    p.add_argument("-o", "--output", type=Path, required=True, help="shuffled BED output")
    p.add_argument("-k", "--keep-chrom", action="store_true",
                   help="shuffle within each chromosome (default: genome-wide)")
    p.add_argument("-s", "--seed", type=int, default=42, help="random seed")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in (args.input, args.genome):
        if not f.exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would shuffle {args.input} -> {args.output}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
