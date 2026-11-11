#!/usr/bin/env python3
"""Shuffle fiber intervals to build a null distribution for FDR peak calling.

Reads fiber alignment intervals (BED3+, e.g. from `bedtools bamtobed` or an
awk conversion of `samtools view`), moves every interval to a random position
within its chromosome (matching how `ft call-peaks` scores real vs shuffled
coverage), and writes the shuffled BED consumed by `ft call-peaks --shuffled`.
Stdlib-only, deterministic with --seed.

Usage:
    python shuffle_bed.py -i fibers.bed -g genome.fa.fai -o shuffled.bed --seed 42
"""
from __future__ import annotations

import argparse
import logging
import random
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence, TextIO

# Copied verbatim from src/common/util/LogUtil.py (module bin scripts are
# self-contained — see modules.md 脚本自包含规范 #0).


def setup_logger(
    logger_name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    stream: Optional[TextIO] = None,
) -> logging.Logger:
    """Create and configure a logger with stream/file handlers.

    Parameters
    ----------
    logger_name : str
        Logger name.
    level : int, default=logging.INFO
        Logging level.
    log_file : Optional[str], default=None
        Optional log file path. If provided, file logging is enabled.
    stream : Optional[TextIO], default=None
        Console output stream for the stream handler. None keeps the
        logging default (sys.stderr); pass sys.stdout to log to stdout.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(logger_name)
    logger.propagate = True
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(level)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


logger = setup_logger("ShuffleBed", stream=sys.stdout)


def read_chrom_sizes(path: str) -> Dict[str, int]:
    """Parse a .fai/.chrom.sizes file into {chrom: length}.

    Skips malformed lines (fewer than 2 tab fields or non-integer length).
    """
    sizes: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                try:
                    sizes[f[0]] = int(f[1])
                except ValueError:
                    continue
    return sizes


def run(args: argparse.Namespace) -> int:
    """Shuffle every input interval and write the null BED.

    Intervals are streamed line by line (fiber BEDs run to millions of
    rows; no whole-file buffering). Intervals on chromosomes absent from
    the sizes file, or at least as long as their chromosome, are skipped
    (they cannot be placed) and counted in the summary log.
    """
    sizes = read_chrom_sizes(str(args.genome))
    chrom_names = list(sizes)  # only needed for genome-wide mode
    rng = random.Random(args.seed)
    n_in, n_out, n_skip = 0, 0, 0

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with open(args.input) as fh, open(tmp, "w") as out:
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
            if args.keep_chrom:
                # uniform within the original chromosome (ft call-peaks
                # compares real vs shuffled coverage per chromosome)
                if length >= gsize:
                    n_skip += 1
                    continue
                new_start = rng.randint(0, gsize - length)
            else:
                # genome-wide: also redraw the chromosome
                chrom = rng.choice(chrom_names)
                gsize = sizes[chrom]
                if length >= gsize:
                    n_skip += 1
                    continue
                new_start = rng.randint(0, gsize - length)
            name = f[3] if len(f) > 3 else f"shuf_{n_out}"
            out.write(f"{chrom}\t{new_start}\t{new_start + length}\t{name}\n")
            n_out += 1
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
    p.add_argument("-k", "--keep-chrom", action="store_true", default=True,
                   help="shuffle within each chromosome (default; matches ft call-peaks)")
    p.add_argument("-K", "--no-keep-chrom", dest="keep_chrom", action="store_false",
                   help="shuffle genome-wide, also redrawing the chromosome")
    p.add_argument("-s", "--seed", type=int, default=42, help="random seed")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.log:
        global logger
        logger = setup_logger("ShuffleBed", log_file=str(args.log), stream=sys.stdout)
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
