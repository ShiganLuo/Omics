#!/usr/bin/env python3
"""Co-actuation analysis: pairs of regulatory elements active on the same fiber.

For every pair of regions within --max-distance, counts fibers covering both
regions and how many carry a FIRE element at each (or both). Tests enrichment
of co-actuation with Fisher's exact test and compares the observed co-actuated
fiber proportion against the independence expectation (paper [1]).

Usage:
    python coactuation.py -p sample.fire_peaks.bed --fire sample.fire.bed.gz \
        --msp sample.msp.bed.gz --nuc sample.nuc.bed.gz -o out/ --sample s1
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import RegionIndex, bh_fdr, read_bed, safe_div, write_tsv  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("Coactuation", log_file=str(log_file))
    else:
        lg = setup_logger("Coactuation")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def fiber_sets(index: RegionIndex, records, chrom: str, start: int, end: int) -> Set[str]:
    return {records[i].name for i in index.query(chrom, start, end)}


def run(args: argparse.Namespace) -> int:
    from scipy.stats import fisher_exact

    peaks = read_bed(str(args.peaks))
    fire = read_bed(str(args.fire))
    msp = read_bed(str(args.msp))
    nuc = read_bed(str(args.nuc))
    fire_i, msp_i, nuc_i = RegionIndex(fire), RegionIndex(msp), RegionIndex(nuc)
    logger.info(f"peaks={len(peaks)} fire={len(fire)} msp={len(msp)} nuc={len(nuc)}")

    pairs = []
    for i, p1 in enumerate(peaks):
        for p2 in peaks[i + 1:]:
            if p2.chrom != p1.chrom:
                continue
            gap = max(p2.start - p1.end, p1.start - p2.end, 0)
            if gap > args.max_distance:
                continue
            pairs.append((p1, p2))
            if len(pairs) >= args.max_pairs:
                break
        if len(pairs) >= args.max_pairs:
            break
    logger.info(f"Testing {len(pairs)} region pairs (max_distance={args.max_distance})")

    rows = []
    pvals = []
    obs_diffs = []
    for p1, p2 in pairs:
        cov1 = fiber_sets(nuc_i, nuc, p1.chrom, p1.start, p1.end) | \
            fiber_sets(msp_i, msp, p1.chrom, p1.start, p1.end)
        cov2 = fiber_sets(nuc_i, nuc, p2.chrom, p2.start, p2.end) | \
            fiber_sets(msp_i, msp, p2.chrom, p2.start, p2.end)
        both_cov = cov1 & cov2
        if len(both_cov) < args.min_fibers:
            continue
        f1 = fiber_sets(fire_i, fire, p1.chrom, p1.start, p1.end) & both_cov
        f2 = fiber_sets(fire_i, fire, p2.chrom, p2.start, p2.end) & both_cov
        b = len(f1 & f2)
        a1_only = len(f1) - b
        a2_only = len(f2) - b
        n = len(both_cov) - b - a1_only - a2_only
        table = [[b, a1_only], [a2_only, n]]
        try:
            _, p = fisher_exact(table, alternative="greater")
        except ValueError:
            p = 1.0
        actual = safe_div(b, len(both_cov))
        expected = safe_div(len(f1), len(both_cov)) * safe_div(len(f2), len(both_cov))
        obs_diffs.append(actual - expected)
        pvals.append(p)
        rows.append([
            p1.chrom, f"{p1.start}-{p1.end}", p1.name,
            f"{p2.start}-{p2.end}", p2.name,
            len(both_cov), len(f1), len(f2), b, a1_only, a2_only, n,
            round(actual, 6), round(expected, 6), round(actual - expected, 6),
            p,
        ])

    if not rows:
        logger.error("no region pair passed min_fibers; no outputs")
        return 1
    qvals = bh_fdr(pvals)
    for row, q in zip(rows, qvals):
        row.append(round(q, 6))
        row.append("yes" if q < args.alpha else "no")

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{args.sample}.coactuation.tsv"
    write_tsv(str(out),
              ["chrom", "site1", "name1", "site2", "name2", "n_covering_fibers",
               "n_fire_site1", "n_fire_site2", "n_both", "n_site1_only",
               "n_site2_only", "n_neither", "actual_frac", "expected_frac",
               "diff", "p_value", "q_value", "significant"], rows)

    n_sig = sum(1 for r in rows if r[17] == "yes")
    mean_diff = sum(obs_diffs) / len(obs_diffs)
    summary = args.outdir / f"{args.sample}.coactuation_summary.tsv"
    write_tsv(str(summary),
              ["metric", "value"],
              [["pairs_tested", len(rows)], ["significant_pairs", n_sig],
               ["mean_actual_minus_expected", round(mean_diff, 6)]])
    logger.info(f"Wrote {out} and {summary} ({n_sig} significant pairs)")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="coactuation.py",
        description="Same-fiber co-actuation analysis of region pairs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-p", "--peaks", type=Path, required=True, help="regions/peaks BED")
    p.add_argument("-f", "--fire", type=Path, required=True,
                   help="per-fiber FIRE BED (from `ft fire --extract`)")
    p.add_argument("-s", "--msp", type=Path, required=True, help="MSP BED12 (.bed.gz)")
    p.add_argument("-u", "--nuc", type=Path, required=True, help="nucleosome BED12 (.bed.gz)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-d", "--max-distance", type=int, default=50000,
                   help="max gap between paired regions (bp)")
    p.add_argument("-i", "--min-fibers", type=int, default=10,
                   help="min fibers covering both regions")
    p.add_argument("-m", "--max-pairs", type=int, default=200000,
                   help="cap on tested pairs")
    p.add_argument("-a", "--alpha", type=float, default=0.05, help="FDR threshold")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in (args.peaks, args.fire, args.msp, args.nuc):
        if not f.exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would test co-actuation pairs -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
