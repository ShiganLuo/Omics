#!/usr/bin/env python3
"""Haplotype-specific accessibility differences across regions.

Compares per-region fiber counts between two haplotype quant tables (from
quant_matrix.py run on haplotype-split BAMs) with Fisher's exact test and
log2FC of percent accessibility, following the haplotype analysis in [1][2].

Usage:
    python haplotype_diff.py -1 hap1.region_quant.tsv -2 hap2.region_quant.tsv \
        -o results/sample/haplotype --sample sample
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import bh_fdr, write_tsv  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("HaplotypeDiff", log_file=str(log_file))
    else:
        lg = setup_logger("HaplotypeDiff")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def read_quant(path: Path) -> Dict[str, Dict[str, int]]:
    """Read region_quant.tsv into region key -> counts dict."""
    out: Dict[str, Dict[str, int]] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d = dict(zip(header, f))
            key = f"{d['chrom']}:{d['start']}-{d['end']}|{d['name']}"
            out[key] = {
                "chrom": d["chrom"], "start": d["start"], "end": d["end"],
                "name": d["name"],
                "n_fibers": int(float(d["n_fibers"])),
                "n_accessible": int(float(d["n_accessible"])),
            }
    return out


def run(args: argparse.Namespace) -> int:
    from scipy.stats import fisher_exact

    h1 = read_quant(args.hap1)
    h2 = read_quant(args.hap2)
    keys = [k for k in h1 if k in h2]
    logger.info(f"hap1 regions={len(h1)}, hap2 regions={len(h2)}, shared={len(keys)}")
    if not keys:
        logger.error("no shared regions between haplotype tables")
        return 1

    rows = []
    pvals = []
    for key in keys:
        a, b = h1[key], h2[key]
        if a["n_fibers"] < args.min_fibers or b["n_fibers"] < args.min_fibers:
            continue
        acc1 = a["n_accessible"]
        acc2 = b["n_accessible"]
        closed1 = a["n_fibers"] - acc1
        closed2 = b["n_fibers"] - acc2
        try:
            _, p = fisher_exact([[acc1, closed1], [acc2, closed2]],
                                alternative="two-sided")
        except ValueError:
            p = 1.0
        frac1 = acc1 / a["n_fibers"]
        frac2 = acc2 / b["n_fibers"]
        l2fc = float(np.log2((frac1 + 1e-9) / (frac2 + 1e-9)))
        pvals.append(p)
        rows.append([
            a["chrom"], a["start"], a["end"], a["name"],
            a["n_fibers"], acc1, round(frac1, 6),
            b["n_fibers"], acc2, round(frac2, 6),
            round(l2fc, 4), p,
        ])

    if not rows:
        logger.error(f"no region passed min_fibers={args.min_fibers}")
        return 1
    qvals = bh_fdr(pvals)
    for row, q in zip(rows, qvals):
        row.append(round(q, 6))
        row.append("yes" if q < args.alpha else "no")

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{args.sample}.haplotype_diff.tsv"
    write_tsv(str(out),
              ["chrom", "start", "end", "name", "hap1_n_fibers", "hap1_n_accessible",
               "hap1_frac", "hap2_n_fibers", "hap2_n_accessible", "hap2_frac",
               "log2FC_hap1_vs_hap2", "p_value", "q_value", "significant"], rows)

    n_sig = sum(1 for r in rows if r[13] == "yes")
    write_tsv(str(args.outdir / f"{args.sample}.haplotype_summary.tsv"),
              ["metric", "value"],
              [["regions_tested", len(rows)], ["significant_regions", n_sig]])
    logger.info(f"Wrote {out} ({n_sig} significant regions)")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="haplotype_diff.py",
        description="Haplotype-specific accessibility comparison per region.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-1", "--hap1", type=Path, required=True,
                   help="haplotype 1 region_quant.tsv")
    p.add_argument("-2", "--hap2", type=Path, required=True,
                   help="haplotype 2 region_quant.tsv")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-i", "--min-fibers", type=int, default=10,
                   help="min fibers per region per haplotype")
    p.add_argument("-a", "--alpha", type=float, default=0.05, help="FDR threshold")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in (args.hap1, args.hap2):
        if not f.exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would compare {args.hap1} vs {args.hap2} -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
