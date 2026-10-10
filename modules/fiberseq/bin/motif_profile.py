#!/usr/bin/env python3
"""Aggregate m6A / accessibility / nucleosome profiles around motif sites.

For every motif occurrence (BED), computes per-position signals across the
window (default +/-250 bp): fiber coverage, m6A fraction, MSP accessibility
fraction, nucleosome occupancy, and a per-site footprint score
(flank maxima minus motif-center signal, GC-robust by flank normalization).

Usage:
    python motif_profile.py --motif-sites ctcf_sites.bed \
        --m6a sample.m6a.bed.gz --msp sample.msp.bed.gz --nuc sample.nuc.bed.gz \
        -o results/sample/motif --sample sample -f 250
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import (  # noqa: E402
    block_signal,
    coverage_signal,
    read_bed,
    safe_div,
    write_tsv,
)


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("MotifProfile", log_file=str(log_file))
    else:
        lg = setup_logger("MotifProfile")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def site_signals(site, m6a, msp, nuc, flank):
    """Windowed per-base signals centered on the motif site midpoint."""
    center = (site.start + site.end) // 2
    ws, we = center - flank, center + flank
    cov = coverage_signal(nuc, msp, site.chrom, ws, we)
    m6a_sig = block_signal(m6a, site.chrom, ws, we, uniq_fibers=True)
    acc = block_signal(msp, site.chrom, ws, we, uniq_fibers=True)
    occ = block_signal(nuc, site.chrom, ws, we, uniq_fibers=True)
    return cov, m6a_sig, acc, occ


def footprint_score(m6a_frac: np.ndarray, flank: int, motif_half: int) -> float:
    """Flank-max minus center-sum score over the motif window (paper [1] style).

    score = max(left flank) + max(right flank) - sum(center motif window)
    using the per-position m6A fraction signal.
    """
    n = len(m6a_frac)
    center = n // 2
    left = m6a_frac[: center - motif_half]
    right = m6a_frac[center + motif_half:]
    center_sig = m6a_frac[center - motif_half: center + motif_half]
    if len(left) == 0 or len(right) == 0:
        return 0.0
    return float(np.max(left) + np.max(right) - np.sum(center_sig))


def run(args: argparse.Namespace) -> int:
    m6a = read_bed(str(args.m6a))
    msp = read_bed(str(args.msp))
    nuc = read_bed(str(args.nuc))
    sites = read_bed(str(args.motif_sites))
    logger.info(f"Loaded {len(sites)} motif sites, m6a={len(m6a)}, "
                f"msp={len(msp)}, nuc={len(nuc)}")

    flank = args.flank
    width = 2 * flank
    sum_cov = np.zeros(width, dtype=np.float64)
    sum_m6a = np.zeros(width, dtype=np.float64)
    sum_acc = np.zeros(width, dtype=np.float64)
    sum_occ = np.zeros(width, dtype=np.float64)
    used = 0
    score_rows = []

    for site in sites:
        cov, m6a_sig, acc, occ = site_signals(site, m6a, msp, nuc, flank)
        if np.max(cov) < args.min_fibers:
            continue
        used += 1
        sum_cov += cov
        sum_m6a += m6a_sig
        sum_acc += acc
        sum_occ += occ
        with np.errstate(divide="ignore", invalid="ignore"):
            m6a_frac = np.where(cov > 0, m6a_sig / np.maximum(cov, 1), 0.0)
        half = max(1, (site.end - site.start) // 2)
        score_rows.append([
            site.chrom, site.start, site.end, site.name,
            round(footprint_score(m6a_frac, flank, half), 6),
            round(float(np.mean(m6a_frac)), 6),
        ])

    logger.info(f"Used {used}/{len(sites)} sites (min_fibers={args.min_fibers})")
    if used == 0:
        logger.error("No motif site passed min_fibers filter; no outputs written")
        return 1

    with np.errstate(divide="ignore", invalid="ignore"):
        m6a_frac = np.where(sum_cov > 0, sum_m6a / np.maximum(sum_cov, 1), np.nan)
        acc_frac = np.where(sum_cov > 0, sum_acc / np.maximum(sum_cov, 1), np.nan)
        occ_frac = np.where(sum_cov > 0, sum_occ / np.maximum(sum_cov, 1), np.nan)

    args.outdir.mkdir(parents=True, exist_ok=True)
    rel = np.arange(-flank, flank)
    rows = [
        [int(rel[i]), int(sum_cov[i]), round(float(m6a_frac[i]), 6),
         round(float(acc_frac[i]), 6), round(float(occ_frac[i]), 6)]
        for i in range(width)
    ]
    profile_tsv = args.outdir / f"{args.sample}.motif_profile.tsv"
    write_tsv(str(profile_tsv),
              ["position", "fiber_coverage", "m6a_frac", "accessibility_frac",
               "nuc_occupancy_frac"], rows)

    scores_tsv = args.outdir / f"{args.sample}.motif_footprint_scores.tsv"
    write_tsv(str(scores_tsv),
              ["chrom", "start", "end", "motif", "footprint_score", "mean_m6a_frac"],
              score_rows)

    _plot(args, rel, m6a_frac, acc_frac, occ_frac)
    logger.info(f"Wrote {profile_tsv}, {scores_tsv} and profile PDF")
    return 0


def _plot(args, rel, m6a_frac, acc_frac, occ_frac) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(rel, m6a_frac, label="m6A fraction", color="#c0392b")
    ax.plot(rel, acc_frac, label="accessibility (MSP)", color="#2980b9")
    ax.plot(rel, occ_frac, label="nucleosome occupancy", color="#27ae60")
    ax.axvline(0, ls="--", lw=0.8, color="grey")
    ax.set_xlabel(f"position relative to motif center (bp)")
    ax.set_ylabel("fraction of fibers")
    ax.set_title(f"{args.sample} motif-centered profile")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out = args.outdir / f"{args.sample}.motif_profile.pdf"
    fig.savefig(out)
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="motif_profile.py",
        description="Aggregate m6A/accessibility/nucleosome profiles around motif sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-m", "--motif-sites", type=Path, required=True,
                   help="motif occurrence BED (e.g. from HOMER findMotifs -find)")
    p.add_argument("-a", "--m6a", type=Path, required=True, help="m6A BED12 (.bed.gz)")
    p.add_argument("-s", "--msp", type=Path, required=True, help="MSP BED12 (.bed.gz)")
    p.add_argument("-u", "--nuc", type=Path, required=True, help="nucleosome BED12 (.bed.gz)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-f", "--flank", type=int, default=250,
                   help="window half-size around motif center (bp)")
    p.add_argument("-i", "--min-fibers", type=int, default=10,
                   help="minimum fiber coverage per site to include")
    p.add_argument("-l", "--log", type=Path, default=None,
                   help="log file path (default: log to stdout)")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in (args.motif_sites, args.m6a, args.msp, args.nuc):
        if not f.exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would profile {args.motif_sites} -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
