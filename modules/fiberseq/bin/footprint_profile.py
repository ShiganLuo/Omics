#!/usr/bin/env python3
"""TF footprinting on single fibers around motif sites.

Classifies every fiber over a motif site into three categories [1]:
  (i)   MSP covering the motif with an m6A-free footprint in the motif core,
  (ii)  MSP covering the motif without a footprint,
  (iii) nucleosome occupying the motif.
Then aggregates m6A fraction profiles per category and compares per-site
footprint scores between categories with a Kolmogorov-Smirnov test.

Usage:
    python footprint_profile.py -m ctcf_sites.bed -a sample.m6a.bed.gz \
        -s sample.msp.bed.gz -u sample.nuc.bed.gz -o out/ --sample s1
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

from src.common.util.FiberlibUtil import (  # noqa: E402
    BedRecord,
    block_signal,
    coverage_signal,
    read_bed,
    write_tsv,
)

CATEGORIES = ("msp_footprint", "msp_no_footprint", "nucleosome")


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("FootprintProfile", log_file=str(log_file))
    else:
        lg = setup_logger("FootprintProfile")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def classify_fibers(site, m6a, msp, nuc) -> Dict[str, set]:
    """Classify fibers covering the motif core into the three categories."""
    core_s, core_e = site.start, site.end
    m6a_by_fiber: Dict[str, List[tuple]] = {}
    for r in m6a:
        if r.chrom != site.chrom:
            continue
        for bs, be in r.blocks:
            if bs < core_e and be > core_s:
                m6a_by_fiber.setdefault(r.name, []).append((bs, be))
    msp_fibers, nuc_fibers = set(), set()
    for r in msp:
        if r.chrom != site.chrom:
            continue
        for bs, be in r.blocks:
            if bs <= core_s and be >= core_e:
                msp_fibers.add(r.name)
            elif bs < core_e and be > core_s:
                msp_fibers.add(r.name)
    for r in nuc:
        if r.chrom != site.chrom:
            continue
        for bs, be in r.blocks:
            if bs < core_e and be > core_s:
                nuc_fibers.add(r.name)

    out = {c: set() for c in CATEGORIES}
    for fiber in msp_fibers:
        has_m6a = any(be > core_s and bs < core_e for bs, be in m6a_by_fiber.get(fiber, []))
        out["msp_no_footprint" if has_m6a else "msp_footprint"].add(fiber)
    for fiber in nuc_fibers:
        if fiber not in out["msp_footprint"] and fiber not in out["msp_no_footprint"]:
            out["nucleosome"].add(fiber)
    return out


def footprint_score(m6a_frac: np.ndarray, motif_half: int) -> float:
    n = len(m6a_frac)
    center = n // 2
    left = m6a_frac[: center - motif_half]
    right = m6a_frac[center + motif_half:]
    center_sig = m6a_frac[center - motif_half: center + motif_half]
    if len(left) == 0 or len(right) == 0:
        return 0.0
    return float(np.max(left) + np.max(right) - np.sum(center_sig))


def run(args: argparse.Namespace) -> int:
    from scipy.stats import ks_2samp

    sites = read_bed(str(args.motif_sites))
    m6a = read_bed(str(args.m6a))
    msp = read_bed(str(args.msp))
    nuc = read_bed(str(args.nuc))
    logger.info(f"sites={len(sites)} m6a={len(m6a)} msp={len(msp)} nuc={len(nuc)}")

    flank = args.flank
    width = 2 * flank
    agg = {c: {"cov": np.zeros(width), "m6a": np.zeros(width)} for c in CATEGORIES}
    score_sets: Dict[str, List[float]] = {c: [] for c in CATEGORIES}
    cat_rows = []

    for site in sites:
        cats = classify_fibers(site, m6a, msp, nuc)
        counts = {c: len(cats[c]) for c in CATEGORIES}
        cat_rows.append([site.chrom, site.start, site.end, site.name,
                         counts["msp_footprint"], counts["msp_no_footprint"],
                         counts["nucleosome"]])
        center = (site.start + site.end) // 2
        ws, we = center - flank, center + flank
        for cat, fibers in cats.items():
            if not fibers:
                continue
            m6a_sub = [r for r in m6a if r.name in fibers]
            msp_sub = [r for r in msp if r.name in fibers]
            nuc_sub = [r for r in nuc if r.name in fibers]
            cov = coverage_signal(nuc_sub, msp_sub, site.chrom, ws, we)
            m6a_sig = block_signal(m6a_sub, site.chrom, ws, we, uniq_fibers=True)
            agg[cat]["cov"] += cov
            agg[cat]["m6a"] += m6a_sig
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = np.where(cov > 0, m6a_sig / np.maximum(cov, 1), 0.0)
            score_sets[cat].append(
                footprint_score(frac, max(1, (site.end - site.start) // 2))
            )

    if sum(len(v) for v in score_sets.values()) == 0:
        logger.error("no fiber classified at any site; no outputs")
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(str(args.outdir / f"{args.sample}.footprint_categories.tsv"),
              ["chrom", "start", "end", "motif", *CATEGORIES], cat_rows)

    rel = np.arange(-flank, flank)
    rows = []
    for i in range(width):
        row = [int(rel[i])]
        for c in CATEGORIES:
            cov = agg[c]["cov"][i]
            frac = agg[c]["m6a"][i] / cov if cov > 0 else "NA"
            row.append(round(float(frac), 6) if frac != "NA" else "NA")
        rows.append(row)
    write_tsv(str(args.outdir / f"{args.sample}.footprint_profiles.tsv"),
              ["position"] + [f"m6a_frac_{c}" for c in CATEGORIES], rows)

    score_rows = []
    pairs = [("msp_footprint", "msp_no_footprint"),
             ("msp_footprint", "nucleosome")]
    for c1, c2 in pairs:
        if score_sets[c1] and score_sets[c2]:
            stat, p = ks_2samp(score_sets[c1], score_sets[c2])
            score_rows.append([c1, c2, len(score_sets[c1]), len(score_sets[c2]),
                               round(float(stat), 6), float(p)])
        else:
            score_rows.append([c1, c2, len(score_sets[c1]), len(score_sets[c2]), "NA", "NA"])
    write_tsv(str(args.outdir / f"{args.sample}.footprint_scores.tsv"),
              ["group1", "group2", "n1", "n2", "ks_stat", "ks_pvalue"], score_rows)

    _plot(args, rel, agg)
    logger.info(f"Wrote footprint tables and plot to {args.outdir}")
    return 0


def _plot(args, rel, agg) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"msp_footprint": "#c0392b", "msp_no_footprint": "#2980b9",
              "nucleosome": "#27ae60"}
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for c in CATEGORIES:
        cov = agg[c]["cov"]
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(cov > 0, agg[c]["m6a"] / np.maximum(cov, 1), np.nan)
        ax.plot(rel, frac, label=c, color=colors[c])
    ax.axvline(0, ls="--", lw=0.8, color="grey")
    ax.set_xlabel("position relative to motif center (bp)")
    ax.set_ylabel("m6A fraction of fibers")
    ax.set_title(f"{args.sample} TF footprint")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.outdir / f"{args.sample}.footprint_profile.pdf")
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="footprint_profile.py",
        description="Single-fiber TF footprinting around motif sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-m", "--motif-sites", type=Path, required=True, help="motif site BED")
    p.add_argument("-a", "--m6a", type=Path, required=True, help="m6A BED12 (.bed.gz)")
    p.add_argument("-s", "--msp", type=Path, required=True, help="MSP BED12 (.bed.gz)")
    p.add_argument("-u", "--nuc", type=Path, required=True, help="nucleosome BED12 (.bed.gz)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-f", "--flank", type=int, default=1000,
                   help="profile half-window around motif center (bp)")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
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
        logger.info(f"[DRY-RUN] would profile footprints -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
