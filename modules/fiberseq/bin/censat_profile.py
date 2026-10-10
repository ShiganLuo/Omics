#!/usr/bin/env python3
"""Centromere/repeat-region focused Fiber-seq and methylation profiles.

For each region class in a region BED (name column, e.g. CDR_core,
alphaSat_HOR, flank_10kb), computes:
  - chromatin accessibility: fraction of region covered by MSPs >= --msp-min-len
    (the >50 bp MSP proportion metric of Nat Genet 2025 [3]), pooled and
    per-fiber mean,
  - mCpG fraction from a modkit pileup bedGraph (recommended, ONT path [3])
    or methylated-CpG mark density from ft extract cpg.bed.gz,
  - rolling binned profiles across region +/- flank for accessibility and
    methylation (e.g. CDR core vs flanking boundaries).

Usage:
    python censat_profile.py --regions censat.bed --msp sample.msp.bed.gz \
        --nuc sample.nuc.bed.gz --pileup sample.modkit.bedGraph \
        -o results/censat --sample sample
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

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
    RegionIndex,
    block_signal,
    read_bed,
    safe_div,
    write_tsv,
)


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("CenSatProfile", log_file=str(log_file))
    else:
        lg = setup_logger("CenSatProfile")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def read_bedgraph(path: str) -> Dict[str, List[Tuple[int, int, float]]]:
    """Read a bedGraph into chrom -> [(start, end, value)]."""
    out: Dict[str, List[Tuple[int, int, float]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out[f[0]].append((int(f[1]), int(f[2]), float(f[3])))
    for chrom in out:
        out[chrom].sort()
    return out


def bedgraph_mean(track, chrom: str, start: int, end: int) -> float:
    vals = [v for s, e, v in track.get(chrom, []) if s < end and e > start]
    return float(np.mean(vals)) if vals else float("nan")


def msp_coverage_fraction(msp, index: RegionIndex, chrom, start, end,
                          msp_min_len: int) -> Tuple[float, float]:
    """(pooled coverage fraction by long MSPs, per-fiber mean coverage fraction)."""
    length = max(1, end - start)
    cover = np.zeros(length, dtype=np.int64)
    per_fiber = defaultdict(lambda: np.zeros(length, dtype=np.int64))
    for i in index.query(chrom, start, end):
        r = msp[i]
        for bs, be in r.blocks:
            if be - bs < msp_min_len:
                continue
            s = max(bs, start) - start
            e = min(be, end) - start
            if s < e:
                cover[s:e] = 1
                per_fiber[r.name][s:e] = 1
    pooled = float(np.mean(cover))
    fiber_means = [float(np.mean(v)) for v in per_fiber.values()] or [0.0]
    return pooled, float(np.mean(fiber_means))


def run(args: argparse.Namespace) -> int:
    regions = read_bed(str(args.regions))
    msp = read_bed(str(args.msp))
    nuc = read_bed(str(args.nuc))
    msp_i, nuc_i = RegionIndex(msp), RegionIndex(nuc)
    pileup = read_bedgraph(str(args.pileup)) if args.pileup else None
    cpg = read_bed(str(args.cpg)) if args.cpg else None
    cpg_i = RegionIndex(cpg) if cpg else None
    logger.info(f"regions={len(regions)} msp={len(msp)} nuc={len(nuc)} "
                f"pileup={'yes' if pileup else 'no'} cpg={len(cpg) if cpg else 0}")

    class_rows = []
    prof_rows = []
    class_acc: Dict[str, List[float]] = defaultdict(list)
    class_mcpg: Dict[str, List[float]] = defaultdict(list)

    for r in regions:
        cls = r.name
        pooled, fiber_mean = msp_coverage_fraction(
            msp, msp_i, r.chrom, r.start, r.end, args.msp_min_len)
        nuc_cov = len({nuc[i].name for i in nuc_i.query(r.chrom, r.start, r.end)})
        msp_cov = len({msp[i].name for i in msp_i.query(r.chrom, r.start, r.end)})
        if pileup is not None:
            mcpg = bedgraph_mean(pileup, r.chrom, r.start, r.end)
        elif cpg_i is not None:
            marks = sum(
                1 for i in cpg_i.query(r.chrom, r.start, r.end)
                for bs, be in cpg[i].blocks
                if bs < r.end and be > r.start
            )
            mcpg = safe_div(marks, max(1, r.end - r.start)) * 1000.0
        else:
            mcpg = float("nan")
        class_rows.append([
            r.chrom, r.start, r.end, cls,
            round(pooled, 6), round(fiber_mean, 6), msp_cov, nuc_cov,
            round(mcpg, 6) if mcpg == mcpg else "NA",
        ])
        class_acc[cls].append(pooled)
        if mcpg == mcpg:
            class_mcpg[cls].append(mcpg)

        # rolling profile across region +/- flank
        ws, we = r.start - args.flank, r.end + args.flank
        acc_sig = block_signal(msp, r.chrom, ws, we, uniq_fibers=True)
        nbins = max(1, (we - ws) // args.bin)
        for b in range(nbins):
            s0, e0 = b * args.bin, min((b + 1) * args.bin, we - ws)
            acc_bin = float(np.mean(acc_sig[s0:e0])) if e0 > s0 else float("nan")
            pos = ws + (s0 + e0) // 2 - (r.start + r.end) // 2
            meth_bin = (bedgraph_mean(pileup, r.chrom, ws + s0, ws + e0)
                        if pileup is not None else float("nan"))
            prof_rows.append([cls, r.start, r.end, pos, round(acc_bin, 6),
                              round(meth_bin, 6) if meth_bin == meth_bin else "NA"])

    if not class_rows:
        logger.error("no regions analyzed")
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(str(args.outdir / f"{args.sample}.censat_regions.tsv"),
              ["chrom", "start", "end", "region_class", "msp_pooled_frac",
               "msp_fiber_mean_frac", "n_msp_fibers", "n_nuc_fibers", "mcpg"],
              class_rows)

    summary_rows = []
    for cls in sorted(class_acc):
        acc = class_acc[cls]
        mc = class_mcpg.get(cls, [])
        summary_rows.append([
            cls, len(acc), round(float(np.mean(acc)), 6),
            round(float(np.mean(mc)), 6) if mc else "NA",
        ])
    write_tsv(str(args.outdir / f"{args.sample}.censat_class_summary.tsv"),
              ["region_class", "n_regions", "mean_msp_frac", "mean_mcpg"],
              summary_rows)

    write_tsv(str(args.outdir / f"{args.sample}.censat_profile.tsv"),
              ["region_class", "region_start", "region_end", "position",
               "accessibility", "mcpg"], prof_rows)
    _plot(args, prof_rows)
    logger.info(f"Wrote censat tables and plot to {args.outdir}")
    return 0


def _plot(args, prof_rows) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_class: Dict[str, Dict[int, List[float]]] = defaultdict(lambda: defaultdict(list))
    for cls, _, _, pos, acc, _ in prof_rows:
        by_class[cls][pos].append(acc)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for cls, posmap in by_class.items():
        xs = sorted(posmap)
        ys = [float(np.nanmean(posmap[x])) for x in xs]
        ax.plot(xs, ys, label=cls)
    ax.axvline(0, ls="--", lw=0.8, color="grey")
    ax.set_xlabel("position relative to region center (bp)")
    ax.set_ylabel("accessibility (MSP fiber fraction)")
    ax.set_title(f"{args.sample} region-class accessibility profile")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(args.outdir / f"{args.sample}.censat_profile.pdf")
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="censat_profile.py",
        description="Centromere/repeat-region accessibility and methylation profiles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-r", "--regions", type=Path, required=True,
                   help="region BED (name column = region class)")
    p.add_argument("-s", "--msp", type=Path, required=True, help="MSP BED12 (.bed.gz)")
    p.add_argument("-u", "--nuc", type=Path, required=True, help="nucleosome BED12 (.bed.gz)")
    p.add_argument("-p", "--pileup", type=Path, default=None,
                   help="modkit pileup bedGraph (mCpG/CpG per position)")
    p.add_argument("-c", "--cpg", type=Path, default=None,
                   help="methylated CpG marks BED12 (.bed.gz) if no pileup")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-m", "--msp-min-len", type=int, default=50,
                   help="min MSP length to count as accessible [3]")
    p.add_argument("-f", "--flank", type=int, default=10000,
                   help="profile flank around region (bp)")
    p.add_argument("-b", "--bin", type=int, default=500, help="profile bin size (bp)")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    for f in filter(None, [args.regions, args.msp, args.nuc, args.pileup, args.cpg]):
        if not Path(f).exists():
            logger.error(f"input not found: {f}")
            return 1
    if args.pileup is None and args.cpg is None:
        logger.warning("no --pileup/--cpg given; methylation columns will be NA")
    if args.dry_run:
        logger.info(f"[DRY-RUN] would profile {args.regions} -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
