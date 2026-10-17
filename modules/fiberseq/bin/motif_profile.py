#!/usr/bin/env python3
"""Aggregate m6A / accessibility / nucleosome profiles around motif sites.

For every motif occurrence (BED), computes per-position signals across the
window (default +/-250 bp): fiber coverage, m6A fraction, MSP accessibility
fraction, nucleosome occupancy, and a per-site footprint score
(flank maxima minus motif-center signal, GC-robust by flank normalization).

Implementation note (vectorized):
    The legacy implementation re-scanned all fiber records overlapping each
    site window separately for every site (~130 ms/site => weeks on full
    motif catalogs). This version streams each element BED exactly once per
    chromosome into per-base coverage arrays, then each site window is a
    numpy slice:

    * m6a / accessibility / nucleosome arrays: plain +1/-1 difference
      arrays. Exact, because `ft extract` BED12 stores one row per fiber
      whose blocks are disjoint, so per-base distinct-fiber coverage equals
      the summed block coverage (the per-window fiber merge of the legacy
      block_signal(uniq_fibers=True) is a no-op for disjoint blocks).
    * fiber coverage (nuc+msp union): per-fiber interval merge across BOTH
      files per chromosome (a fiber carrying both a nucleosome and an MSP
      over the same base must count once), matching the legacy
      coverage_signal() semantics exactly.

    Site order, min_fibers filtering, arithmetic, and rounding are
    unchanged, so outputs are identical to the legacy implementation for
    the same (chromosome-grouped) inputs.

Usage:
    python motif_profile.py --motif-sites ctcf_sites.bed \
        --m6a sample.m6a.bed.gz --msp sample.msp.bed.gz --nuc sample.nuc.bed.gz \
        -o results/sample/motif --sample sample -f 250 [--motif-filter CTCF@,MYC@]
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from array import array
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
_ROOT_DIR = _BIN_DIR.parent.parent.parent
if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402

from src.common.util.FiberlibUtil import stream_bed_records, write_tsv  # noqa: E402
from _fiberprofile_core import ChromBlockStream, IdMapper  # noqa: E402


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


def _load_sites_grouped(path: str, motif_filter: Optional[set]) -> Dict[str, List[Tuple[int, int, str]]]:
    """Stream the sites BED into {chrom: [(start, end, name), ...]} in file order."""
    by_chrom: Dict[str, List[Tuple[int, int, str]]] = {}
    n_all = 0
    for r in stream_bed_records(path):
        n_all += 1
        if motif_filter and r.name.split("@", 1)[0] not in motif_filter:
            continue
        by_chrom.setdefault(r.chrom, []).append((r.start, r.end, r.name))
    if motif_filter:
        kept = sum(len(v) for v in by_chrom.values())
        logger.info(f"motif-filter kept {kept}/{n_all} sites")
    return by_chrom


def run(args: argparse.Namespace) -> int:
    t0 = time.time()
    motif_filter = None
    if args.motif_filter:
        motif_filter = {k.strip() for k in args.motif_filter.split(",") if k.strip()}
    by_chrom = _load_sites_grouped(str(args.motif_sites), motif_filter)
    n_sites = sum(len(v) for v in by_chrom.values())
    logger.info(f"Loaded {n_sites} motif sites")

    flank = args.flank
    width = 2 * flank
    if n_sites == 0:
        logger.error("no motif sites to profile; no outputs written")
        return 1
    if n_sites > 1_000_000 and not motif_filter:
        logger.warning(
            f"{n_sites} sites and no --motif-filter: the per-site scores table "
            f"will be huge; consider --motif-filter"
        )

    paths = {"m6a": str(args.m6a), "msp": str(args.msp), "nuc": str(args.nuc)}
    # background parser per BED: gzip decompression overlaps with the scan
    rec_iters = {p: iter(ChromBlockStream(path)) for p, path in paths.items()}
    pending = {p: next(rec_iters[p], None) for p in paths}
    passed: Dict[str, set] = {p: set() for p in paths}
    # ONE shared namespace for all three files: classification compares m6a
    # fiber ids against msp/nuc fiber ids, so ids must be globally consistent
    mappers = {"all": IdMapper()}
    mapper_of = {"m6a": "all", "msp": "all", "nuc": "all"}

    sum_cov = np.zeros(width, dtype=np.float64)
    sum_m6a = np.zeros(width, dtype=np.float64)
    sum_acc = np.zeros(width, dtype=np.float64)
    sum_occ = np.zeros(width, dtype=np.float64)
    used = 0
    SITE_INCR = 100_000
    next_report = SITE_INCR

    args.outdir.mkdir(parents=True, exist_ok=True)
    scores_path = args.outdir / f"{args.sample}.motif_footprint_scores.tsv"
    scores_tmp = Path(f"{scores_path}.tmp")
    scores_fh = scores_tmp.open("w", encoding="utf-8")
    scores_fh.write("\t".join(
        ["chrom", "start", "end", "motif", "footprint_score", "mean_m6a_frac"]) + "\n")

    for chrom, idx_sites in by_chrom.items():
        # window extent for this chromosome: only [0, size) is ever read
        size = 0
        for s, e, _n in idx_sites:
            we = ((s + e) // 2) + flank
            if we > size:
                size = we
        # per-fiber block collection; EVERY file gets a per-fiber union merge
        # (legacy block_signal uniq_fibers=True) before building coverage
        file_parts = {"m6a": [], "msp": [], "nuc": []}

        n_rec = 0
        for p in ("m6a", "msp", "nuc"):
            item = pending[p]
            while item is not None and item[0] != chrom:
                passed[p].add(item[0])
                item = next(rec_iters[p], None)  # skip chromosomes without sites
            pending[p] = item
            if chrom in passed[p]:
                raise ValueError(
                    f"{paths[p]}: chromosome {chrom} records were already passed; "
                    f"site BED and element BED chromosome orders diverge "
                    f"(run `sort -k1,1` on every input)")
            if item is None:
                continue
            _chrom, b_fid, b_start, b_end, new_pairs = item
            mp_key = mapper_of[p]
            mappers[mp_key].extend(p, new_pairs)
            fid_arr = mappers[mp_key].arr(p)[b_fid]
            clip = np.minimum(b_end, size)
            keep = (b_start < clip) & (b_end > b_start)
            n_rec += int(keep.sum())
            file_parts[p].append((fid_arr[keep], b_start[keep], clip[keep]))
            pending[p] = next(rec_iters[p], None)

        # per-fiber union merge for EVERY file (legacy block_signal
        # uniq_fibers=True), then +1/-1 difference arrays -- identical to the
        # legacy coverage_signal()/block_signal() semantics.
        def merged_diff(parts):
            n_blocks = sum(f.size for f, _s, _e in parts)
            d = np.zeros(size + 1, dtype=np.int64)
            if not n_blocks:
                return d
            fid_a = np.concatenate([f for f, _s, _e in parts])
            st_a = np.concatenate([s for _f, s, _e in parts])
            en_a = np.concatenate([e for _f, _s, e in parts])
            order = np.lexsort((st_a, fid_a))
            fid_s = fid_a[order]
            st_s = st_a[order]
            en_s = en_a[order]
            first = np.empty(n_blocks, dtype=bool)
            first[0] = True
            np.not_equal(fid_s[1:], fid_s[:-1], out=first[1:])
            # segmented running max of block ends within each fiber: a
            # per-fiber base offset makes the global cummax segment-safe
            base = fid_s * (int(en_s.max()) + 1)
            seg_max = np.maximum.accumulate(en_s + base) - base
            prev_max = np.empty(n_blocks, dtype=np.int64)
            prev_max[0] = np.iinfo(np.int64).min
            prev_max[1:] = seg_max[:-1]
            starts = first | (st_s > prev_max)
            seg_idx = np.flatnonzero(starts)
            seg_starts = st_s[seg_idx]
            seg_ends = np.empty(seg_idx.shape[0], dtype=np.int64)
            seg_ends[:-1] = seg_max[seg_idx[1:] - 1]
            seg_ends[-1] = seg_max[-1]
            d = (np.bincount(seg_starts, minlength=size + 1)
                 - np.bincount(seg_ends, minlength=size + 1)).astype(np.int64)
            return d

        # nuc+msp cross-file union for fiber coverage
        un_parts = file_parts["msp"] + file_parts["nuc"]
        cov_m6a = np.cumsum(merged_diff(file_parts["m6a"])[:size])
        cov_msp = np.cumsum(merged_diff(file_parts["msp"])[:size])
        cov_nuc = np.cumsum(merged_diff(file_parts["nuc"])[:size])
        cov_un = np.cumsum(merged_diff(un_parts)[:size])
        del file_parts, un_parts

        for s, e, name in idx_sites:
            center = (s + e) // 2
            ws, we = center - flank, center + flank
            pad_l = max(0, -ws)
            pad_r = max(0, we - size)
            lo = max(0, ws)
            hi = min(size, we)
            if hi - lo == width:
                w_m6a = cov_m6a[lo:hi]
                w_msp = cov_msp[lo:hi]
                w_nuc = cov_nuc[lo:hi]
                w_un = cov_un[lo:hi]
            else:  # window clipped at contig bounds: pad with zero coverage
                w_m6a = np.pad(cov_m6a[lo:hi], (pad_l, pad_r))
                w_msp = np.pad(cov_msp[lo:hi], (pad_l, pad_r))
                w_nuc = np.pad(cov_nuc[lo:hi], (pad_l, pad_r))
                w_un = np.pad(cov_un[lo:hi], (pad_l, pad_r))

            cov = w_un.astype(np.float64)
            m6a_sig = w_m6a.astype(np.float64)
            acc = w_msp.astype(np.float64)
            occ = w_nuc.astype(np.float64)
            if cov.max() < args.min_fibers:
                continue
            used += 1
            sum_cov += cov
            sum_m6a += m6a_sig
            sum_acc += acc
            sum_occ += occ
            with np.errstate(divide="ignore", invalid="ignore"):
                m6a_frac = np.where(cov > 0, m6a_sig / np.maximum(cov, 1), 0.0)
            half = max(1, (e - s) // 2)
            scores_fh.write(
                f"{chrom}\t{s}\t{e}\t{name}\t"
                f"{round(footprint_score(m6a_frac, flank, half), 6)}\t"
                f"{round(float(np.mean(m6a_frac)), 6)}\n")
            if used >= next_report:
                el = time.time() - t0
                rate = used / el
                eta = (n_sites - used) / rate if rate else 0.0
                logger.info(
                    f"profiled {used}/{n_sites} sites ({used / n_sites * 100:.0f}%), "
                    f"{el:.0f}s elapsed, {rate:.0f} sites/s, ETA {eta / 60:.0f} min")
                next_report += SITE_INCR
        logger.info(f"chrom {chrom}: {len(idx_sites)} sites, {n_rec} fiber records, "
                    f"{time.time() - t0:.0f}s elapsed")

    scores_fh.close()
    os.replace(scores_tmp, scores_path)

    logger.info(f"Used {used}/{n_sites} sites (min_fibers={args.min_fibers})")
    if used == 0:
        logger.error("No motif site passed min_fibers filter; no outputs written")
        return 1

    with np.errstate(divide="ignore", invalid="ignore"):
        m6a_frac = np.where(sum_cov > 0, sum_m6a / np.maximum(sum_cov, 1), np.nan)
        acc_frac = np.where(sum_cov > 0, sum_acc / np.maximum(sum_cov, 1), np.nan)
        occ_frac = np.where(sum_cov > 0, sum_occ / np.maximum(sum_cov, 1), np.nan)

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

    _plot(args, rel, m6a_frac, acc_frac, occ_frac)
    logger.info(f"Wrote {profile_tsv}, {scores_path} and profile PDF "
                f"(total {time.time() - t0:.0f}s)")
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
    ax.set_xlabel("position relative to motif center (bp)")
    ax.set_ylabel("fraction of fibers")
    ax.set_title(f"{args.sample} motif-centered profile")
    ax.legend(frameon=False, fontsize=8)
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
                   help="motif occurrence BED")
    p.add_argument("-a", "--m6a", type=Path, required=True, help="m6A BED12 (.bed.gz)")
    p.add_argument("-s", "--msp", type=Path, required=True, help="MSP BED12 (.bed.gz)")
    p.add_argument("-u", "--nuc", type=Path, required=True, help="nucleosome BED12 (.bed.gz)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--sample", type=str, default="sample", help="sample name prefix")
    p.add_argument("-f", "--flank", type=int, default=250,
                   help="window half-size around motif center (bp)")
    p.add_argument("-i", "--min-fibers", type=int, default=10,
                   help="minimum fiber coverage per site to include")
    p.add_argument("--motif-filter", type=str, default=None,
                   help="comma-separated motif names to keep (part before '@')")
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
