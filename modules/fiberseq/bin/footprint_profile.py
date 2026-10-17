#!/usr/bin/env python3
"""TF footprinting on single fibers around motif sites.

Classifies every fiber over a motif site into three categories [1]:
  (i)   MSP covering the motif with an m6A-free footprint in the motif core,
  (ii)  MSP covering the motif without a footprint,
  (iii) nucleosome occupying the motif.
Then aggregates m6A fraction profiles per category and compares per-site
footprint scores between categories with a Kolmogorov-Smirnov test.

Implementation note (vectorized):
    The legacy implementation re-scanned all fiber records overlapping each
    site window separately for every site (hundreds of ms/site => weeks on
    full motif catalogs). This version streams each element BED exactly once
    per chromosome into per-fiber MERGED span lists (a fiber's disjoint
    BED12 blocks are union-merged per file, plus the msp+nuc cross-file
    union), position-sorted. Each site then

    * classifies fibers with searchsorted slices + ``np.isin`` (no per-fiber
      Python loop) -- a merged span overlaps the core iff some raw block of
      that fiber does (the union is contiguous), so classification is exact;
    * derives per-category window signals by clipping merged spans of
      category fibers into a window difference array -- exactly the legacy
      block_signal(uniq_fibers=True) / coverage_signal() semantics
      (positions before the contig start carry zero coverage, as in the
      legacy block clipping).

    Site order, arithmetic and rounding are unchanged, so outputs are
    identical to the legacy run for the same chromosome-grouped inputs.

Usage:
    python footprint_profile.py -m ctcf_sites.bed -a sample.m6a.bed.gz \
        -s sample.msp.bed.gz -u sample.nuc.bed.gz -o out/ --sample s1 \
        [--motif-filter CTCF@,MYC@]
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


def footprint_score(m6a_frac: np.ndarray, motif_half: int) -> float:
    """Flank-max minus center-sum score over the motif window."""
    n = len(m6a_frac)
    center = n // 2
    left = m6a_frac[: center - motif_half]
    right = m6a_frac[center + motif_half:]
    center_sig = m6a_frac[center - motif_half: center + motif_half]
    if len(left) == 0 or len(right) == 0:
        return 0.0
    return float(np.max(left) + np.max(right) - np.sum(center_sig))


def _merged_spans(fid: array, start: array, end: array,
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Union-merge blocks per fiber; return (fid, start, end) sorted by start.

    Same semantics as FiberlibUtil.block_signal's per-fiber span merge,
    vectorized via a segmented running-max.
    """
    n = len(fid)
    if n == 0:
        e_f = np.empty(0, dtype=np.int64)
        return e_f, e_f.copy(), e_f.copy()
    fid_a = np.frombuffer(fid, dtype=np.int64)
    st_a = np.frombuffer(start, dtype=np.int64)
    en_a = np.frombuffer(end, dtype=np.int64)
    order = np.lexsort((st_a, fid_a))
    f = fid_a[order]
    s = st_a[order]
    e = en_a[order]
    first = np.empty(n, dtype=bool)
    first[0] = True
    np.not_equal(f[1:], f[:-1], out=first[1:])
    # segmented running max of block ends within each fiber: a per-fiber
    # base offset makes the global cummax segment-safe
    base = f * (int(e.max()) + 1)
    run_max = np.maximum.accumulate(e + base) - base
    prev_max = np.empty(n, dtype=np.int64)
    prev_max[0] = np.iinfo(np.int64).min
    prev_max[1:] = run_max[:-1]
    starts = first | (s > prev_max)
    idx = np.flatnonzero(starts)
    m_start = s[idx]
    m_end = np.empty(idx.shape[0], dtype=np.int64)
    m_end[:-1] = run_max[idx[1:] - 1]
    m_end[-1] = run_max[-1]
    m_fid = f[idx]
    pos = np.argsort(m_start, kind="stable")
    return m_fid[pos], m_start[pos], m_end[pos]


class _SpanSet:
    """Position-sorted per-fiber merged spans of one element set."""

    __slots__ = ("fid", "start", "end", "max_span")

    def __init__(self, fid: np.ndarray, start: np.ndarray, end: np.ndarray):
        self.fid = fid
        self.start = start
        self.end = end
        self.max_span = int(np.max(end - start)) if len(end) else 1

    def window(self, ws: int, we: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Merged spans overlapping [ws, we), clipped to the window (rel. coords)."""
        lo = np.searchsorted(self.start, ws - self.max_span, side="left")
        hi = np.searchsorted(self.start, we, side="right")
        s = self.start[lo:hi]
        e = self.end[lo:hi]
        f = self.fid[lo:hi]
        keep = (s < we) & (e > ws)
        s = np.clip(s[keep], ws, we) - ws
        e = np.clip(e[keep], ws, we) - ws
        return f[keep], s, e


def _coverage_window(spans: _SpanSet, keep_fids: Optional[np.ndarray],
                     ws: int, we: int) -> np.ndarray:
    """Per-base distinct-fiber coverage over [ws, we) from merged spans.

    keep_fids: sorted fiber ids to restrict to (None = all). Positions
    outside the contig carry zero coverage, matching legacy block clipping.
    """
    width = we - ws
    diff = np.zeros(width + 1, dtype=np.int64)
    f, s, e = spans.window(ws, we)
    if keep_fids is not None and len(f):
        pos = np.searchsorted(keep_fids, f)
        pos[pos == len(keep_fids)] = 0
        m = keep_fids[pos] == f
        f, s, e = f[m], s[m], e[m]
    if len(s):
        diff += np.bincount(s, minlength=width + 1)
        diff -= np.bincount(e, minlength=width + 1)
    return np.cumsum(diff[:width])


def _overlapping_fids(spans: _SpanSet, cs: int, ce: int) -> np.ndarray:
    """Sorted unique fiber ids with any merged span overlapping [cs, ce)."""
    lo = np.searchsorted(spans.start, cs - spans.max_span, side="left")
    hi = np.searchsorted(spans.start, ce, side="right")
    s = spans.start[lo:hi]
    e = spans.end[lo:hi]
    f = spans.fid[lo:hi]
    hit = f[(s < ce) & (e > cs)]
    return np.unique(hit)


def _load_sites_grouped(path: str, motif_filter: Optional[set]) -> Dict[str, List[Tuple[int, int, str]]]:
    """Stream the sites BED into {chrom: [(start, end, name), ...]} (file order)."""
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
    from scipy.stats import ks_2samp

    t0 = time.time()
    motif_filter = None
    if args.motif_filter:
        motif_filter = {k.strip() for k in args.motif_filter.split(",") if k.strip()}
    by_chrom = _load_sites_grouped(str(args.motif_sites), motif_filter)
    n_sites = sum(len(v) for v in by_chrom.values())
    logger.info(f"Loaded {n_sites} motif sites")
    if n_sites == 0:
        logger.error("no motif sites; no outputs")
        return 1
    if n_sites > 1_000_000 and not motif_filter:
        logger.warning(
            f"{n_sites} sites and no --motif-filter: outputs will be huge; "
            f"consider --motif-filter")

    flank = args.flank
    width = 2 * flank

    paths = {"m6a": str(args.m6a), "msp": str(args.msp), "nuc": str(args.nuc)}
    # background parser per BED: gzip decompression overlaps with the scan
    rec_iters = {p: iter(ChromBlockStream(path)) for p, path in paths.items()}
    pending = {p: next(rec_iters[p], None) for p in paths}
    passed: Dict[str, set] = {p: set() for p in paths}
    # ONE shared namespace for all three files: classification compares m6a
    # fiber ids against msp/nuc fiber ids, so ids must be globally consistent
    mappers = {"all": IdMapper()}
    mapper_of = {"m6a": "all", "msp": "all", "nuc": "all"}

    args.outdir.mkdir(parents=True, exist_ok=True)
    cats_path = args.outdir / f"{args.sample}.footprint_categories.tsv"
    cats_fh = Path(f"{cats_path}.tmp").open("w", encoding="utf-8")
    cats_fh.write("\t".join(["chrom", "start", "end", "motif", *CATEGORIES]) + "\n")

    agg = {c: {"cov": np.zeros(width), "m6a": np.zeros(width)} for c in CATEGORIES}
    score_sets: Dict[str, array] = {c: array("d") for c in CATEGORIES}

    SITE_INCR = 100_000
    next_report = SITE_INCR
    done = 0

    for chrom, idx_sites in by_chrom.items():
        size = 0
        for s, e, _n in idx_sites:
            we = ((s + e) // 2) + flank
            if we > size:
                size = we

        raw: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for p in paths:
            item = pending[p]
            while item is not None and item[0] != chrom:
                passed[p].add(item[0])
                item = next(rec_iters[p], None)
            pending[p] = item
            if chrom in passed[p]:
                raise ValueError(
                    f"{paths[p]}: chromosome {chrom} records were already passed; "
                    f"site BED and element BED chromosome orders diverge "
                    f"(run `sort -k1,1` on every input)")
            if item is None:
                raw[p] = (np.empty(0, dtype=np.int64),) * 3
                continue
            _chrom, b_fid, b_start, b_end, new_pairs = item
            mp_key = mapper_of[p]
            mappers[mp_key].extend(p, new_pairs)
            fid_arr = mappers[mp_key].arr(p)[b_fid]
            clip = np.minimum(b_end, size)
            keep = (b_start < clip) & (b_end > b_start)
            raw[p] = (fid_arr[keep], b_start[keep], clip[keep])
            pending[p] = next(rec_iters[p], None)

        S_m6a = _SpanSet(*_merged_spans(*raw["m6a"]))
        S_msp = _SpanSet(*_merged_spans(*raw["msp"]))
        S_nuc = _SpanSet(*_merged_spans(*raw["nuc"]))
        # a fiber's MSP and nucleosome spans may overlap: union-merge per fiber
        # (legacy coverage_signal counts such a fiber once per base)
        S_un = _SpanSet(*_merged_spans(
            np.concatenate([S_msp.fid, S_nuc.fid]),
            np.concatenate([S_msp.start, S_nuc.start]),
            np.concatenate([S_msp.end, S_nuc.end])))
        del raw

        for s, e, name in idx_sites:
            core_s, core_e = s, e
            center = (s + e) // 2
            ws, we = center - flank, center + flank

            msp_fids = _overlapping_fids(S_msp, core_s, core_e)
            nuc_fids = _overlapping_fids(S_nuc, core_s, core_e)
            m6a_fids = _overlapping_fids(S_m6a, core_s, core_e)

            if len(msp_fids) and len(m6a_fids):
                in_m6a = np.isin(msp_fids, m6a_fids, assume_unique=True)
                fp_fids = msp_fids[~in_m6a]
                nf_fids = msp_fids[in_m6a]
            else:
                fp_fids = msp_fids
                nf_fids = msp_fids[:0]
            if len(nuc_fids) and len(msp_fids):
                nuc_cat = nuc_fids[~np.isin(nuc_fids, msp_fids, assume_unique=True)]
            else:
                nuc_cat = nuc_fids

            cats_fh.write(f"{chrom}\t{s}\t{e}\t{name}\t"
                          f"{len(fp_fids)}\t{len(nf_fids)}\t{len(nuc_cat)}\n")
            done += 1
            if done >= next_report:
                el = time.time() - t0
                rate = done / el
                logger.info(
                    f"processed {done}/{n_sites} sites ({done / n_sites * 100:.0f}%), "
                    f"{el:.0f}s elapsed, {rate:.0f} sites/s, "
                    f"ETA {(n_sites - done) / rate / 60:.0f} min")
                next_report += SITE_INCR

            for cat, keep in (("msp_footprint", fp_fids),
                              ("msp_no_footprint", nf_fids),
                              ("nucleosome", nuc_cat)):
                if not keep.size:
                    continue
                cov = _coverage_window(S_un, keep, ws, we).astype(np.float64)
                m6a_sig = _coverage_window(S_m6a, keep, ws, we).astype(np.float64)
                agg[cat]["cov"] += cov
                agg[cat]["m6a"] += m6a_sig
                with np.errstate(divide="ignore", invalid="ignore"):
                    frac = np.where(cov > 0, m6a_sig / np.maximum(cov, 1), 0.0)
                score_sets[cat].append(
                    footprint_score(frac, max(1, (e - s) // 2)))

        logger.info(f"chrom {chrom}: {len(idx_sites)} sites, "
                    f"{len(S_un.start)} union spans, {time.time() - t0:.0f}s elapsed")

    cats_fh.close()
    os.replace(Path(f"{cats_path}.tmp"), cats_path)

    if sum(len(v) for v in score_sets.values()) == 0:
        logger.error("no fiber classified at any site; no outputs")
        return 1

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
        a = np.asarray(score_sets[c1], dtype=np.float64)
        b = np.asarray(score_sets[c2], dtype=np.float64)
        if a.size and b.size:
            stat, p = ks_2samp(a, b)
            score_rows.append([c1, c2, int(a.size), int(b.size),
                               round(float(stat), 6), float(p)])
        else:
            score_rows.append([c1, c2, int(a.size), int(b.size), "NA", "NA"])
    write_tsv(str(args.outdir / f"{args.sample}.footprint_scores.tsv"),
              ["group1", "group2", "n1", "n2", "ks_stat", "ks_pvalue"], score_rows)

    _plot(args, rel, agg)
    logger.info(f"Wrote footprint tables and plot to {args.outdir} "
                f"(total {time.time() - t0:.0f}s)")
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
    p.add_argument("--motif-filter", type=str, default=None,
                   help="comma-separated motif names to keep (part before '@')")
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
