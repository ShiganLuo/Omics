#!/usr/bin/env python3
"""Cluster single chromatin fibers by their accessibility pattern at target regions.

For each region (e.g. a marker-gene TSS window), extract every Fiber-seq fiber
(ft extract BED12, one row per read) overlapping the region from one or more
samples, binarize the per-bin accessibility (MSP) pattern, cluster the fibers
(hierarchical or k-means), and report per-cluster patterns plus per-sample
cluster proportions.

Clustering unit = ONE fiber = one chromatin molecule from ONE nucleus, so a
cluster proportion approximates the fraction of (alleles of) cells carrying
that chromatin state. Bulk Fiber-seq carries no cell barcode: clusters are
chromatin-state subpopulations, NOT identified cell types.

BED12 semantics (fibertools-rs v0.13.1): the FIRST and LAST blocks of every
record are boundary sentinels and are dropped; only interior blocks are real
element calls (reference implementation: modules/fibertools/bin/distributed_plot.py).
m6A/MSP is the enzyme-deposited ACCESSIBILITY mark, not a native modification.

Usage:
    python fiber_cluster.py \
        -m luanchao-21310=/data/.../4_elements/luanchao-21310/luanchao-21310.msp.bed.gz \
        -m luanchao-11236=/data/.../4_elements/luanchao-11236/luanchao-11236.msp.bed.gz \
        -u luanchao-21310=/data/.../luanchao-21310.nuc.bed.gz \
        -u luanchao-11236=/data/.../luanchao-11236.nuc.bed.gz \
        -f luanchao-21310=/data/.../luanchao-21310.fire_peaks.bed \
        -f luanchao-11236=/data/.../luanchao-11236.fire_peaks.bed \
        -r regions.bed -o fiber_cluster_out/ \
        -b 100 -k 3 -e hierarchical -s luanchao-21310,luanchao-11236

    # same via inline loci instead of a regions BED (name defaults to coordinates)
    python fiber_cluster.py -m S1=... -l 'NC_133418.1:75646761-75656761|FSHR' ...
"""
from __future__ import annotations

import argparse
import gzip
import logging
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# setup_logger copied verbatim from src/common/util/LogUtil.py (module bin
# scripts are self-contained; provenance: Omics repo LogUtil, 2026-11).
# ---------------------------------------------------------------------------
from typing import Literal, TextIO  # noqa: E402

PlotFormat = Literal["png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf", "raw", "rgba"]


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


logger = setup_logger("FiberCluster", stream=sys.stdout)


Region = Tuple[str, str, int, int]  # (name, chrom, start, end)
Blocks = List[Tuple[int, int]]  # absolute [start, end) intervals
FiberRec = Tuple[str, Blocks, int, int]  # (read_name, element_blocks, fiber_start, fiber_end)
FiberChrom = Dict[str, List[FiberRec]]  # chrom -> fiber records
FiberSample = Dict[str, FiberChrom]  # sample -> chrom -> fiber records
FireSample = Dict[str, Dict[str, Blocks]]  # sample -> chrom -> merged significant FIRE peaks


class FiberRow(NamedTuple):
    """One clustering unit: a fiber passing the span filter at one region."""
    read: str
    sample: str
    span_frac: float
    open: np.ndarray


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_sample_paths(entries: Sequence[str], flag: str) -> Dict[str, Path]:
    """Parse repeated '<sample>=<path>' CLI entries into a sample->path map.

    Parameters:
        entries: raw CLI strings in '<sample>=<path>' form.
        flag: flag name used in error messages.

    Returns:
        Mapping of sample id to validated existing path.

    Errors:
        ValueError on malformed entry or missing file.
    """
    out: Dict[str, Path] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"{flag}: expected '<sample>=<path>', got {entry!r}")
        sample, path_str = entry.split("=", 1)
        if not sample:
            raise ValueError(f"{flag}: empty sample id in {entry!r}")
        path = Path(path_str)
        if not path.is_file():
            raise ValueError(f"{flag}: file not found for sample {sample!r}: {path}")
        if sample in out:
            raise ValueError(f"{flag}: duplicate sample {sample!r}")
        out[sample] = path
    return out


def parse_loci(loci: Sequence[str]) -> List[Region]:
    """Parse '<chrom>:<start>-<end>|<name>' strings into regions.

    Parameters:
        loci: raw CLI locus strings; the '|name' part is optional.

    Returns:
        List of (name, chrom, start, end) regions.

    Errors:
        ValueError on unparseable strings.
    """
    regions: List[Region] = []
    for loc in loci:
        name = ""
        coord = loc
        if "|" in loc:
            coord, name = loc.split("|", 1)
        m = re.fullmatch(r"([^:]+):(\d+)-(\d+)", coord)
        if not m:
            raise ValueError(f"--locus: unparseable {loc!r} (want 'chr:start-end|name')")
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        if end <= start:
            raise ValueError(f"--locus: end <= start in {loc!r}")
        regions.append((name or f"{chrom}_{start}_{end}", chrom, start, end))
    return regions


def read_regions_bed(path: Path) -> List[Region]:
    """Read a regions BED (>=4 columns; col4 = region name).

    Parameters:
        path: BED file path.

    Returns:
        List of (name, chrom, start, end) regions.
    """
    regions: List[Region] = []
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                raise ValueError(f"{path}:{line_no}: fewer than 3 columns")
            name = f[3] if len(f) > 3 and f[3] else f"{f[0]}_{f[1]}_{f[2]}"
            regions.append((name, f[0], int(f[1]), int(f[2])))
    return regions


# ---------------------------------------------------------------------------
# BED12 streaming extraction
# ---------------------------------------------------------------------------

def merged_intervals(regions: Sequence[Region]) -> Dict[str, List[Tuple[int, int]]]:
    """Per-chromosome sorted disjoint [start, end) intervals of all regions."""
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for _, chrom, start, end in regions:
        by_chrom[chrom].append((start, end))
    out: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, spans in by_chrom.items():
        spans.sort()
        merged: List[List[int]] = [list(spans[0])]
        for s, e in spans[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[chrom] = [(s, e) for s, e in merged]
    return out


def _overlaps_any(spans: Sequence[Tuple[int, int]], start: int, end: int) -> bool:
    """True if [start, end) intersects any interval in the sorted disjoint list."""
    lo, hi = 0, len(spans)
    while lo < hi:  # first span with span_end > start
        mid = (lo + hi) // 2
        if spans[mid][1] <= start:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(spans) and spans[lo][0] < end


def extract_bed12(path: Path, query: Dict[str, List[Tuple[int, int]]]) -> Dict[str, List[FiberRec]]:
    """Stream a gzipped ft-extract BED12 and keep region-overlapping fibers.

    One pass per file regardless of region count. Sentinel handling: the first
    and last block of each BED12 record are boundary sentinels (ft extract
    v0.13.1) and are dropped; interior blocks are real element calls.

    Parameters:
        path: bgzip/gzip BED12 element file (m6a/msp/nuc/cpg).
        query: per-chromosome merged query intervals.

    Returns:
        chrom -> list of (read_name, blocks, fiber_start, fiber_end); blocks
        are absolute [start, end) intervals clipped to the query intervals.
    """
    out: Dict[str, List[FiberRec]] = defaultdict(list)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] not in query:
                continue
            start, end = int(f[1]), int(f[2])
            spans = query[f[0]]
            if not _overlaps_any(spans, start, end):
                continue
            name = f[3]
            n_blocks = int(f[9])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            if len(sizes) != n_blocks or len(offsets) != n_blocks:
                logger.warning("malformed BED12 block lists at %s:%d-%d; record skipped", f[0], start, end)
                continue
            # Drop the two boundary sentinel blocks; keep interior calls only.
            blocks: Blocks = []
            for size, off in zip(sizes[1:-1], offsets[1:-1]):
                if size <= 0:
                    continue
                b_start, b_end = start + off, start + off + size
                if _overlaps_any(spans, b_start, b_end):
                    blocks.append((b_start, b_end))
            out[f[0]].append((name, blocks, start, end))
    return out


def extract_fire_peaks(path: Path, query: Dict[str, List[Tuple[int, int]]],
                       fdr_max: float) -> Dict[str, Blocks]:
    """Load significant FIRE peaks overlapping the query intervals.

    fire_peaks.bed columns: #chrom peak_start peak_end peak_max FDR coverage
    fire_coverage score ... (fibertools `ft fire --extract` peak table).
    Peaks with FDR <= fdr_max are kept and merged per chromosome. The
    per-fiber fire.bed.gz is NOT used here: it carries all MSP intervals, so
    significance can only be filtered at the peak level.

    Parameters:
        path: fire_peaks.bed (plain or gzipped).
        query: per-chromosome merged query intervals.
        fdr_max: keep peaks with FDR <= this value.

    Returns:
        chrom -> merged significant peak intervals [start, end).
    """
    raw: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if chrom not in query or not _overlaps_any(query[chrom], start, end):
                continue
            try:
                fdr = float(f[4])
            except ValueError:
                continue
            if fdr <= fdr_max:
                raw[chrom].append((start, end))

    out: Dict[str, Blocks] = {}
    for chrom, pts in raw.items():
        pts.sort()
        merged: List[List[int]] = [list(pts[0])]
        for s, e in pts[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[chrom] = [(s, e) for s, e in merged]
    return out


def intersect_blocks(blocks: Blocks, regions_b: Blocks) -> Blocks:
    """Intersect two sorted interval lists.

    Parameters:
        blocks: first interval list (absolute coords).
        regions_b: second interval list (absolute coords).

    Returns:
        Non-empty intersection intervals.
    """
    out: Blocks = []
    i = j = 0
    while i < len(blocks) and j < len(regions_b):
        lo = max(blocks[i][0], regions_b[j][0])
        hi = min(blocks[i][1], regions_b[j][1])
        if hi > lo:
            out.append((lo, hi))
        if blocks[i][1] <= regions_b[j][1]:
            i += 1
        else:
            j += 1
    return out


# ---------------------------------------------------------------------------
# Binarization and clustering
# ---------------------------------------------------------------------------

def bin_covered_bp(blocks: Blocks, start: int, end: int, bin_size: int) -> np.ndarray:
    """Per-bin covered base pairs of a fiber's element blocks.

    Parameters:
        blocks: absolute [start, end) element intervals of one fiber.
        start, end: region span.
        bin_size: bin width in bp.

    Returns:
        Covered bp per bin (bins fully outside the blocks are 0).
    """
    n_bins = max(1, math.ceil((end - start) / bin_size))
    cov = np.zeros(n_bins, dtype=float)
    for b_start, b_end in blocks:
        lo = max(b_start, start)
        hi = min(b_end, end)
        if hi <= lo:
            continue
        i_lo = (lo - start) // bin_size
        i_hi = min(n_bins - 1, (hi - 1 - start) // bin_size)
        for i in range(i_lo, i_hi + 1):
            bin_lo = start + i * bin_size
            bin_hi = min(end, bin_lo + bin_size)
            cov[i] += min(hi, bin_hi) - max(lo, bin_lo)
    return cov


def binarize(cov_bp: np.ndarray, start: int, end: int, bin_size: int,
             mode: str, bin_threshold: float, min_bp: int) -> np.ndarray:
    """Convert per-bin covered bp into a binary open/closed vector.

    Parameters:
        cov_bp: covered bp per bin (from bin_covered_bp).
        start, end: region span.
        bin_size: bin width in bp.
        mode: 'any' -> bin open when covered bp >= min_bp; 'frac' -> bin open
            when covered fraction >= bin_threshold.
        bin_threshold: fraction cutoff for mode 'frac'.
        min_bp: bp cutoff for mode 'any'.

    Returns:
        Binary (0/1) vector of length n_bins.
    """
    n_bins = len(cov_bp)
    widths = np.array([min(end, start + (i + 1) * bin_size) - (start + i * bin_size)
                       for i in range(n_bins)], dtype=float)
    if mode == "any":
        return (cov_bp >= min_bp).astype(int)
    return (cov_bp / widths >= bin_threshold).astype(int)


def cluster_matrix(mat: np.ndarray, k: int, method: str, seed: int) -> np.ndarray:
    """Cluster binary fiber x bin patterns.

    Parameters:
        mat: (n_fibers, n_bins) binary accessibility matrix.
        k: requested cluster count.
        method: 'hierarchical' (average-linkage, Hamming) or 'kmeans'.
        seed: RNG seed for k-means init.

    Returns:
        1-based cluster labels per fiber, relabeled so cluster 1 has the
        highest mean openness (deterministic label semantics).

    Errors:
        ValueError when clustering cannot produce k non-empty clusters.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist

    n = mat.shape[0]
    if n < k:
        raise ValueError(f"only {n} fibers for k={k}")
    if method == "hierarchical":
        dist = pdist(mat.astype(float), metric="hamming")
        z = linkage(dist, method="average")
        labels = fcluster(z, t=k, criterion="maxclust")
    elif method == "kmeans":
        from scipy.cluster.vq import kmeans2

        np.random.seed(seed)
        _, labels = kmeans2(mat.astype(float), k, minit="++", iter=50)
        labels = labels + 1
    else:
        raise ValueError(f"unknown method {method!r}")

    # Relabel by mean openness (descending) so labels are stable across runs.
    order = sorted(set(int(x) for x in labels),
                   key=lambda c: -float(mat[labels == c].mean()))
    remap = {old: new for new, old in enumerate(order, start=1)}
    out = np.array([remap[int(x)] for x in labels], dtype=int)
    if len(set(out.tolist())) < k:
        logger.warning("clustering produced %d non-empty clusters, expected %d "
                       "(identical/tied patterns in this channel)", len(set(out.tolist())), k)
    return out


# ---------------------------------------------------------------------------
# Region-level pipeline
# ---------------------------------------------------------------------------

def run_region(region: Region, msp: FiberSample, nuc: FiberSample,
               fire: FireSample, sample_order: List[str], args: argparse.Namespace,
               outdir: Path) -> None:
    """Build the fiber matrix, cluster, and write tables + figure for one region.

    Parameters:
        region: (name, chrom, start, end).
        msp, nuc: per-sample chrom -> fiber BED12 records.
        fire: per-sample chrom -> merged significant FIRE peaks.
        sample_order: display order of samples.
        args: CLI namespace (bin_size, k, method, thresholds, seed).
        outdir: per-region output directory.

    Outputs (inside outdir):
        fiber_matrix.tsv, cluster_assignments.tsv, cluster_patterns.tsv,
        cluster_summary.tsv, <region>.clusters.png
    """
    name, chrom, start, end = region
    bin_size = args.bin_size
    n_bins = max(1, math.ceil((end - start) / bin_size))
    reg_len = end - start

    # ---- collect fibers per sample (MSP is the clustering channel) --------
    rows: List[FiberRow] = []
    fiber_blocks: List[Blocks] = []  # MSP blocks per kept fiber (FIRE overlay source)
    for sample in sample_order:
        for read_name, blocks, f_start, f_end in msp.get(sample, {}).get(chrom, []):
            overlap = min(f_end, end) - max(f_start, start)
            if overlap <= 0:
                continue
            span_frac = overlap / reg_len
            if span_frac < args.min_span_frac:
                continue
            if args.cluster_channel == "fire":
                peaks_b = fire.get(sample, {}).get(chrom, [])
                channel_blocks = intersect_blocks(sorted(blocks), peaks_b)
            else:
                channel_blocks = blocks
            cov_bp = bin_covered_bp(channel_blocks, start, end, bin_size)
            rows.append(FiberRow(read_name, sample, span_frac,
                                 binarize(cov_bp, start, end, bin_size, args.bin_mode,
                                          args.bin_threshold, args.min_msp_bp)))
            fiber_blocks.append(blocks)
    n_fibers = len(rows)
    if n_fibers < args.min_fibers:
        raise ValueError(f"too_few_fibers: {n_fibers} < --min-fibers {args.min_fibers}")

    open_mat = np.vstack([r.open for r in rows]).astype(int)
    labels = cluster_matrix(open_mat, args.clusters, args.method, args.seed)

    # ---- overlay matrices (nuc occupancy, FIRE presence) -----------------
    def nuc_matrix(recs: FiberSample) -> np.ndarray:
        mat = np.zeros((n_fibers, n_bins), dtype=float)
        for i, r in enumerate(rows):
            hits: Blocks = []
            for r_name, blocks, _f_start, _f_end in recs.get(r.sample, {}).get(chrom, []):
                if r_name == r.read:
                    hits.extend(blocks)
            cov_bp = bin_covered_bp(sorted(hits), start, end, bin_size)
            mat[i] = binarize(cov_bp, start, end, bin_size, "frac",
                              args.bin_threshold, args.min_msp_bp)
        return mat

    def fire_matrix(peaks: FireSample) -> np.ndarray:
        # Per-fiber FIRE = this fiber's MSP blocks intersected with the
        # sample's significant FIRE peaks (peak-level FDR filtering).
        mat = np.zeros((n_fibers, n_bins), dtype=float)
        for i, r in enumerate(rows):
            pk = peaks.get(r.sample, {}).get(chrom, [])
            if not pk:
                continue
            inter = intersect_blocks(sorted(fiber_blocks[i]), pk)
            cov_bp = bin_covered_bp(inter, start, end, bin_size)
            mat[i] = binarize(cov_bp, start, end, bin_size, args.bin_mode,
                              args.bin_threshold, args.min_msp_bp)
        return mat

    nuc_mat = nuc_matrix(nuc) if nuc else np.zeros((n_fibers, n_bins))
    fire_mat = fire_matrix(fire) if fire else np.zeros((n_fibers, n_bins))

    # ---- tables ----------------------------------------------------------
    outdir.mkdir(parents=True, exist_ok=True)
    bin_starts = [start + i * bin_size for i in range(n_bins)]

    with open(outdir / "cluster_assignments.tsv", "w") as fh:
        fh.write("region\tread_name\tsample\tcluster\tspan_frac\tn_open_bins\n")
        for r, lab in zip(rows, labels):
            fh.write(f"{name}\t{r.read}\t{r.sample}\t{lab}\t{r.span_frac:.3f}\t{int(r.open.sum())}\n")

    with open(outdir / "cluster_patterns.tsv", "w") as fh:
        fh.write("region\tcluster\tbin_start\tbin_end\topen_frac\tnuc_frac\tfire_frac\tn_fibers\n")
        for c in sorted(set(labels.tolist())):
            sel = labels == c
            o = open_mat[sel].mean(axis=0)
            nu = nuc_mat[sel].mean(axis=0)
            fi = fire_mat[sel].mean(axis=0)
            for i in range(n_bins):
                fh.write(f"{name}\t{c}\t{bin_starts[i]}\t{min(end, bin_starts[i] + bin_size)}\t"
                         f"{o[i]:.3f}\t{nu[i]:.3f}\t{fi[i]:.3f}\t{int(sel.sum())}\n")

    with open(outdir / "cluster_summary.tsv", "w") as fh:
        fh.write("region\tcluster\tsample\tn_fibers\tfrac_of_sample\n")
        for sample in sample_order:
            n_s = sum(1 for r in rows if r.sample == sample)
            for c in sorted(set(labels.tolist())):
                n_c = sum(1 for r, lab in zip(rows, labels) if lab == c and r.sample == sample)
                frac = n_c / n_s if n_s else 0.0
                fh.write(f"{name}\t{c}\t{sample}\t{n_c}\t{frac:.4f}\n")

    with open(outdir / "fiber_matrix.tsv", "w") as fh:
        fh.write("read_name\tsample\tcluster\tspan_frac\t" +
                 "\t".join(f"bin_{s}" for s in bin_starts) + "\n")
        for r, lab in zip(rows, labels):
            fh.write(f"{r.read}\t{r.sample}\t{lab}\t{r.span_frac:.3f}\t" +
                     "\t".join(str(x) for x in r.open) + "\n")

    # ---- figure ----------------------------------------------------------
    plot_region(region, rows, labels, open_mat, nuc_mat, fire_mat, sample_order, args,
                outdir / f"{_safe(name)}.clusters.png")


def _safe(name: str) -> str:
    """Sanitize a region name for use in filenames."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def plot_region(region: Region, rows: List[FiberRow], labels: np.ndarray,
                open_mat: np.ndarray, nuc_mat: np.ndarray, fire_mat: np.ndarray,
                sample_order: List[str], args: argparse.Namespace, out_png: Path) -> None:
    """Render the per-region cluster figure.

    Layout: fiber x bin accessibility heatmap (rows grouped by cluster, sample
    side strip), per-bin aggregate tracks (open/nuc/fire), and per-cluster
    per-sample composition bars. Youth/reference samples stack on top
    (--sample-order).

    Parameters:
        region: (name, chrom, start, end).
        rows: per-fiber clustering rows.
        labels: cluster labels per fiber.
        open_mat, nuc_mat, fire_mat: binary fiber x bin matrices.
        sample_order: display order of samples.
        args: CLI namespace (bin_size, k, method).
        out_png: output image path.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Rectangle

    name, chrom, start, end = region
    n_bins = open_mat.shape[1]
    sample_rank = {s: i for i, s in enumerate(sample_order)}

    # Sort rows: cluster, then sample (top = first in --sample-order), then openness.
    order = sorted(range(len(rows)),
                   key=lambda i: (int(labels[i]), sample_rank.get(rows[i].sample, 99),
                                  -int(rows[i].open.sum())))
    labels_s = labels[order]
    rows_s = [rows[i] for i in order]
    open_s = open_mat[order]

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[4.2, 1.3], height_ratios=[4.2, 1.0],
                          hspace=0.28, wspace=0.12)
    ax_hm = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_tr = fig.add_subplot(gs[1, 0])

    cmap = ListedColormap(["#f2f2f2", "#b2182b"])
    ax_hm.imshow(open_s, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    bounds = np.cumsum([0] + [int((labels_s == c).sum()) for c in sorted(set(labels_s.tolist()))])
    for b in bounds[1:-1]:
        ax_hm.axhline(b - 0.5, color="black", lw=1.0)
    ax_hm.set_yticks([bounds[i] + (bounds[i + 1] - bounds[i]) / 2 for i in range(len(bounds) - 1)])
    ax_hm.set_yticklabels([f"cluster {c}\n(n={int((labels_s == c).sum())})"
                           for c in sorted(set(labels_s.tolist()))])
    # sample side strip
    for i, r in enumerate(rows_s):
        ax_hm.add_patch(Rectangle((-0.06 * n_bins - 1, i - 0.5), 0.05 * n_bins, 1.0,
                                  color=plt.cm.tab10(sample_rank.get(r.sample, 9) % 10),
                                  clip_on=False, lw=0))
    tick_pos = np.linspace(0, n_bins - 1, 6).astype(int)
    ax_hm.set_xticks(tick_pos)
    ax_hm.set_xticklabels([f"{(start + p * args.bin_size) / 1000:.1f}" for p in tick_pos])
    ax_hm.set_xlabel(f"position on {chrom} (kb)")
    ax_hm.set_title(f"{name}  {chrom}:{start}-{end}   method={args.method}, k={args.clusters}, "
                    f"bin={args.bin_size}bp({args.bin_mode}), channel={args.cluster_channel}, "
                    f"n_fibers={len(rows_s)}\n"
                    "red = accessible (MSP/m6A, enzyme-tagged), grey = closed",
                    fontsize=10)

    # aggregate tracks
    x = np.arange(n_bins)
    ax_tr.plot(x, open_mat.mean(axis=0), color="#b2182b", lw=1.2, label="open frac")
    ax_tr.plot(x, nuc_mat.mean(axis=0), color="#2166ac", lw=1.2, label="nuc occupancy")
    ax_tr.plot(x, fire_mat.mean(axis=0), color="#1a9850", lw=1.2, label="FIRE frac")
    ax_tr.set_xlim(-0.5, n_bins - 0.5)
    ax_tr.set_ylim(0, 1)
    ax_tr.set_xticks(tick_pos)
    ax_tr.set_xticklabels([f"{(start + p * args.bin_size) / 1000:.1f}" for p in tick_pos])
    ax_tr.set_ylabel("fraction of fibers")
    ax_tr.legend(fontsize=8, loc="upper right", ncol=3)

    # per-cluster per-sample composition bars
    clusters = sorted(set(labels.tolist()))
    y_pos = np.arange(len(clusters))
    height = 0.8 / max(1, len(sample_order))
    for j, sample in enumerate(sample_order):
        n_s = sum(1 for r in rows if r.sample == sample)
        fracs = []
        for c in clusters:
            n_c = sum(1 for r, lab in zip(rows, labels) if lab == c and r.sample == sample)
            fracs.append(n_c / n_s if n_s else 0.0)
        ax_bar.barh(y_pos + (j - (len(sample_order) - 1) / 2) * height, fracs, height=height * 0.92,
                    color=plt.cm.tab10(sample_rank[sample] % 10), label=f"{sample} (n={n_s})")
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels([f"cluster {c}" for c in clusters])
    ax_bar.invert_yaxis()  # cluster 1 on top, matching the heatmap row order
    ax_bar.set_xlabel("fraction within sample")
    ax_bar.set_xlim(0, 1)
    ax_bar.legend(fontsize=7, loc="lower right")
    ax_bar.set_title("cluster proportions per sample", fontsize=9)

    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser.

    Returns:
        Configured argparse.ArgumentParser.
    """
    p = argparse.ArgumentParser(
        description="Cluster single Fiber-seq fibers by accessibility pattern at target regions.")
    p.add_argument("-r", "--regions", type=Path, default=None,
                   help="regions BED (>=4 cols; col4 = region name)")
    p.add_argument("-c", "--locus", action="append", default=[],
                   help="inline locus 'chr:start-end|name' (repeatable)")
    p.add_argument("-m", "--msp", action="append", required=True, metavar="S=P",
                   help="sample=path to msp.bed.gz (accessibility; clustering channel, repeatable)")
    p.add_argument("-u", "--nuc", action="append", default=[], metavar="S=P",
                   help="sample=path to nuc.bed.gz (overlay, repeatable)")
    p.add_argument("-f", "--fire", action="append", default=[], metavar="S=P",
                   help="sample=path to fire_peaks.bed (peak table with FDR col5; "
                        "FIRE overlay = fiber MSP blocks inside peaks passing --fire-fdr, repeatable)")
    p.add_argument("-C", "--cluster-channel", choices=["msp", "fire"], default="msp",
                   help="binarization source for clustering: all MSP blocks (msp, default) "
                        "or only MSP blocks inside significant FIRE peaks (fire)")
    p.add_argument("-q", "--fire-fdr", type=float, default=0.1,
                   help="keep FIRE peaks with FDR <= this (default 0.1)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("-b", "--bin-size", type=int, default=100, help="bin width in bp (default 100)")
    p.add_argument("-B", "--bin-mode", choices=["any", "frac"], default="any",
                   help="bin open when MSP overlap >= --min-msp-bp (any, default) "
                        "or covered fraction >= --bin-threshold (frac)")
    p.add_argument("-P", "--min-msp-bp", type=int, default=10,
                   help="min MSP bp overlap for --bin-mode any (default 10)")
    p.add_argument("-k", "--clusters", type=int, default=3, help="cluster count (default 3)")
    p.add_argument("-e", "--method", choices=["hierarchical", "kmeans"], default="hierarchical",
                   help="clustering method (default hierarchical)")
    p.add_argument("-n", "--min-fibers", type=int, default=10,
                   help="skip regions with fewer spanning fibers (default 10)")
    p.add_argument("-s", "--min-span-frac", type=float, default=0.8,
                   help="fiber must span this fraction of the region (default 0.8)")
    p.add_argument("-t", "--bin-threshold", type=float, default=0.5,
                   help="bin counts as covered when block fraction >= this (default 0.5)")
    p.add_argument("-S", "--sample-order", type=str, default="",
                   help="comma-separated display order; first sample stacks on top")
    p.add_argument("-x", "--seed", type=int, default=0, help="RNG seed for k-means (default 0)")
    p.add_argument("-w", "--overwrite", action="store_true",
                   help="recompute regions whose outputs already exist")
    p.add_argument("-D", "--dry-run", action="store_true", help="log the plan and exit")
    p.add_argument("-l", "--log", type=Path, default=None, help="optional log file path")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the fiber-clustering pipeline over all regions.

    Parameters:
        argv: CLI arguments (None = sys.argv[1:]).

    Returns:
        Process exit code: 0 on full success, 1 when any region failed
        (collected in <outdir>/regions_failed.txt).
    """
    global logger
    args = build_parser().parse_args(argv)
    if args.log is not None:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        logger = setup_logger("FiberCluster", log_file=str(args.log), stream=sys.stdout)

    try:
        msp_map = parse_sample_paths(args.msp, "--msp")
        nuc_map = parse_sample_paths(args.nuc, "--nuc")
        fire_map = parse_sample_paths(args.fire, "--fire")
        regions = ([] if args.regions is None else read_regions_bed(args.regions)) + parse_loci(args.locus)
        if not regions:
            raise ValueError("no regions: pass -r/--regions and/or -c/--locus")
    except ValueError as exc:
        logger.error(str(exc))
        return 2

    if args.sample_order:
        sample_order = [s for s in args.sample_order.split(",") if s]
        unknown = set(sample_order) - set(msp_map)
        if unknown:
            logger.error("--sample-order names unknown samples: %s", sorted(unknown))
            return 2
    else:
        sample_order = list(msp_map)
    for sample in sample_order:
        if sample not in nuc_map:
            logger.warning("sample %s has no --nuc overlay", sample)
        if sample not in fire_map:
            logger.warning("sample %s has no --fire overlay", sample)

    logger.info("regions: %d | samples: %s | bin=%d k=%d method=%s",
                len(regions), sample_order, args.bin_size, args.clusters, args.method)

    if args.dry_run:
        for name, chrom, start, end in regions:
            logger.info("[dry-run] region %s %s:%d-%d", name, chrom, start, end)
            for sample in sample_order:
                logger.info("[dry-run]   msp=%s nuc=%s fire=%s", msp_map[sample],
                            nuc_map.get(sample, "-"), fire_map.get(sample, "-"))
        logger.info("[dry-run] outdir=%s", args.outdir)
        return 0

    # ---- one extraction pass per element file (multi-region filter) -------
    query = merged_intervals(regions)
    logger.info("extracting MSP channels ...")
    msp_all = {s: extract_bed12(p, query) for s, p in msp_map.items()}
    nuc_all = {s: extract_bed12(p, query) for s, p in nuc_map.items()}
    fire_all = {s: extract_fire_peaks(p, query, args.fire_fdr) for s, p in fire_map.items()}

    args.outdir.mkdir(parents=True, exist_ok=True)
    failed: List[Tuple[str, str]] = []
    done = 0
    for region in regions:
        name = region[0]
        region_dir = args.outdir / _safe(name)
        if (region_dir / "cluster_assignments.tsv").is_file() and not args.overwrite:
            logger.info("skip existing region %s (use --overwrite to recompute)", name)
            continue
        try:
            run_region(region,
                       {s: msp_all[s] for s in sample_order},
                       {s: nuc_all.get(s, {}) for s in sample_order},
                       {s: fire_all.get(s, {}) for s in sample_order},
                       sample_order, args, region_dir)
            logger.info("region %s done -> %s", name, region_dir)
            done += 1
        except ValueError as exc:
            logger.warning("region %s failed: %s", name, exc)
            failed.append((name, str(exc)))

    if failed:
        with open(args.outdir / "regions_failed.txt", "w") as fh:
            fh.write("region\treason\n")
            for name, reason in failed:
                fh.write(f"{name}\t{reason}\n")
    logger.info("finished: %d regions done, %d failed", done, len(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
