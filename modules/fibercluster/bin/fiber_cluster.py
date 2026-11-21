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
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple

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
FireSample = Dict[str, Dict[str, List[Tuple[int, int, float]]]]  # sample -> chrom -> peaks (start, end, fdr)
FireSeg = Tuple[int, int, int]  # (start, end, tier) tier: 2=q<=0.01, 1=0.01<q<=0.05, 0=0.05<q<=0.1


class FiberRow(NamedTuple):
    """One clustering unit: a fiber passing the span filter at one region."""
    read: str
    sample: str
    span_frac: float
    open: np.ndarray
    f_start: int
    f_end: int
    msp_blocks: Blocks
    nuc_blocks: Blocks
    fire_segs: List[FireSeg]


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


def parse_bed12_line(line: str, query: Dict[str, List[Tuple[int, int]]],
                     chrom: str) -> Optional[FiberRec]:
    """Parse one BED12 line into a fiber record, dropping sentinel blocks.

    fibertools-rs v0.13.1 writes two sentinel blocks (first size 0, last
    size 1) at the read boundaries; only interior blocks are real calls.

    Parameters:
        line: raw BED12 text line.
        query: per-chrom merged query intervals (for block clipping tests).
        chrom: chromosome of the record (already matched against query).

    Returns:
        (read_name, blocks, fiber_start, fiber_end) or None when malformed.
    """
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        return None
    start, end = int(f[1]), int(f[2])
    spans = query[chrom]
    if not _overlaps_any(spans, start, end):
        return None
    n_blocks = int(f[9])
    sizes = [int(x) for x in f[10].rstrip(",").split(",")]
    offsets = [int(x) for x in f[11].rstrip(",").split(",")]
    if len(sizes) != n_blocks or len(offsets) != n_blocks:
        logger.warning("malformed BED12 block lists at %s:%d-%d; record skipped",
                       chrom, start, end)
        return None
    blocks: Blocks = []
    for size, off in zip(sizes[1:-1], offsets[1:-1]):
        if size <= 0:
            continue
        b_start, b_end = start + off, start + off + size
        if _overlaps_any(spans, b_start, b_end):
            blocks.append((b_start, b_end))
    return f[3], blocks, start, end


def _query_key(query: Dict[str, List[Tuple[int, int]]]) -> str:
    """Stable md5 of the merged query intervals (cache validity key)."""
    import hashlib
    h = hashlib.md5()
    for chrom in sorted(query):
        for s, e in query[chrom]:
            h.update(f"{chrom}:{s}-{e}\n".encode())
    return h.hexdigest()


def extract_bed12(path: Path, query: Dict[str, List[Tuple[int, int]]],
                  cache_dir: Optional[Path] = None) -> Dict[str, List[FiberRec]]:
    """Collect region-overlapping fibers from one gzipped BED12 element file.

    Access order (fastest first):
    1. tabix index (beside the file or inside cache_dir) -> per-region jump
       reads, seconds regardless of file size;
    2. region cache (<cache_dir>/<name>.subset.bed.gz, keyed by the query
       md5) -> decompress only the cached subset;
    3. full-file stream (gzip) -> and, when cache_dir is given, write the
       subset cache atomically so later runs take path 2.

    Parameters:
        path: bgzip/gzip BED12 element file (m6a/msp/nuc/cpg).
        query: per-chromosome merged query intervals.
        cache_dir: optional cache directory for the subset and index files.

    Returns:
        chrom -> list of (read_name, blocks, fiber_start, fiber_end).
    """
    import time

    out: Dict[str, List[FiberRec]] = defaultdict(list)
    t0 = time.time()
    tbi_local = Path(str(path) + ".tbi")
    tbi_cached = cache_dir / f"{path.name}.tbi" if cache_dir is not None else None

    def parse_into(fh_lines) -> int:
        n = 0
        for line in fh_lines:
            f = line.split("\t", 1)
            if not f or len(f) < 2:
                continue
            chrom = f[0]
            if chrom not in query:
                continue
            rec = parse_bed12_line(line, query, chrom)
            if rec is not None:
                out[chrom].append(rec)
                n += 1
        return n

    if tbi_local.is_file() or (tbi_cached is not None and tbi_cached.is_file()):
        import pysam  # bundled in fiberseq.sif  # pyright: ignore[reportMissingImports]

        idx = str(tbi_cached) if (tbi_cached is not None and tbi_cached.is_file()) else None
        with pysam.TabixFile(str(path), index=idx) as tbx:
            n = 0
            for chrom, spans in query.items():
                if chrom not in tbx.contigs:
                    continue
                for s, e in spans:
                    n += parse_into(tbx.fetch(chrom, s, e))
        logger.info("extract %s: %d records via tabix in %.1fs", path.name, n, time.time() - t0)
        return out

    if cache_dir is not None:
        cache_file = cache_dir / f"{path.name}.subset.bed.gz"
        key_file = cache_dir / f"{path.name}.subset.key"
        key = _query_key(query)
        if cache_file.is_file() and key_file.is_file() and key_file.read_text().strip() == key:
            with gzip.open(cache_file, "rt") as fh:  # type: ignore[operator]
                n = parse_into(fh)
            logger.info("extract %s: %d records via cache in %.1fs", path.name, n, time.time() - t0)
            return out

    # full-file stream (and refresh the cache along the way)
    kept_lines: List[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        for line in fh:
            f = line.split("\t", 1)
            if not f or len(f) < 2:
                continue
            chrom = f[0]
            if chrom not in query:
                continue
            rec = parse_bed12_line(line, query, chrom)
            if rec is not None:
                out[chrom].append(rec)
                kept_lines.append(line if line.endswith("\n") else line + "\n")
    logger.info("extract %s: %d records via full stream in %.1fs",
                path.name, sum(len(v) for v in out.values()), time.time() - t0)

    if cache_dir is not None and kept_lines:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = cache_dir / f"{path.name}.subset.bed.gz"
        key_file = cache_dir / f"{path.name}.subset.key"
        tmp = cache_file.with_suffix(".tmp")
        with gzip.open(tmp, "wt") as fh:  # type: ignore[operator]
            fh.writelines(kept_lines)
        os.replace(tmp, cache_file)
        key_file.write_text(_query_key(query))
        logger.info("cache written: %s (%d records)", cache_file, len(kept_lines))
    return out


def extract_fire_peaks(path: Path, query: Dict[str, List[Tuple[int, int]]],
                       fdr_max: float) -> Dict[str, List[Tuple[int, int, float]]]:
    """Load significant FIRE peaks overlapping the query intervals.

    fire_peaks.bed columns: #chrom peak_start peak_end peak_max FDR coverage
    fire_coverage score ... (fibertools `ft fire --extract` peak table).
    Peaks with FDR <= fdr_max are kept per chromosome, sorted by start, with
    their FDR retained (needed for the paper-style q-tier colors). The
    per-fiber fire.bed.gz is NOT used here: it carries all MSP intervals, so
    significance can only be filtered at the peak level.

    Parameters:
        path: fire_peaks.bed (plain or gzipped).
        query: per-chromosome merged query intervals.
        fdr_max: keep peaks with FDR <= this value.

    Returns:
        chrom -> sorted (start, end, fdr) peak records.
    """
    raw: Dict[str, List[Tuple[int, int, float]]] = defaultdict(list)
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
                raw[chrom].append((start, end, fdr))
    return {chrom: sorted(pts) for chrom, pts in raw.items()}


def fire_tier(fdr: float) -> int:
    """Map a peak FDR to the paper's three-tier coding.

    Parameters:
        fdr: peak-level false discovery rate.

    Returns:
        2 for q<=0.01 (RE/P, dark brown), 1 for 0.01<q<=0.05 (bright red),
        0 for 0.05<q<=0.10 (orange).
    """
    if fdr <= 0.01:
        return 2
    if fdr <= 0.05:
        return 1
    return 0


def intersect_peaks(blocks: Blocks, peaks: Sequence[Tuple[int, int, float]]) -> List[FireSeg]:
    """Intersect a fiber's MSP blocks with significant peaks, keeping tiers.

    Parameters:
        blocks: sorted MSP blocks of one fiber.
        peaks: sorted (start, end, fdr) peaks.

    Returns:
        (start, end, tier) segments where this fiber's MSPs fall inside a
        significant FIRE peak; tier follows fire_tier(fdr).
    """
    out: List[FireSeg] = []
    i = j = 0
    while i < len(blocks) and j < len(peaks):
        lo = max(blocks[i][0], peaks[j][0])
        hi = min(blocks[i][1], peaks[j][1])
        if hi > lo:
            out.append((lo, hi, fire_tier(peaks[j][2])))
        if blocks[i][1] <= peaks[j][1]:
            i += 1
        else:
            j += 1
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


def subtract_blocks(blocks: Blocks, cuts: Sequence[Tuple[int, ...]]) -> Blocks:
    """Intervals of `blocks` not covered by `cuts` (cuts may carry extra fields).

    Parameters:
        blocks: sorted [start, end) intervals.
        cuts: sorted cut intervals; only fields [0], [1] are read.

    Returns:
        The parts of `blocks` outside every cut, sorted.
    """
    out: Blocks = []
    for b0, b1 in blocks:
        cur = b0
        for cut in cuts:
            c0, c1 = cut[0], cut[1]
            if c1 <= cur or c0 >= b1:
                continue
            if c0 > cur:
                out.append((cur, min(c0, b1)))
            cur = max(cur, c1)
            if cur >= b1:
                break
        if cur < b1:
            out.append((cur, b1))
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


def cluster_matrix(mat: np.ndarray, k: int, method: str, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    """Cluster binary fiber x bin patterns.

    Parameters:
        mat: (n_fibers, n_bins) binary accessibility matrix.
        k: requested cluster count.
        method: 'hierarchical' (average-linkage, Hamming) or 'kmeans'.
        seed: RNG seed for k-means init.

    Returns:
        (labels, order): 1-based cluster labels per fiber relabeled so
        cluster 1 has the highest mean openness, and a display order that
        keeps clusters contiguous with tree (dendrogram leaf) ordering
        within each cluster for hierarchical, or centroid-distance ordering
        for k-means.

    Errors:
        ValueError when there are fewer fibers than clusters.
    """
    from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
    from scipy.spatial.distance import pdist

    n = mat.shape[0]
    if n < k:
        raise ValueError(f"only {n} fibers for k={k}")
    if method == "hierarchical":
        dist = pdist(mat.astype(float), metric="hamming")
        z = linkage(dist, method="average")
        labels = fcluster(z, t=k, criterion="maxclust")
        order = leaves_list(z)
    elif method == "kmeans":
        from scipy.cluster.vq import kmeans2

        np.random.seed(seed)
        centroids, labels = kmeans2(mat.astype(float), k, minit="++", iter=50)
        labels = labels + 1
        d = np.array([np.linalg.norm(mat[i].astype(float) - centroids[labels[i] - 1])
                      for i in range(n)])
        order = np.argsort(d, kind="stable")
    else:
        raise ValueError(f"unknown method {method!r}")

    # Relabel by mean openness (descending) so labels are stable across runs.
    sorted_c = sorted(set(int(x) for x in labels),
                      key=lambda c: -float(mat[labels == c].mean()))
    remap = {old: new for new, old in enumerate(sorted_c, start=1)}
    out = np.array([remap[int(x)] for x in labels], dtype=int)
    if len(set(out.tolist())) < k:
        logger.warning("clustering produced %d non-empty clusters, expected %d "
                       "(identical/tied patterns in this channel)", len(set(out.tolist())), k)
    # Display order: clusters contiguous (in label order), tree/centroid order within.
    rank = {int(idx): r for r, idx in enumerate(order)}
    display = sorted(range(n), key=lambda i: (int(out[i]), rank[i]))
    return out, np.array(display, dtype=int)


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
    nuc_lookup: Dict[str, Dict[str, Blocks]] = {
        s: {r_name: blocks for r_name, blocks, _fs, _fe in recs.get(chrom, [])}
        for s, recs in nuc.items()}
    rows: List[FiberRow] = []
    for sample in sample_order:
        peaks_b = fire.get(sample, {}).get(chrom, [])
        for read_name, blocks, f_start, f_end in msp.get(sample, {}).get(chrom, []):
            overlap = min(f_end, end) - max(f_start, start)
            if overlap <= 0:
                continue
            span_frac = overlap / reg_len
            if span_frac < args.min_span_frac:
                continue
            fire_segs = intersect_peaks(sorted(blocks), peaks_b) if peaks_b else []
            channel_blocks: Blocks = ([(s, e) for s, e, _t in fire_segs]
                                      if args.cluster_channel == "fire" else blocks)
            cov_bp = bin_covered_bp(channel_blocks, start, end, bin_size)
            rows.append(FiberRow(
                read_name, sample, span_frac,
                binarize(cov_bp, start, end, bin_size, args.bin_mode,
                         args.bin_threshold, args.min_msp_bp),
                f_start, f_end, blocks,
                nuc_lookup.get(sample, {}).get(read_name, []),
                fire_segs))
    n_fibers = len(rows)
    if n_fibers < args.min_fibers:
        raise ValueError(f"too_few_fibers: {n_fibers} < --min-fibers {args.min_fibers}")

    open_mat = np.vstack([r.open for r in rows]).astype(int)
    labels, disp_order = cluster_matrix(open_mat, args.clusters, args.method, args.seed)

    # ---- overlay matrices (nuc occupancy, FIRE presence) -----------------
    def blocks_matrix(get_blocks: Callable[[FiberRow], Blocks], mode: str) -> np.ndarray:
        mat = np.zeros((n_fibers, n_bins), dtype=float)
        for i, r in enumerate(rows):
            cov_bp = bin_covered_bp(sorted(get_blocks(r)), start, end, bin_size)
            mat[i] = binarize(cov_bp, start, end, bin_size, mode,
                              args.bin_threshold, args.min_msp_bp)
        return mat

    nuc_mat = blocks_matrix(lambda r: r.nuc_blocks, "frac") if nuc else np.zeros((n_fibers, n_bins))
    fire_mat = (blocks_matrix(lambda r: [(s, e) for s, e, _t in r.fire_segs], args.bin_mode)
                if fire else np.zeros((n_fibers, n_bins)))

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
    plot_region(region, rows, labels, disp_order, sample_order, args,
                outdir / f"{_safe(name)}.clusters.png")


def _safe(name: str) -> str:
    """Sanitize a region name for use in filenames."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


# ---------------------------------------------------------------------------
# fiberbrowser palette (copied verbatim from modules/fiberbrowser/bin/
# fiberbrowser.py; colors sampled from the paper figure at 500 dpi pixels)
# ---------------------------------------------------------------------------
NUC_LINE = "#1A1A1A"      # black molecule line = nucleosome / neutral DNA
LINKER = "#8A78B6"        # purple rectangle = linker DNA (non-FIRE MSP)
FIRE_ORANGE = "#F0A246"   # MSP inside a FIRE peak, 0.05 < q <= 0.10
FIRE_RED = "#E01B24"      # MSP inside a FIRE peak, 0.01 < q <= 0.05
FIRE_DKRED = "#7D3219"    # MSP inside a FIRE peak, q <= 0.01 (RE/P)
LABEL_BOX = "#D8ECD8"     # pale green sample label box
LABEL_EC = "#7A9A7A"
COV_NUC = "#B0B1B5"       # grey   = nucleosome coverage (stack top)
COV_LNK = "#8A78B6"       # purple = linker (MSP) coverage (stack middle)
COV_FIRE = "#C0392B"      # red    = FIRE coverage (stack bottom)
TIER_COLORS = {0: FIRE_ORANGE, 1: FIRE_RED, 2: FIRE_DKRED}


def plot_region(region: Region, rows: List[FiberRow], labels: np.ndarray,
                disp_order: np.ndarray, sample_order: List[str], args: argparse.Namespace,
                out_png: Path) -> None:
    """Render the per-region figure as fiberbrowser-style single-fiber tracks.

    Left panel (paper style, per modules/fiberbrowser): one row per fiber
    drawn as a black molecule line (nucleosome / neutral DNA) with purple
    linker (MSP) rectangles and FIRE-tier-colored MSP rectangles on top
    (orange 0.05<q<=0.10, red 0.01<q<=0.05, dark brown q<=0.01). Rows are
    grouped into cluster bands and within each band into per-sample
    sub-stacks, annotated fiberbrowser-style: bracket (tinted per cluster) +
    pale-green rounded label box with the sample id and its fiber count in
    this cluster. The cluster band id and total n sit at the far left. All
    text lives outside the frame; no title (region identity is in the
    output filename).

    Right sub-plots (one per cluster, optional via --right-tracks): per-base
    coverage stacks in the fiberbrowser coverage-track style — red FIRE
    coverage at the bottom, purple linker (non-FIRE MSP) coverage in the
    middle, grey nucleosome coverage on top — drawn INSIDE the cluster
    band's own vertical span (same row coordinate system, so each track is
    y-aligned with exactly the rows it summarizes). Coverage = number of
    class fibers carrying the element at that base; the stack is scaled to
    the band height (max total coverage = band top). Each layer can be
    enabled individually or disabled with 'none'.

    Row order: (cluster, sample, clustering tree leaf order). Fibers beyond
    --max-fibers are subsampled proportionally per sub-stack (display only;
    every printed count uses the full set).

    Parameters:
        region: (name, chrom, start, end).
        rows: per-fiber clustering rows with element blocks.
        labels: cluster labels per fiber.
        disp_order: tree-ordered fiber indices from cluster_matrix.
        sample_order: display order of samples (reference first).
        args: CLI namespace (max_fibers, right_tracks, ...).
        out_png: output image path.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    _name, chrom, start, end = region
    sample_rank = {s: i for i, s in enumerate(sample_order)}

    # ---- display rows: (cluster, sample, tree order), optional subsample --
    rank = {int(idx): r for r, idx in enumerate(disp_order)}
    order_all = sorted(range(len(rows)),
                       key=lambda i: (int(labels[i]), sample_rank.get(rows[i].sample, 99),
                                      rank[i]))
    subsampled = len(order_all) > args.max_fibers
    if subsampled:
        groups: Dict[Tuple[int, str], List[int]] = {}
        for i in order_all:
            groups.setdefault((int(labels[i]), rows[i].sample), []).append(i)
        keep = set()
        for members in groups.values():
            take = max(1, round(len(members) * args.max_fibers / len(order_all)))
            step = len(members) / take
            keep.update(members[int(k * step)] for k in range(take))
        order = [i for i in order_all if i in keep]
    else:
        order = order_all
    n_show = len(order)

    # sub-stack runs: contiguous (cluster, sample) blocks in display order
    runs: List[List] = []  # each: [cluster, sample, y_lo, y_hi]
    for j, i in enumerate(order):
        key = (int(labels[i]), rows[i].sample)
        if runs and (runs[-1][0], runs[-1][1]) == key:
            runs[-1][3] = j
        else:
            runs.append([key[0], key[1], j, j])

    fig = plt.figure(figsize=(13.5, 8.2))
    # right sub-plot width >= half of the main panel (visual balance)
    ax = fig.add_axes((0.26, 0.10, 0.46, 0.84))
    axr = fig.add_axes((0.74, 0.10, 0.25, 0.84), sharey=ax)
    ax.set_xlim(start, end)
    ax.set_ylim(n_show - 0.2, -0.6)  # fiberbrowser row convention
    ax.set_yticks([])
    # no frame: keep only the coordinate axes (bottom x axis)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    axr.set_xlim(start, end)
    axr.set_yticks([])
    axr.set_xticks([])
    for spine in ("top", "right", "left", "bottom"):
        axr.spines[spine].set_visible(False)

    def cluster_color(c: int) -> Tuple[float, float, float, float]:
        """Cluster id -> tab10 color (same mapping as cluster_proportions.png)."""
        return plt.cm.tab10((c - 1) % 10)

    def y2f(y: float) -> float:
        """Row coordinate -> axes fraction measured from the BOTTOM.

        The data axis is row-inverted (y=0 at top), but transAxes fractions
        always grow upward — map explicitly or every label/bracket mirrors.
        """
        return (n_show - 0.2 - y) / (n_show + 0.4)

    # ---- fiber rows: black molecule line, purple linkers, tiered FIRE MSPs
    for j, i in enumerate(order):
        r = rows[i]
        lo_f, hi_f = max(r.f_start, start), min(r.f_end, end)
        if hi_f > lo_f:
            ax.plot([lo_f, hi_f], [j, j], lw=0.9, color=NUC_LINE, zorder=1,
                    solid_capstyle="butt")
        for s0, e0 in r.msp_blocks:
            lo_b, hi_b = max(s0, start), min(e0, end)
            if hi_b > lo_b:
                ax.add_patch(Rectangle((lo_b, j - 0.25), hi_b - lo_b, 0.50,
                                       facecolor=LINKER, ec="none", zorder=2))
        for s0, e0, tier in r.fire_segs:
            lo_b, hi_b = max(s0, start), min(e0, end)
            if hi_b > lo_b:
                ax.add_patch(Rectangle((lo_b, j - 0.28), hi_b - lo_b, 0.56,
                                       facecolor=TIER_COLORS[tier], ec="none", zorder=3))

    # ---- per-run bracket (tinted per cluster) + green label box ----------
    for c, sample, y_lo, y_hi in runs:
        f_lo, f_hi = y2f(y_lo - 0.5), y2f(y_hi + 0.5)
        cc = cluster_color(c)
        ax.plot([-0.012, -0.012], [f_lo, f_hi], transform=ax.transAxes,
                color=cc, lw=1.2, clip_on=False, zorder=5)
        ax.plot([-0.012, -0.004], [f_lo, f_lo], transform=ax.transAxes,
                color=cc, lw=1.2, clip_on=False, zorder=5)
        ax.plot([-0.012, -0.004], [f_hi, f_hi], transform=ax.transAxes,
                color=cc, lw=1.2, clip_on=False, zorder=5)
        n_full = sum(1 for r, lab in zip(rows, labels)
                     if int(lab) == c and r.sample == sample)
        n_here = y_hi - y_lo + 1
        box = f"{sample}\n(n={n_full}" + (f", show {n_here})" if subsampled else ")")
        ax.text(-0.022, (f_lo + f_hi) / 2, box, transform=ax.transAxes,
                fontsize=6.5, ha="right", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.35", facecolor=LABEL_BOX,
                          edgecolor=LABEL_EC, lw=0.6))

    # ---- cluster band labels at the far left (outside, dodged) -----------
    band_center: Dict[int, List[float]] = {}
    for c, _s, y_lo, y_hi in runs:
        band_center.setdefault(c, []).append((y_lo + y_hi) / 2)
    centers = [float(np.mean(v)) for v in band_center.values()]
    lab_frac: List[float] = []
    prev_f = None
    for cy in centers:  # clusters top -> bottom; dodge pushes labels downward
        cf = y2f(cy)
        y_lab = cf if prev_f is None else min(cf, prev_f - 0.075)
        lab_frac.append(y_lab)
        prev_f = y_lab
    if lab_frac and lab_frac[-1] < 0.03:
        shift_l = 0.03 - lab_frac[-1]
        lab_frac = [y + shift_l for y in lab_frac]
    for (c, _centers), y_lab in zip(band_center.items(), lab_frac):
        n_c = int((labels == c).sum())
        ax.text(-0.22, y_lab, f"cluster {c} (n={n_c})", transform=ax.transAxes,
                fontsize=8, ha="right", va="center", zorder=5, color=cluster_color(c))

    # ---- right sub-plots: per-cluster coverage tracks (fiberbrowser style) --
    wanted = [t.strip() for t in args.right_tracks.split(",")
              if t.strip() and t.strip() != "none"]
    x_pos = np.arange(start, end)
    cov_colors = {"fire": COV_FIRE, "lnk": COV_LNK, "nuc": COV_NUC}
    spans: Dict[int, List[int]] = {}
    for c, _s, y_lo, y_hi in runs:
        spans.setdefault(c, [y_lo, y_hi])
        spans[c][0] = min(spans[c][0], y_lo)
        spans[c][1] = max(spans[c][1], y_hi)
    for c in sorted(spans):
        y_lo, y_hi = spans[c]
        y_t, y_b = y_lo - 0.4, y_hi + 0.4  # band top/bottom in row coordinates
        cov = {k: np.zeros(end - start, dtype=float) for k in ("fire", "lnk", "nuc")}

        def paint(arr: np.ndarray, blocks: Blocks) -> None:
            """Add per-base coverage of one fiber's element blocks."""
            for s0, e0 in blocks:
                lo_b, hi_b = max(s0, start), min(e0, end)
                if hi_b > lo_b:
                    arr[lo_b - start:hi_b - start] += 1

        for r, lab in zip(rows, labels):
            if int(lab) != c:
                continue
            paint(cov["fire"], [(s, e) for s, e, _t in r.fire_segs])
            paint(cov["lnk"], subtract_blocks(r.msp_blocks, r.fire_segs))
            paint(cov["nuc"], r.nuc_blocks)
        # fiberbrowser coverage stack (red FIRE bottom, purple linker middle,
        # grey nucleosome top): draw the widest cumulative band first so the
        # inner layers stay visible
        order = [k for k in ("fire", "lnk", "nuc") if k in wanted]
        cum = np.zeros(end - start, dtype=float)
        layers: List[Tuple[np.ndarray, str]] = []
        for k in order:
            cum = cum + cov[k]
            layers.append((cum.copy(), cov_colors[k]))
        scale = max(float(cum.max()), 1.0) if layers else 1.0
        for cum_arr, color in reversed(layers):
            yv = y_b - (cum_arr / scale) * (y_b - y_t)  # max coverage = band top
            axr.fill_between(x_pos, y_b, yv, color=color, lw=0)

    # ---- bottom coordinate ruler + element legend ------------------------
    ticks = np.linspace(start, end, 6)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t / 1000:.1f}" for t in ticks], fontsize=7)
    ax.set_xlabel(f"position on {chrom} (kb)", fontsize=8)
    # ---- bottom legends: one per sub-plot, each under its own panel ------
    handles_left = [
        Line2D([0], [0], color=NUC_LINE, lw=1.2, label="nucleosome / neutral DNA"),
        Rectangle((0, 0), 1, 1, facecolor=LINKER, edgecolor="none", label="linker (MSP)"),
        Rectangle((0, 0), 1, 1, facecolor=FIRE_ORANGE, edgecolor="none",
                  label="FIRE 0.05<q<=0.10"),
        Rectangle((0, 0), 1, 1, facecolor=FIRE_RED, edgecolor="none",
                  label="FIRE 0.01<q<=0.05"),
        Rectangle((0, 0), 1, 1, facecolor=FIRE_DKRED, edgecolor="none",
                  label="FIRE q<=0.01 (RE/P)")]
    handles_cov = []
    for key, lab in (("fire", "coverage: FIRE (in significant RE)"),
                     ("lnk", "coverage: linker (MSP)"),
                     ("nuc", "coverage: nucleosome")):
        if key in wanted:
            handles_cov.append(Rectangle((0, 0), 1, 1, facecolor=cov_colors[key],
                                         edgecolor="none", label=lab))
    fig.legend(handles=handles_left, loc="lower center", bbox_to_anchor=(0.49, -0.01),
               ncol=5, fontsize=7, frameon=False)
    if handles_cov:
        fig.legend(handles=handles_cov, loc="lower center", bbox_to_anchor=(0.865, -0.02),
                   ncol=1, fontsize=7, frameon=False)

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
    p.add_argument("-a", "--max-fibers", type=int, default=300,
                   help="display subsample cap per figure (statistics always use all "
                        "fibers; default 300)")
    p.add_argument("-T", "--right-tracks", type=str, default="fire,lnk,nuc",
                   help="comma list of layers in the per-cluster right coverage "
                        "track (fiberbrowser stack, bottom->top fire/lnk/nuc): "
                        "fire (FIRE coverage), lnk (linker MSP coverage), nuc "
                        "(nucleosome coverage); 'none' disables (default fire,lnk,nuc)")
    p.add_argument("-y", "--cache-dir", type=Path, default=None,
                   help="extraction cache dir: region subsets are cached here "
                        "(keyed by the region set) and tabix indexes looked up "
                        "here or beside the source files; makes reruns seconds-fast")
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
    msp_all = {s: extract_bed12(p, query, args.cache_dir) for s, p in msp_map.items()}
    nuc_all = {s: extract_bed12(p, query, args.cache_dir) for s, p in nuc_map.items()}
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
