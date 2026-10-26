#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""IGV-style stacked coverage tracks for selected TE loci.

For each genomic region, renders one figure with one filled coverage
track per sample (stacked, shared x axis), optional TE/gene feature
annotation lane at the bottom with a dedicated label band, and a
group legend strip. A burst locus shows as a coverage peak in the
bursting samples.

Input data
----------
``--track`` : ``NAME=PATH`` (repeatable), one per sample. The reader
is chosen by extension:

* ``.bam``/``.cram`` (pysam): per-base depth over the region from
  ``count_coverage`` (base quality >= 15). Requires a ``.bai``/``.crai``
  index. ``--normalize cpm`` rescales depth by ``1e6 / mapped_reads``.
* ``.bw``/``.bigWig`` (pyBigWig): interval values; NaN -> 0. Assumed
  pre-normalized (e.g. bamCoverage --normalizeUsing CPM).
* ``.bedGraph``/``.bg``/``.bedGraph.gz``/``.bg.gz``: 4-column
  ``chrom start end value`` intervals (0-based half-open). Assumed
  pre-normalized.

``--region`` : ``[NAME=]chr:start-end`` (repeatable; 1-based
inclusive like IGV), or ``--regions-file`` TSV with header columns
``name``, ``chrom``, ``start``, ``end``. One figure per region.

``--annotation`` : optional TSV with header ``sample``/``group``
(optional ``order``) mapping track names to groups; track fill colors
follow the group and a legend strip is added.

``--features-gtf`` / ``--features-bed`` : optional feature models for
the bottom annotation lane.

* GTF: rmsk TE GTF (``*_GENCODE_rmsk_TE.gtf``) or any GTF with
  ``gene_id``/``family_id``/``class_id`` attributes.
* BED: columns ``chrom start end strand gene_id family_id class_id``
  (header optional), e.g. ``*_GENCODE_rmsk_TE.bed``.

Examples
--------
python te_locus_coverage.py \\
    --track mESC1=mESC1.bam --track ciTotiSC1-1=ciTotiSC1-1.bam \\
    --annotation samples.tsv \\
    --region MT2_locus=chr5:12345678-12347000 \\
    --features-gtf GRCm39_GENCODE_rmsk_TE.gtf --bin-size 50 \\
    -o outdir -f png
"""

import argparse
import concurrent.futures
import gzip
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "modules", "function", "bin"),
)
from LogUtil import setup_logger  # noqa: E402

logger = setup_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - environment dependent
    raise SystemExit(
        f"matplotlib is required for {os.path.basename(__file__)}: {exc}. "
        "Install it into the current environment (no auto-install)."
    )

PlotFormat = str  # any matplotlib savefig format: png, pdf, svg, ...

GENE_NAME_RE = re.compile(r'gene_id "([^"]*)"')
FAMILY_ID_RE = re.compile(r'family_id "([^"]*)"')
CLASS_ID_RE = re.compile(r'class_id "([^"]*)"')
GROUP_PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#64B5CD", "#CCB974"]
FEATURE_FILL = "#A0C4E8"
FEATURE_EDGE = "#2B4C7E"
DEFAULT_TRACK_COLOR = "#4C72B0"


def parse_region(region: str) -> Tuple[str, str, int, int]:
    """Parse ``[NAME=]chr:start-end`` into 0-based half-open coordinates.

    Parameters
    ----------
    region : str
        Region string, e.g. ``MT2_locus=chr5:12345678-12347000`` or
        ``chr5:12345678-12347000`` (name defaults to the raw string
        with ``:``/``-`` replaced by ``_``).

    Returns
    -------
    Tuple[str, str, int, int]
        ``(name, chrom, start0, end0)`` with ``start0`` 0-based and
        ``end0`` exclusive.

    Raises
    ------
    ValueError
        If the region string is malformed or end <= start.
    """
    name = region
    if "=" in region:
        name, region = region.split("=", 1)
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
    if match is None:
        raise ValueError(f"Malformed region {region!r}, expected [NAME=]chr:start-end")
    chrom = match.group(1)
    start1, end1 = int(match.group(2)), int(match.group(3))
    if end1 <= start1:
        raise ValueError(f"Region end <= start: {region!r}")
    safe_name = name if name != region else re.sub(r"[:\-]", "_", region)
    return safe_name, chrom, start1 - 1, end1


def bin_means(values: np.ndarray, bin_size: int) -> np.ndarray:
    """Average per-base values into fixed-size bins.

    Parameters
    ----------
    values : np.ndarray
        1-D per-base values.
    bin_size : int
        Bin width in bp. ``<= 1`` returns ``values`` unchanged. The
        tail is zero-padded to a full bin.

    Returns
    -------
    np.ndarray
        Bin-mean values.
    """
    if bin_size <= 1:
        return values
    n_bins = int(np.ceil(values.size / bin_size))
    padded = np.zeros(n_bins * bin_size, dtype=float)
    padded[: values.size] = values
    return padded.reshape(n_bins, bin_size).mean(axis=1)


def bam_coverage(
    bam_path: str,
    chrom: str,
    start0: int,
    end0: int,
    bin_size: int,
    normalize: str = "cpm",
    min_base_quality: int = 15,
) -> np.ndarray:
    """Extract binned coverage from an indexed BAM/CRAM.

    Parameters
    ----------
    bam_path : str
        BAM/CRAM path (``.bai``/``.crai`` index required).
    chrom : str
        Chromosome name.
    start0, end0 : int
        0-based half-open interval.
    bin_size : int
        Bin width in bp (``<= 1`` = per-base).
    normalize : str, optional
        ``"cpm"`` rescales depth by ``1e6 / mapped_reads``; ``"none"``
        keeps raw depth. Default is ``"cpm"``.
    min_base_quality : int, optional
        Minimum base quality counted by ``count_coverage``. Default
        is ``15``.

    Returns
    -------
    np.ndarray
        Binned coverage values.

    Raises
    ------
    ValueError
        If the BAM index is missing or ``normalize`` is invalid.
    """
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit(
            f"pysam is required to read BAM/CRAM {bam_path}: {exc}. "
            "Install it into the current environment (no auto-install)."
        )
    if normalize not in ("cpm", "none"):
        raise ValueError(f"normalize must be cpm/none, got {normalize!r}")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        try:
            total_mapped = sum(stat.mapped for stat in bam.get_index_statistics())
        except ValueError as exc:
            raise ValueError(f"BAM index missing for {bam_path}: {exc}") from exc
        counts = bam.count_coverage(
            chrom, start0, end0, quality_threshold=min_base_quality
        )
    depth = np.sum(np.array(counts, dtype=float), axis=0)
    if normalize == "cpm":
        if total_mapped <= 0:
            raise ValueError(f"No mapped reads in {bam_path}, cannot CPM-normalize")
        depth = depth * 1e6 / total_mapped
    return bin_means(depth, bin_size)


def bigwig_coverage(
    bw_path: str,
    chrom: str,
    start0: int,
    end0: int,
    bin_size: int,
) -> np.ndarray:
    """Extract binned coverage from a bigWig file.

    Parameters
    ----------
    bw_path : str
        bigWig path (``.bw``/``.bigWig``).
    chrom : str
        Chromosome name.
    start0, end0 : int
        0-based half-open interval.
    bin_size : int
        Bin width in bp (``<= 1`` = per-base).

    Returns
    -------
    np.ndarray
        Binned coverage values (NaN bins -> 0).

    Raises
    ------
    ValueError
        If the chromosome is absent from the bigWig.
    """
    try:
        import pyBigWig
    except ImportError as exc:
        raise SystemExit(
            f"pyBigWig is required to read {bw_path}: {exc}. "
            "Install it into the current environment (no auto-install)."
        )
    with pyBigWig.open(bw_path) as bw:
        if chrom not in bw.chroms():
            raise ValueError(f"Chromosome {chrom} not in {bw_path}")
        values = np.asarray(bw.values(chrom, start0, end0, numpy=True), dtype=float)
    values = np.nan_to_num(values, nan=0.0)
    return bin_means(values, bin_size)


def bedgraph_coverage(
    bg_path: str,
    chrom: str,
    start0: int,
    end0: int,
    bin_size: int,
) -> np.ndarray:
    """Extract binned coverage from a bedGraph file.

    Parameters
    ----------
    bg_path : str
        bedGraph path (``.bedGraph``/``.bg``, optionally ``.gz``).
        4 columns: ``chrom start end value`` (0-based half-open).
    chrom : str
        Chromosome name.
    start0, end0 : int
        0-based half-open interval.
    bin_size : int
        Bin width in bp (``<= 1`` = per-base).

    Returns
    -------
    np.ndarray
        Binned coverage values (uncovered positions -> 0).

    Raises
    ------
    ValueError
        If the file has no usable rows.
    """
    opener = gzip.open if bg_path.endswith(".gz") else open
    with opener(bg_path, "rt") as handle:  # type: ignore[operator]
        table = pd.read_csv(
            handle, sep="\t", comment="#", header=None,
            names=["chrom", "start", "end", "value"],
            usecols=[0, 1, 2, 3],
        )
    if table.empty:
        raise ValueError(f"No rows in bedGraph {bg_path}")
    table = table[table["chrom"] == chrom]
    region_len = end0 - start0
    depth = np.zeros(region_len, dtype=float)
    for row in table.itertuples(index=False):
        lo = max(int(row.start), start0) - start0
        hi = min(int(row.end), end0) - start0
        if hi > lo:
            depth[lo:hi] = float(row.value)
    return bin_means(depth, bin_size)


def load_track_coverage(
    path: str,
    chrom: str,
    start0: int,
    end0: int,
    bin_size: int,
    normalize: str = "cpm",
) -> np.ndarray:
    """Dispatch to the reader matching the track file extension.

    Parameters
    ----------
    path : str
        Track file path (BAM/CRAM/bigWig/bedGraph).
    chrom : str
        Chromosome name.
    start0, end0 : int
        0-based half-open interval.
    bin_size : int
        Bin width in bp.
    normalize : str, optional
        BAM normalization (``"cpm"``/``"none"``). bigWig/bedGraph are
        taken as pre-normalized. Default is ``"cpm"``.

    Returns
    -------
    np.ndarray
        Binned coverage values.

    Raises
    ------
    ValueError
        If the extension is unsupported or the file is missing.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Track file not found: {path}")
    lower = path.lower()
    if lower.endswith((".bam", ".cram")):
        values = bam_coverage(path, chrom, start0, end0, bin_size, normalize=normalize)
    elif lower.endswith((".bw", ".bigwig")):
        values = bigwig_coverage(path, chrom, start0, end0, bin_size)
    elif lower.endswith((".bedgraph", ".bg", ".bedgraph.gz", ".bg.gz")):
        values = bedgraph_coverage(path, chrom, start0, end0, bin_size)
    else:
        raise ValueError(
            f"Unsupported track format: {path} "
            "(expected .bam/.cram/.bw/.bigWig/.bedGraph/.bg[.gz])"
        )
    logger.info(f"Track {os.path.basename(path)} {chrom}:{start0 + 1}-{end0}: {values.size} bins")
    return values


def parse_features(path: str) -> pd.DataFrame:
    """Parse a feature annotation file (rmsk TE GTF or 7-column BED).

    Parameters
    ----------
    path : str
        GTF (``gene_id``/``family_id``/``class_id`` attributes) or BED
        (``chrom start end strand gene_id family_id class_id``,
        optional header).

    Returns
    -------
    pd.DataFrame
        Columns: ``chrom``, ``start`` (1-based), ``end``, ``name``,
        ``family``, ``class``.

    Raises
    ------
    ValueError
        If the format is unrecognized or no features parse.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Feature file not found: {path}")

    rows: List[Dict[str, Any]] = []
    if path.lower().endswith((".gtf", ".gff", ".gff2", ".gff3")):
        with open(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                attrs = fields[8]
                name_m = GENE_NAME_RE.search(attrs)
                family_m = FAMILY_ID_RE.search(attrs)
                class_m = CLASS_ID_RE.search(attrs)
                rows.append(
                    {
                        "chrom": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "name": name_m.group(1) if name_m else "NA",
                        "family": family_m.group(1) if family_m else "NA",
                        "class": class_m.group(1) if class_m else "NA",
                    }
                )
    else:
        with open(path) as handle:
            first = handle.readline()
        header: Optional[List[str]] = None
        if first.lower().startswith("chrom"):
            header = first.rstrip("\n").split("\t")
        table = pd.read_csv(path, sep="\t", header=0 if header else None)
        if header is None:
            table.columns = ["chrom", "start", "end", "strand", "gene_id", "family_id", "class_id"][
                : table.shape[1]
            ]
        rename = {"gene_id": "name", "family_id": "family", "class_id": "class"}
        table = table.rename(columns=rename)
        needed = ["chrom", "start", "end"]
        missing = [c for c in needed if c not in table.columns]
        if missing:
            raise ValueError(f"Feature BED {path} missing columns {missing}")
        for col in ("name", "family", "class"):
            if col not in table.columns:
                table[col] = "NA"
        rows = table[["chrom", "start", "end", "name", "family", "class"]].to_dict("records")

    features = pd.DataFrame(rows)
    if features.empty:
        raise ValueError(f"No features parsed from {path}")
    logger.info(f"Parsed {len(features)} features from {path}")
    return features


def load_track_groups(
    annotation_tsv: Optional[str],
    track_names: List[str],
) -> Dict[str, str]:
    """Map track names to group labels.

    Parameters
    ----------
    annotation_tsv : str, optional
        TSV with header ``sample``/``group`` (optional ``order``).
        If None, every track is its own group.
    track_names : List[str]
        Track names in plotting order.

    Returns
    -------
    Dict[str, str]
        Track name -> group label.

    Raises
    ------
    ValueError
        If the annotation lacks required columns or no annotated
        sample matches any track.
    """
    if annotation_tsv is None:
        return {name: name for name in track_names}
    annot = pd.read_csv(annotation_tsv, sep="\t")
    for col in ("sample", "group"):
        if col not in annot.columns:
            raise ValueError(f"annotation column '{col}' missing in {annotation_tsv}")
    annot = annot[annot["sample"].isin(track_names)].copy()
    if annot.empty:
        raise ValueError(f"No annotation rows match track names {track_names}")
    mapping = dict(zip(annot["sample"].astype(str), annot["group"].astype(str)))
    for name in track_names:
        mapping.setdefault(name, "unannotated")
    return mapping


def group_color_map(mapping: Dict[str, str]) -> Dict[str, str]:
    """Assign a stable palette color to each group.

    Parameters
    ----------
    mapping : Dict[str, str]
        Track name -> group label.

    Returns
    -------
    Dict[str, str]
        Group label -> hex color (order of first appearance).
    """
    groups = list(dict.fromkeys(mapping.values()))
    return {group: GROUP_PALETTE[i % len(GROUP_PALETTE)] for i, group in enumerate(groups)}


def plan_feature_layout(
    overlap: pd.DataFrame,
    start0: int,
    end0: int,
    data_per_pt: float,
    label_fontsize: float,
) -> Tuple[List[Tuple[int, float, float]], List[Tuple[int, float, str]]]:
    """Pre-compute feature lane rows and label band rows.

    Parameters
    ----------
    overlap : pd.DataFrame
        Features overlapping the window (``start``/``end`` 1-based).
    start0, end0 : int
        Window coordinates (0-based half-open).
    data_per_pt : float
        Data units per point of text width (estimated from the
        figure width; conservative).
    label_fontsize : float
        Label font size in points.

    Returns
    -------
    Tuple[List[...], List[...]]
        ``(rects, labels)`` where ``rects`` is a list of
        ``(lane_row, x0, x1)`` and ``labels`` is a list of
        ``(label_row, center_x, text)`` (row 0 = topmost label row).
    """
    rects: List[Tuple[int, float, float]] = []
    lane_ends: List[float] = []
    label_candidates: List[Tuple[float, float, str]] = []
    char_pt = label_fontsize * 0.72  # verified text-width factor

    for row in overlap.sort_values(["start", "end"]).itertuples(index=False):
        x0 = max(float(row.start) - 1.0, float(start0))
        x1 = min(float(row.end), float(end0))
        if x1 <= x0:
            continue
        lane = None
        for idx, lane_end in enumerate(lane_ends):
            if x0 > lane_end:
                lane = idx
                lane_ends[idx] = x1
                break
        if lane is None:
            lane_ends.append(x1)
            lane = len(lane_ends) - 1
        rects.append((lane, x0, x1))
        center_x = (x0 + x1) / 2.0
        half_width = len(str(row.name)) * char_pt * data_per_pt / 2.0 * 1.15
        label_candidates.append((center_x, half_width, str(row.name)))

    # greedy row packing of labels (row 0 = topmost)
    rows_end: List[float] = []
    placed: List[Tuple[float, float, str, int]] = []
    for center_x, half_width, text in sorted(label_candidates, key=lambda item: item[0]):
        row = None
        for idx, row_end in enumerate(rows_end):
            if center_x - half_width > row_end:
                row = idx
                rows_end[idx] = center_x + half_width
                break
        if row is None:
            rows_end.append(center_x + half_width)
            row = len(rows_end) - 1
        placed.append((center_x, half_width, text, row))
    n_label_rows = max((row for _, _, _, row in placed), default=-1) + 1
    labels = [(n_label_rows - 1 - row, center_x, text) for center_x, _hw, text, row in placed]
    return rects, labels


def plot_locus_coverage(
    tracks: List[Tuple[str, str]],
    region: Tuple[str, str, int, int],
    out_base: str,
    image_formats: Optional[List[PlotFormat]] = None,
    bin_size: int = 50,
    normalize: str = "cpm",
    share_y: str = "max",
    track_groups: Optional[Dict[str, str]] = None,
    features: Optional[pd.DataFrame] = None,
    flank: int = 0,
    figsize: Tuple[float, float] = (10.0, 2.2),
    dpi: int = 300,
    fill_alpha: float = 0.55,
    line_width: float = 0.8,
    feature_lane_height: float = 0.22,
    label_band_unit: float = 0.3,
    label_fontsize: float = 7.5,
    legend_row_height: float = 0.5,
) -> None:
    """Render the stacked coverage figure for one region.

    Layout (top to bottom): group legend strip (only when track
    groups are given), coverage tracks (one axes each, shared x), the
    feature lane (rectangles only), the label band (feature names,
    row-packed), and the coordinate axis on the bottom-most axes.
    Every band is a dedicated axes, so elements cannot overlap across
    bands.

    Parameters
    ----------
    tracks : List[Tuple[str, str]]
        ``(name, path)`` in plotting order (top to bottom).
    region : Tuple[str, str, int, int]
        ``(name, chrom, start0, end0)`` from :func:`parse_region`.
    out_base : str
        Output path without extension (region name is appended).
    image_formats : List[str], optional
        Save formats. Default is ``["png"]``.
    bin_size : int, optional
        Bin width in bp. Default is ``50``.
    normalize : str, optional
        BAM normalization (``"cpm"``/``"none"``). Default ``"cpm"``.
    share_y : str, optional
        ``"max"`` (default): all tracks share the global maximum;
        ``"track"``: per-track autoscale (per-track max is printed in
        the track label).
    track_groups : Dict[str, str], optional
        Track name -> group label for fill colors + legend strip.
    features : pd.DataFrame, optional
        Feature table from :func:`parse_features`; overlapping
        features are drawn in the bottom annotation lane.
    flank : int, optional
        Padding in bp around the region. Default is ``0``.
    figsize : Tuple[float, float], optional
        Size per track row in inches. Default is ``(10.0, 2.2)``.
    dpi : int, optional
        Output resolution. Default is ``300``.
    fill_alpha : float, optional
        Coverage fill alpha. Default is ``0.55``.
    line_width : float, optional
        Coverage outline width. Default is ``0.8``.
    feature_lane_height : float, optional
        Feature lane height as a fraction of one track row.
        Default is ``0.22``.
    label_band_unit : float, optional
        Label band height per label row as a fraction of one track
        row. Default is ``0.3``.
    label_fontsize : float, optional
        Feature label font size. Default is ``7.5``.
    legend_row_height : float, optional
        Legend strip height as a fraction of one track row. Default
        is ``0.5``.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If ``share_y`` is invalid or no tracks are given.
    """
    if image_formats is None:
        image_formats = ["png"]
    if share_y not in ("max", "track"):
        raise ValueError(f"share_y must be max/track, got {share_y!r}")
    if not tracks:
        raise ValueError("No tracks given")

    name, chrom, start0, end0 = region
    start0 = max(start0 - flank, 0)
    end0 = end0 + flank

    coverages: Dict[str, np.ndarray] = {}
    for track_name, path in tracks:
        coverages[track_name] = load_track_coverage(
            path, chrom, start0, end0, bin_size, normalize=normalize
        )

    n_bins = next(iter(coverages.values())).size
    if bin_size > 1:
        x = start0 + (np.arange(n_bins) + 0.5) * bin_size
    else:
        x = start0 + np.arange(n_bins) + 0.5

    if features is not None:
        overlap = features[
            (features["chrom"] == chrom) & (features["start"] <= end0) & (features["end"] >= start0)
        ].copy()
    else:
        overlap = pd.DataFrame(columns=["chrom", "start", "end", "name", "family", "class"])
    show_features = not overlap.empty

    # ---- layout planning (before any drawing) ----
    axes_width_pt = figsize[0] * 72.0 * 0.77  # conservative default-margin estimate
    data_per_pt = (end0 - start0) / max(axes_width_pt, 1.0)
    if show_features:
        rects, labels = plan_feature_layout(
            overlap, start0, end0, data_per_pt, label_fontsize
        )
        n_lanes = max((lane for lane, _, _ in rects), default=-1) + 1
        n_label_rows = max((row for row, _, _ in labels), default=-1) + 1
    else:
        rects, labels, n_lanes, n_label_rows = [], [], 0, 0

    n_tracks = len(tracks)
    has_legend = track_groups is not None and len(set(track_groups.values())) > 1
    height_ratios: List[float] = [1.0] * n_tracks
    if show_features:
        height_ratios.append(feature_lane_height * max(n_lanes, 1))
        height_ratios.append(label_band_unit * max(n_label_rows, 1) + 0.15)
    legend_height = legend_row_height if has_legend else 0.0
    total_height = figsize[1] * (legend_height + sum(height_ratios))

    gs_rows = n_tracks + (2 if show_features else 0)
    fig = plt.figure(figsize=(figsize[0], total_height))
    gs = fig.add_gridspec(
        gs_rows + (1 if has_legend else 0), 1,
        height_ratios=(([legend_height] + height_ratios) if has_legend else height_ratios),
        hspace=0.15,
        top=0.88 if has_legend else 0.90,
    )

    ax_leg: Any = None
    row_offset = 0
    if has_legend:
        ax_leg = fig.add_subplot(gs[0])
        ax_leg.axis("off")
        row_offset = 1

    group_colors = group_color_map(track_groups) if track_groups else {}
    axes: List[Any] = []
    global_max = max(float(v.max()) for v in coverages.values())
    for idx, (track_name, _path) in enumerate(tracks):
        ax = fig.add_subplot(gs[row_offset + idx])
        axes.append(ax)
        values = coverages[track_name]
        color = (
            group_colors.get(track_groups[track_name], DEFAULT_TRACK_COLOR)
            if track_groups
            else DEFAULT_TRACK_COLOR
        )
        ax.fill_between(x, 0.0, values, step="mid", color=color, alpha=fill_alpha, linewidth=0.0)
        ax.step(x, values, where="mid", color=color, linewidth=line_width)

        if share_y == "track":
            v_max = float(values.max()) if values.size else 0.0
            ax.set_ylim(0.0, v_max * 1.15 if v_max > 0 else 1.0)
            ax.set_ylabel(
                f"{track_name}\n(max {v_max:.1f})", rotation=0, ha="right", va="center", fontsize=9
            )
        else:
            ax.set_ylim(0.0, global_max * 1.15 if global_max > 0 else 1.0)
            ax.set_ylabel(track_name, rotation=0, ha="right", va="center", fontsize=9)
        ax.set_xlim(start0, end0)
        ax.tick_params(axis="x", labelbottom=False, length=2)

    # ---- feature lane: rectangles only, dedicated axes ----
    ax_feat: Any = None
    ax_lab: Any = None
    if show_features:
        ax_feat = fig.add_subplot(gs[row_offset + n_tracks])
        ax_feat.axis("off")
        ax_feat.set_xlim(start0, end0)
        ax_feat.set_ylim(0.0, max(n_lanes, 1) + 0.1)
        for lane, x0, x1 in rects:
            ax_feat.add_patch(
                plt.Rectangle(
                    (x0, lane + 0.25), x1 - x0, 0.5,
                    facecolor=FEATURE_FILL, edgecolor=FEATURE_EDGE, linewidth=0.6, zorder=2,
                )
            )

        ax_lab = fig.add_subplot(gs[row_offset + n_tracks + 1])
        for spine in ax_lab.spines.values():
            spine.set_visible(False)
        ax_lab.set_yticks([])
        ax_lab.set_xlim(start0, end0)
        ax_lab.set_ylim(-0.15, max(n_label_rows, 1))
        for row_idx, center_x, text in labels:
            ax_lab.text(
                center_x, row_idx + 0.5, text,
                ha="center", va="center", fontsize=label_fontsize,
                color="#333333", clip_on=True, zorder=3,
            )

    # ---- coordinate axis on the bottom-most axes ----
    bottom_ax = axes[-1] if ax_lab is None else ax_lab
    tick_positions = np.linspace(start0, end0, 5)
    if end0 - start0 >= 2000:
        tick_labels = [f"{v / 1000.0:.1f}k" for v in tick_positions]
    else:
        tick_labels = [f"{v:.0f}" for v in tick_positions]
    bottom_ax.tick_params(axis="x", labelbottom=True, length=3)
    bottom_ax.set_xticks(tick_positions)
    bottom_ax.set_xticklabels(tick_labels)
    bottom_ax.set_xlabel(f"{chrom} position")

    axes[0].set_title(f"{name}  {chrom}:{start0 + 1}-{end0:,}  (bin={bin_size} bp)", fontsize=11)

    if ax_leg is not None and track_groups is not None:
        handles = [
            plt.Line2D([0], [0], color=group_colors[group], linewidth=6, label=group)
            for group in dict.fromkeys(track_groups.values())
        ]
        ax_leg.legend(handles=handles, loc="center left", ncol=len(handles), frameon=False, fontsize=9)

    for fmt in image_formats:
        out_path = f"{out_base}.{name}.{fmt}"
        fig.savefig(out_path, dpi=dpi)
        logger.info(f"Saved {out_path}")
    plt.close(fig)


def _plot_region_worker(payload: Dict[str, Any]) -> str:
    """Plot one region inside a worker process.

    Parameters
    ----------
    payload : Dict[str, Any]
        Keys: ``tracks``, ``region``, ``out_base``, ``image_formats``,
        ``bin_size``, ``normalize``, ``share_y``, ``track_groups``,
        ``features_path``, ``flank``, ``figsize``, ``dpi``.

    Returns
    -------
    str
        The region name rendered.
    """
    features = parse_features(payload["features_path"]) if payload["features_path"] else None
    plot_locus_coverage(
        tracks=payload["tracks"],
        region=payload["region"],
        out_base=payload["out_base"],
        image_formats=payload["image_formats"],
        bin_size=payload["bin_size"],
        normalize=payload["normalize"],
        share_y=payload["share_y"],
        track_groups=payload["track_groups"],
        features=features,
        flank=payload["flank"],
        figsize=payload["figsize"],
        dpi=payload["dpi"],
    )
    return str(payload["region"][0])


def parse_track_arg(items: List[str]) -> List[Tuple[str, str]]:
    """Parse ``NAME=PATH`` track arguments.

    Parameters
    ----------
    items : List[str]
        Raw ``NAME=PATH`` strings.

    Returns
    -------
    List[Tuple[str, str]]
        ``(name, path)`` pairs in input order.

    Raises
    ------
    ValueError
        If an item has no ``=`` or the path does not exist.
    """
    tracks: List[Tuple[str, str]] = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --track {item!r}, expected NAME=PATH")
        name, path = item.split("=", 1)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Track file not found: {path}")
        tracks.append((name, path))
    return tracks


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser.
    """
    parser = argparse.ArgumentParser(
        description="Stacked IGV-style coverage tracks for TE loci (BAM/bigWig/bedGraph).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--track", action="append", dest="tracks", metavar="NAME=PATH", required=True,
                        help="Coverage track. Repeatable.")
    parser.add_argument("--region", action="append", dest="regions", metavar="[NAME=]CHR:START-END",
                        help="Region (1-based inclusive). Repeatable.")
    parser.add_argument("--regions-file", default=None,
                        help="TSV with header columns name, chrom, start, end (1-based inclusive).")
    parser.add_argument("--annotation", default=None,
                        help="sample/group TSV for track colors + legend strip.")
    parser.add_argument("--features-gtf", default=None, help="Feature GTF for the annotation lane.")
    parser.add_argument("--features-bed", default=None, help="Feature BED for the annotation lane.")
    parser.add_argument("--bin-size", type=int, default=50, help="Bin width in bp (<=1 = per-base).")
    parser.add_argument("--normalize", choices=["cpm", "none"], default="cpm",
                        help="BAM normalization; bigWig/bedGraph are taken as pre-normalized.")
    parser.add_argument("--share-y", choices=["max", "track"], default="max",
                        help="Shared y scale across tracks or per-track autoscale.")
    parser.add_argument("--flank", type=int, default=0, help="Padding around the region in bp.")
    parser.add_argument("--figsize", default="10,2.2", help="Per-track size in inches, 'w,h'.")
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution.")
    parser.add_argument("--fill-alpha", type=float, default=0.55, help="Coverage fill alpha.")
    parser.add_argument("--line-width", type=float, default=0.8, help="Coverage outline width.")
    parser.add_argument("--feature-lane-height", type=float, default=0.22,
                        help="Feature lane height (fraction of one track row).")
    parser.add_argument("--label-band-unit", type=float, default=0.3,
                        help="Label band height per label row (fraction of one track row).")
    parser.add_argument("--label-fontsize", type=float, default=7.5, help="Feature label font size.")
    parser.add_argument("-j", "--threads", type=int, default=1,
                        help="Parallel workers over regions (ProcessPoolExecutor).")
    parser.add_argument("-f", "--format", action="append", dest="formats", metavar="FMT",
                        help="Image format (png/pdf/svg/...). Repeatable. Default: png.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory.")
    parser.add_argument("--out-prefix", default="te_locus_coverage", help="Output file prefix.")
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    """Command-line entry point.

    Parameters
    ----------
    argv : List[str], optional
        Argument list (defaults to ``sys.argv[1:]``).

    Returns
    -------
    None
    """
    args = build_parser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    tracks = parse_track_arg(args.tracks)
    if args.features_gtf and args.features_bed:
        raise SystemExit("Pass only one of --features-gtf / --features-bed")

    regions: List[Tuple[str, str, int, int]] = []
    if args.regions_file:
        table = pd.read_csv(args.regions_file, sep="\t")
        for col in ("name", "chrom", "start", "end"):
            if col not in table.columns:
                raise SystemExit(f"--regions-file missing column {col}")
        for row in table.itertuples(index=False):
            regions.append((str(row.name), str(row.chrom), int(row.start) - 1, int(row.end)))
    for region_str in args.regions or []:
        regions.append(parse_region(region_str))
    if not regions:
        raise SystemExit("No regions: pass --region and/or --regions-file")
    logger.info(f"{len(tracks)} tracks x {len(regions)} regions")

    track_groups = load_track_groups(args.annotation, [name for name, _ in tracks])
    w, h = (float(v) for v in args.figsize.split(","))
    features_path = args.features_gtf or args.features_bed or ""

    payloads = [
        {
            "tracks": tracks,
            "region": region,
            "out_base": os.path.join(args.outdir, args.out_prefix),
            "image_formats": args.formats if args.formats else ["png"],
            "bin_size": args.bin_size,
            "normalize": args.normalize,
            "share_y": args.share_y,
            "track_groups": track_groups,
            "features_path": features_path,
            "flank": args.flank,
            "figsize": (w, h),
            "dpi": args.dpi,
        }
        for region in regions
    ]

    if args.threads > 1 and len(payloads) > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
            for rendered in executor.map(_plot_region_worker, payloads):
                logger.info(f"Region rendered: {rendered}")
    else:
        for payload in payloads:
            logger.info(f"Region rendered: {_plot_region_worker(payload)}")


if __name__ == "__main__":
    main()
