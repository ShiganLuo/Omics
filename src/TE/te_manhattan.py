#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genome-wide Manhattan-style scatter of per-locus TE activity.

Visualizes a "TE burst" as a genome-wide elevation of per-locus TE
values (e.g. DESeq2 log2FoldChange of locus-level TE counts): the x
axis lays out all TE loci by genomic coordinate (chromosomes placed
end-to-end), the y axis shows the chosen statistic.

Input data
----------
``--values`` : TSV, one row per TE locus, with header. Two modes:

* Coordinate mode (default): the table itself carries the locus
  coordinates. Column names are configurable via ``--chrom-col``,
  ``--start-col``, ``--end-col`` (defaults: ``chrom``, ``start``,
  ``end``). Optional ``--name-col``/``--family-col``/``--class-col``
  columns (defaults: ``name``, ``family``, ``class``) label points.
  Typical source: featureCounts run locus-level on
  ``<ref>_GENCODE_rmsk_TE.bed``/``.gtf`` followed by DESeq2.

* Keyed mode (``--te-gtf`` given): the table carries a locus id
  column (``--locus-col``, default ``locus_id``) whose values match
  the ``transcript_id`` attribute of the rmsk TE GTF (e.g.
  ``/data/pub/zhousha/Database/Reference/<species>/GENCODE/<ref>/<ref>_GENCODE_rmsk_TE.gtf``).
  Coordinates and name/family/class are resolved from the GTF.

A numeric value column (``--value-col``, default ``value``) holds the
plotted statistic. An optional significance column (``--signif-col``,
e.g. ``padj``) combines with ``--value-thresh`` to define "burst"
loci.

Notes
-----
* TEcount/TEtranscripts output is FAMILY-level and has no genomic
  coordinates -- it cannot feed this plot. Locus-level counting
  (featureCounts ``-g transcript_id`` / TElocal) is required.
* Chromosome order: ``--chrom-order`` file (one name per line) or
  comma list; default natural order chr1..chrN, chrX, chrY, chrM,
  then remaining names sorted.

Examples
--------
python te_manhattan.py --values locus_deseq2.tsv --te-gtf GRCm39_GENCODE_rmsk_TE.gtf \\
    --locus-col gene_id --value-col log2FoldChange --signif-col padj \\
    --value-thresh 1 --highlight-family MT2_Mm --highlight-family IAP \\
    -o outdir -f png -f pdf
"""

import argparse
import math
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
TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]*)"')
FAMILY_ID_RE = re.compile(r'family_id "([^"]*)"')
CLASS_ID_RE = re.compile(r'class_id "([^"]*)"')

CHROM_COLOR_ODD = "#4C72B0"
CHROM_COLOR_EVEN = "#8DA0CB"
BURST_COLOR = "#C44E52"
HIGHLIGHT_PALETTE = [
    "#DD8452",
    "#55A868",
    "#CCB974",
    "#64B5CD",
    "#8172B3",
    "#937860",
    "#DA8BC3",
    "#8C8C8C",
]


def parse_te_gtf(gtf: str) -> pd.DataFrame:
    """Parse a RepeatMasker-derived TE GTF into a per-locus table.

    Parameters
    ----------
    gtf : str
        Path to the TE GTF (e.g. ``*_GENCODE_rmsk_TE.gtf``). One exon
        line per TE locus with ``gene_id``/``transcript_id``/
        ``family_id``/``class_id`` attributes.

    Returns
    -------
    pd.DataFrame
        Columns: ``locus_id`` (transcript_id), ``chrom``, ``start``,
        ``end``, ``name`` (gene_id), ``family``, ``class``.

    Raises
    ------
    FileNotFoundError
        If ``gtf`` does not exist.
    ValueError
        If no valid feature lines are parsed.
    """
    if not os.path.isfile(gtf):
        raise FileNotFoundError(f"TE GTF not found: {gtf}")

    rows: List[Dict[str, object]] = []
    with open(gtf) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            locus_m = TRANSCRIPT_ID_RE.search(attrs)
            name_m = GENE_NAME_RE.search(attrs)
            family_m = FAMILY_ID_RE.search(attrs)
            class_m = CLASS_ID_RE.search(attrs)
            if locus_m is None:
                logger.warning(f"{gtf}:{line_no}: no transcript_id attribute, skipped")
                continue
            rows.append(
                {
                    "locus_id": locus_m.group(1),
                    "chrom": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "name": name_m.group(1) if name_m else locus_m.group(1),
                    "family": family_m.group(1) if family_m else "NA",
                    "class": class_m.group(1) if class_m else "NA",
                }
            )
    loci = pd.DataFrame(rows)
    if loci.empty:
        raise ValueError(f"No feature lines parsed from TE GTF: {gtf}")
    logger.info(f"Parsed {len(loci)} TE loci from {gtf}")
    return loci


def load_locus_values(
    values_tsv: str,
    value_col: str,
    te_gtf: Optional[str] = None,
    locus_col: str = "locus_id",
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    name_col: str = "name",
    family_col: str = "family",
    class_col: str = "class",
    signif_col: Optional[str] = None,
) -> pd.DataFrame:
    """Load the per-locus value table and standardize its columns.

    Parameters
    ----------
    values_tsv : str
        Path to the value TSV (see module docstring for the contract).
    value_col : str
        Column holding the plotted numeric statistic.
    te_gtf : str, optional
        If given, keyed mode: resolve coordinates via ``locus_col``.
    locus_col, chrom_col, start_col, end_col, name_col, family_col, class_col : str
        Input column names (coordinate mode / key column).
    signif_col : str, optional
        Optional significance column carried through (e.g. ``padj``).

    Returns
    -------
    pd.DataFrame
        Standardized columns: ``chrom``, ``start``, ``end``, ``name``,
        ``family``, ``class``, ``value`` (+ ``signif`` when given).

    Raises
    ------
    ValueError
        If required columns are missing or no rows survive.
    """
    values = pd.read_csv(values_tsv, sep="\t")
    if value_col not in values.columns:
        raise ValueError(f"--value-col '{value_col}' not in {values_tsv}: {list(values.columns)}")

    if te_gtf is not None:
        if locus_col not in values.columns:
            raise ValueError(f"--locus-col '{locus_col}' not in {values_tsv}: {list(values.columns)}")
        loci = parse_te_gtf(te_gtf)
        merged = values.merge(loci, left_on=locus_col, right_on="locus_id", how="inner")
        n_lost = len(values) - len(merged)
        if n_lost > 0:
            logger.warning(f"{n_lost} value rows had no {locus_col} match in {te_gtf}, dropped")
    else:
        needed = [chrom_col, start_col, end_col]
        missing = [c for c in needed if c not in values.columns]
        if missing:
            raise ValueError(
                f"Coordinate columns {missing} not in {values_tsv}; "
                f"pass --te-gtf for keyed mode or fix --chrom-col/--start-col/--end-col"
            )
        merged = values.rename(
            columns={
                chrom_col: "chrom",
                start_col: "start",
                end_col: "end",
            }
        ).copy()
        for src, dst in ((name_col, "name"), (family_col, "family"), (class_col, "class")):
            merged[dst] = merged[src] if src in merged.columns else "NA"

    merged["value"] = pd.to_numeric(merged[value_col], errors="coerce")
    merged = merged[~merged["value"].isna()].copy()
    if signif_col is not None:
        if signif_col not in values.columns:
            raise ValueError(f"--signif-col '{signif_col}' not in {values_tsv}")
        merged["signif"] = pd.to_numeric(merged[signif_col], errors="coerce")
    if merged.empty:
        raise ValueError(f"No usable rows in {values_tsv}")
    logger.info(f"Loaded {len(merged)} loci from {values_tsv}")
    return merged


def resolve_chrom_order(chroms: List[str], chrom_order: Optional[str]) -> List[str]:
    """Return the plotting order of chromosome names.

    Parameters
    ----------
    chroms : List[str]
        Chromosome names present in the data.
    chrom_order : str, optional
        Path to a file with one name per line, or a comma-separated
        list. If None, natural order chr1..chrN, chrX, chrY, chrM,
        then remaining names sorted.

    Returns
    -------
    List[str]
        Chromosome names in plotting order (only names present in
        ``chroms`` are kept).

    Raises
    ------
    FileNotFoundError
        If ``chrom_order`` names a non-existing file.
    """
    present = set(chroms)
    if chrom_order:
        if os.path.isfile(chrom_order):
            with open(chrom_order) as handle:
                wanted = [line.strip() for line in handle if line.strip()]
        else:
            wanted = [name.strip() for name in chrom_order.split(",") if name.strip()]
        ordered = [name for name in wanted if name in present]
        leftover = sorted(present - set(ordered))
        if leftover:
            logger.warning(f"Chromosomes not listed in --chrom-order appended at end: {leftover}")
            ordered.extend(leftover)
        return ordered

    def natural_key(name: str) -> Tuple[int, str]:
        core = name[3:] if name.lower().startswith("chr") else name
        special = {"X": 10 ** 6, "Y": 10 ** 6 + 1, "M": 10 ** 6 + 2, "MT": 10 ** 6 + 2}
        if core.upper() in special:
            return (special[core.upper()], name)
        return (int(core), name) if core.isdigit() else (10 ** 7, name)

    return sorted(present, key=natural_key)


def compute_genome_offsets(
    loci: pd.DataFrame,
    chrom_order: List[str],
    gap_frac: float = 0.02,
) -> pd.DataFrame:
    """Assign cumulative x positions to loci, chromosomes end-to-end.

    Parameters
    ----------
    loci : pd.DataFrame
        Standardized locus table (see :func:`load_locus_values`).
    chrom_order : List[str]
        Chromosome plotting order.
    gap_frac : float, optional
        Gap between chromosomes as a fraction of one chromosome's
        span. Default is ``0.02``.

    Returns
    -------
    pd.DataFrame
        ``loci`` with an added ``xpos`` column and a ``chrom_center``
        column (per-chromosome tick position).

    Raises
    ------
    ValueError
        If ``chrom_order`` does not cover all chromosomes in data.
    """
    unknown = set(loci["chrom"]) - set(chrom_order)
    if unknown:
        raise ValueError(f"chromosomes missing from order: {sorted(unknown)}")

    offsets: Dict[str, float] = {}
    centers: Dict[str, float] = {}
    cursor = 0.0
    for chrom in chrom_order:
        sub = loci.loc[loci["chrom"] == chrom, ["start", "end"]]
        c_min = float(sub[["start", "end"]].min().min())
        c_max = float(sub[["start", "end"]].max().max())
        span = max(c_max - c_min, 1.0)
        offsets[chrom] = cursor - c_min
        centers[chrom] = cursor + span / 2.0
        cursor += span * (1.0 + gap_frac)

    out = loci.copy()
    out["xpos"] = out["start"].astype(float) + out["chrom"].map(offsets)
    out["chrom_center"] = out["chrom"].map(centers)
    return out


def select_highlights(
    loci: pd.DataFrame,
    highlight_families: Optional[List[str]],
    burst_thresh: Optional[float],
    signif_thresh: Optional[float],
    burst_sign: str = "up",
) -> pd.DataFrame:
    """Flag highlight-family and burst loci.

    Parameters
    ----------
    loci : pd.DataFrame
        Standardized locus table.
    highlight_families : List[str], optional
        Substring patterns matched case-insensitively against
        ``name``/``family``/``class``.
    burst_thresh : float, optional
        Value threshold for burst loci (>= if ``burst_sign='up'``,
        <= if ``burst_sign='down'``). None disables the burst filter.
    signif_thresh : float, optional
        If given (and a ``signif`` column exists), burst loci must
        additionally satisfy ``signif <= signif_thresh``.
    burst_sign : str, optional
        ``"up"``, ``"down"`` or ``"both"``. Default is ``"up"``.

    Returns
    -------
    pd.DataFrame
        ``loci`` with boolean ``is_highlight`` and ``is_burst`` columns.

    Raises
    ------
    ValueError
        If ``burst_sign`` is not one of ``up``/``down``/``both``.
    """
    if burst_sign not in ("up", "down", "both"):
        raise ValueError(f"burst_sign must be up/down/both, got {burst_sign!r}")

    out = loci.copy()
    haystack = (
        out["name"].astype(str).str.lower()
        + "|"
        + out["family"].astype(str).str.lower()
        + "|"
        + out["class"].astype(str).str.lower()
    )
    if highlight_families:
        mask = pd.Series(False, index=out.index)
        for pattern in highlight_families:
            mask |= haystack.str.contains(pattern.lower(), regex=False)
        out["is_highlight"] = mask
    else:
        out["is_highlight"] = False

    if burst_thresh is None:
        out["is_burst"] = False
    else:
        if burst_sign == "up":
            burst = out["value"] >= burst_thresh
        elif burst_sign == "down":
            burst = out["value"] <= burst_thresh
        else:
            burst = out["value"].abs() >= abs(burst_thresh)
        if signif_thresh is not None and "signif" in out.columns:
            burst = burst & (out["signif"] <= signif_thresh)
        out["is_burst"] = burst.fillna(False)
    logger.info(
        f"Selected {int(out['is_highlight'].sum())} highlight loci, "
        f"{int(out['is_burst'].sum())} burst loci"
    )
    return out


def dodge_labels(
    ax: "plt.Axes",
    items: List[Tuple[float, float, str]],
    fontsize: float = 8.0,
    max_iter: int = 60,
    push_pt: float = 2.0,
) -> None:
    """Draw horizontal point labels with iterative collision push.

    Labels start just above their point; colliding labels are pushed
    vertically (alternating up/down preference by y rank) until no
    pair overlaps or ``max_iter`` is reached.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes (data coordinates).
    items : List[Tuple[float, float, str]]
        ``(x, y, text)`` entries in data coordinates.
    fontsize : float, optional
        Label font size in points. Default is ``8.0``.
    max_iter : int, optional
        Maximum push iterations. Default is ``60``.
    push_pt : float, optional
        Vertical push per iteration in points. Default is ``2.0``.

    Returns
    -------
    None
    """
    if not items:
        return
    fig = ax.figure
    fig.canvas.draw()  # resolve annotation positions before measuring extents
    renderer = fig.canvas.get_renderer()
    texts = []
    for x, y, label in items:
        text = ax.annotate(
            label,
            xy=(x, y),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            color="#333333",
            clip_on=True,
            zorder=5,
        )
        texts.append(text)

    # deterministic stacking: cluster labels by x overlap, then stack each
    # cluster bottom-up (no push iterations -> no oscillation), flipping the
    # cluster below its anchor when stacking would overflow the axes top.
    base_boxes = [t.get_window_extent(renderer=renderer) for t in texts]
    px_per_pt = fig.dpi / 72.0
    gap_px = push_pt * px_per_pt
    n = len(texts)
    offsets = [0.0] * n

    order = sorted(range(n), key=lambda i: base_boxes[i].x0)
    x_margin = 10.0  # px: labels within this horizontal gap stack vertically too
    clusters: List[List[int]] = []
    current: List[int] = []
    current_x1 = -float("inf")
    for i in order:
        if current and base_boxes[i].x0 >= current_x1 + x_margin:
            clusters.append(current)
            current = []
            current_x1 = -float("inf")
        current.append(i)
        current_x1 = max(current_x1, base_boxes[i].x1)
    if current:
        clusters.append(current)

    axes_top = ax.get_window_extent(renderer=renderer).y1
    anchor_offset_pt = 6.0
    for cluster in clusters:
        cluster.sort(key=lambda i: base_boxes[i].y0)  # bottom-most first
        stacked_max = max(base_boxes[i].y1 + (anchor_offset_pt * px_per_pt) for i in cluster)
        flip = stacked_max > axes_top
        if flip:
            # stack downward from a below-point anchor instead
            drop_pt = anchor_offset_pt + fontsize * 1.4
            for i in cluster:
                texts[i].set_position((texts[i].get_position()[0], -drop_pt))
            fig.canvas.draw()
            base_boxes = [t.get_window_extent(renderer=renderer) for t in texts]
            cluster.sort(key=lambda i: base_boxes[i].y1, reverse=True)  # top-most first
            prev_bottom = None
            for i in cluster:
                b = base_boxes[i]
                if prev_bottom is not None and b.y1 > prev_bottom - gap_px:
                    offsets[i] = (prev_bottom - gap_px) - b.y1
                prev_bottom = b.y0 + offsets[i]
        else:
            prev_top = None
            for i in cluster:
                b = base_boxes[i]
                if prev_top is not None and b.y0 < prev_top + gap_px:
                    offsets[i] = (prev_top + gap_px) - b.y0
                prev_top = b.y1 + offsets[i]

    for idx, text in enumerate(texts):
        base_x, base_y = text.get_position()
        text.set_position((base_x, base_y + offsets[idx] / px_per_pt))
    logger.debug(f"Label dodge finished for {len(items)} labels in {len(clusters)} clusters")


def thin_tick_labels(
    positions: List[float],
    labels: List[str],
    span: float,
    axes_width_pt: float,
    rotation_deg: float = 45.0,
    fontsize: float = 8.0,
    gap_pt: float = 4.0,
) -> Tuple[List[float], List[str]]:
    """Drop tick labels whose horizontal footprint collides with a kept one.

    Chromosomes with tiny observed spans (unplaced scaffolds) have
    nearly identical centers, so their 45-degree labels pile up. The
    labels are thinned greedily left-to-right; all tick positions are
    preserved by the caller as minor ticks.

    Parameters
    ----------
    positions : List[float]
        Tick positions in data coordinates.
    labels : List[str]
        Tick label strings (same length as ``positions``).
    span : float
        Data range covered by the axes (``xlim`` width).
    axes_width_pt : float
        Axes width in points (for the data-per-point conversion).
    rotation_deg : float, optional
        Label rotation in degrees. Default is ``45.0``.
    fontsize : float, optional
        Label font size in points. Default is ``8.0``.
    gap_pt : float, optional
        Extra clear gap between labels in points. Default ``4.0``.

    Returns
    -------
    Tuple[List[float], List[str]]
        Kept positions and labels.

    Notes
    -----
    For rotated collinear labels the binding constraint is the
    perpendicular distance between neighboring label baselines:
    required center distance = ``(height + gap) / sin(rotation)``.
    The axis-aligned footprint (``width*cos + height*sin``)
    overestimates and over-thins. Near-horizontal labels
    (``rotation < 15``) pack along x by width instead.
    """
    theta = math.radians(rotation_deg)
    data_per_pt = span / max(axes_width_pt, 1.0)
    height_pt = fontsize * 1.2
    if abs(rotation_deg) < 15.0:
        required_pt = max((len(label) * fontsize * 0.72 for label in labels), default=fontsize) + gap_pt
    else:
        required_pt = (height_pt + gap_pt) / max(abs(math.sin(theta)), 0.25)
    required = required_pt * data_per_pt
    kept_pos: List[float] = []
    kept_lab: List[str] = []
    last_pos = -float("inf")
    for pos, label in sorted(zip(positions, labels)):
        if pos - last_pos > required:
            kept_pos.append(pos)
            kept_lab.append(label)
            last_pos = pos
    return kept_pos, kept_lab


def plot_manhattan(
    loci: pd.DataFrame,
    out_base: str,
    image_formats: Optional[List[PlotFormat]] = None,
    figsize: Tuple[float, float] = (14.0, 5.5),
    dpi: int = 300,
    point_size: float = 3.0,
    point_alpha: float = 0.35,
    highlight_size: float = 9.0,
    max_points: int = 300_000,
    seed: int = 0,
    ylabel: str = "value",
    title: Optional[str] = None,
    gap_frac: float = 0.02,
    label_top: int = 0,
    hlines: Optional[List[float]] = None,
    chrom_order: Optional[str] = None,
) -> None:
    """Render the genome-wide TE locus Manhattan plot.

    Parameters
    ----------
    loci : pd.DataFrame
        Output of :func:`select_highlights` (needs ``xpos``,
        ``chrom_center`` via :func:`compute_genome_offsets`).
    out_base : str
        Output path without extension.
    image_formats : List[str], optional
        Save formats. Default is ``["png"]``.
    figsize : Tuple[float, float], optional
        Figure size in inches. Default is ``(14.0, 5.5)``.
    dpi : int, optional
        Output resolution. Default is ``300``.
    point_size : float, optional
        Background point size. Default is ``3.0``.
    point_alpha : float, optional
        Background point alpha. Default is ``0.35``.
    highlight_size : float, optional
        Burst/highlight point size. Default is ``9.0``.
    max_points : int, optional
        Rendering cap for background points (burst/highlight points
        are always drawn in full). ``<= 0`` disables the cap.
        Default is ``300000`` (3.8M rmsk loci would otherwise take
        >15 min to raster).
    seed : int, optional
        RNG seed for background subsampling. Default is ``0``.
    ylabel : str, optional
        Y axis label. Default is ``"value"``.
    title : str, optional
        Figure title. Default is ``None`` (no title).
    gap_frac : float, optional
        Unused here (offsets precomputed); kept for CLI symmetry.
    label_top : int, optional
        Label the N loci with the largest ``|value|`` among the
        marked loci (burst + highlight). Default is ``0`` (no
        labels); keep N small (<= 30) -- one text artist per label.
    hlines : List[float], optional
        Horizontal reference lines (e.g. ``[1, -1]``).
    chrom_order : str, optional
        Unused here (order precomputed); kept for CLI symmetry.

    Returns
    -------
    None
    """
    if image_formats is None:
        image_formats = ["png"]

    n_highlight_families = loci.loc[loci["is_highlight"], "family"].nunique() if "is_highlight" in loci else 0
    has_legend = bool(loci["is_burst"].any()) or n_highlight_families > 0
    fig = plt.figure(figsize=figsize)
    ax_leg: Any = None
    if has_legend:
        gs = fig.add_gridspec(2, 1, height_ratios=[8.0, 1.0], hspace=0.25)
        ax = fig.add_subplot(gs[0])
        ax_leg = fig.add_subplot(gs[1])
        ax_leg.axis("off")
    else:
        ax = fig.add_subplot()

    chrom_names = resolve_chrom_order(list(loci["chrom"].unique()), chrom_order)
    background = loci[~loci["is_highlight"] & ~loci["is_burst"]]
    if max_points > 0 and len(background) > max_points:
        logger.info(
            f"Background loci downsampled for rendering: {len(background)} -> {max_points} "
            "(burst/highlight points kept in full)"
        )
        background = background.sample(n=max_points, random_state=seed)
    for idx, chrom in enumerate(chrom_names):
        sub = background[background["chrom"] == chrom]
        color = CHROM_COLOR_ODD if idx % 2 == 0 else CHROM_COLOR_EVEN
        ax.scatter(
            sub["xpos"], sub["value"], s=point_size, c=color, alpha=point_alpha,
            linewidths=0, rasterized=True, zorder=2,
        )

    legend_handles: List[object] = []
    highlight_family_names = sorted(loci.loc[loci["is_highlight"], "family"].astype(str).unique())
    for idx, fam_name in enumerate(highlight_family_names):
        sub = loci[loci["is_highlight"] & (loci["family"].astype(str) == fam_name)]
        handle = ax.scatter(
            sub["xpos"], sub["value"], s=highlight_size,
            c=HIGHLIGHT_PALETTE[idx % len(HIGHLIGHT_PALETTE)], alpha=0.9,
            linewidths=0, zorder=3, label=fam_name,
        )
        legend_handles.append(handle)

    if loci["is_burst"].any():
        burst = loci[loci["is_burst"]]
        handle = ax.scatter(
            burst["xpos"], burst["value"], s=highlight_size + 1.0, c=BURST_COLOR,
            alpha=0.95, linewidths=0, zorder=4, label="burst locus",
        )
        legend_handles.append(handle)

    centers = loci.groupby("chrom")["chrom_center"].first()
    ax.set_ylabel(ylabel)
    x_max = float(loci["xpos"].max()) * 1.01
    ax.set_xlim(0, x_max)
    all_pos = [float(centers[c]) for c in chrom_names if c in centers]
    all_lab = [c for c in chrom_names if c in centers]
    axes_width_pt = figsize[0] * 72.0 * 0.65  # conservative (y labels eat left margin)
    kept_pos, kept_lab = thin_tick_labels(all_pos, all_lab, x_max, axes_width_pt)
    if len(kept_lab) < len(all_lab):
        logger.info(
            f"Chromosome tick labels thinned to fit: kept {len(kept_lab)}/{len(all_lab)} "
            "(all chromosomes keep a minor tick; use --chrom-order/--values to restrict scope)"
        )
    ax.set_xticks(all_pos, minor=True)
    ax.tick_params(axis="x", which="minor", length=2)
    ax.set_xticks(kept_pos)
    ax.set_xticklabels(kept_lab, fontsize=8, rotation=45, ha="right")
    if title:
        ax.set_title(title, fontsize=12)

    if hlines:
        for level in hlines:
            ax.axhline(level, color="#555555", linestyle="--", linewidth=0.8, zorder=1)
            ax.annotate(
                f"{level:g}",
                xy=(1.0, level),
                xycoords=("axes fraction", "data"),
                ha="right",
                va="bottom" if level >= 0 else "top",
                fontsize=8,
                color="#555555",
                annotation_clip=False,
                zorder=6,
            )

    label_items: List[Tuple[float, float, str]] = []
    if label_top > 0:
        marked = loci[loci["is_burst"] | loci["is_highlight"]]
        top = marked.reindex(
            marked["value"].abs().sort_values(ascending=False).index
        ).head(label_top)
        for _, row in top.iterrows():
            label_items.append((float(row["xpos"]), float(row["value"]), str(row["name"])))
    dodge_labels(ax, label_items)

    if ax_leg is not None:
        ax_leg.legend(
            handles=legend_handles,
            loc="center left",
            ncol=max(1, len(legend_handles)),
            frameon=False,
            fontsize=9,
        )

    fig.align_ylabels()
    for fmt in image_formats:
        out_path = f"{out_base}.{fmt}"
        fig.savefig(out_path, dpi=dpi)
        logger.info(f"Saved {out_path}")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser.
    """
    parser = argparse.ArgumentParser(
        description="Genome-wide Manhattan-style scatter of per-locus TE activity.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--values", required=True, help="Per-locus value TSV (see module docstring).")
    parser.add_argument("--te-gtf", default=None,
                        help="rmsk TE GTF for keyed mode (coordinates resolved via --locus-col).")
    parser.add_argument("--locus-col", default="locus_id", help="Locus id column (keyed mode).")
    parser.add_argument("--chrom-col", default="chrom", help="Chromosome column (coordinate mode).")
    parser.add_argument("--start-col", default="start", help="Start column (coordinate mode, 1-based GTF-style).")
    parser.add_argument("--end-col", default="end", help="End column (coordinate mode).")
    parser.add_argument("--name-col", default="name", help="Locus name column (optional).")
    parser.add_argument("--family-col", default="family", help="TE family column (optional).")
    parser.add_argument("--class-col", default="class", help="TE class column (optional).")
    parser.add_argument("--value-col", default="value", help="Numeric column to plot on y.")
    parser.add_argument("--signif-col", default=None, help="Optional significance column (e.g. padj).")
    parser.add_argument("--signif-thresh", type=float, default=None,
                        help="Significance cutoff applied to --signif-col for burst calls.")
    parser.add_argument("--value-thresh", type=float, default=None,
                        help="Value cutoff for burst loci (with --burst-sign).")
    parser.add_argument("--burst-sign", choices=["up", "down", "both"], default="up",
                        help="Direction of the burst value cutoff.")
    parser.add_argument("--highlight-family", action="append", dest="highlight_families",
                        metavar="PATTERN", help="Substring matched on name/family/class. Repeatable.")
    parser.add_argument("--label-top", type=int, default=0,
                        help="Label the N marked loci (burst+highlight) with largest |value|. Keep <= 30.")
    parser.add_argument("--hline", action="append", type=float, dest="hlines",
                        metavar="Y", help="Horizontal reference line. Repeatable.")
    parser.add_argument("--chrom-order", default=None,
                        help="Chromosome order: file (one per line) or comma list.")
    parser.add_argument("--ylabel", default="value", help="Y axis label.")
    parser.add_argument("--title", default=None, help="Figure title (default: none).")
    parser.add_argument("--gap-frac", type=float, default=0.02,
                        help="Gap between chromosomes as a fraction of chromosome span.")
    parser.add_argument("--figsize", default="14,5.5", help="Figure size in inches, 'w,h'.")
    parser.add_argument("--point-size", type=float, default=3.0, help="Background point size.")
    parser.add_argument("--point-alpha", type=float, default=0.35, help="Background point alpha.")
    parser.add_argument("--highlight-size", type=float, default=9.0, help="Burst/highlight point size.")
    parser.add_argument("--max-points", type=int, default=300000,
                        help="Rendering cap for background points (burst/highlight kept in full). <=0 disables.")
    parser.add_argument("--seed", type=int, default=0, help="RNG seed for background subsampling.")
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution.")
    parser.add_argument("-f", "--format", action="append", dest="formats", metavar="FMT",
                        help="Image format (png/pdf/svg/...). Repeatable. Default: png.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory.")
    parser.add_argument("--out-prefix", default="te_manhattan", help="Output file prefix.")
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

    loci = load_locus_values(
        values_tsv=args.values,
        value_col=args.value_col,
        te_gtf=args.te_gtf,
        locus_col=args.locus_col,
        chrom_col=args.chrom_col,
        start_col=args.start_col,
        end_col=args.end_col,
        name_col=args.name_col,
        family_col=args.family_col,
        class_col=args.class_col,
        signif_col=args.signif_col,
    )
    loci = select_highlights(
        loci,
        highlight_families=args.highlight_families,
        burst_thresh=args.value_thresh,
        signif_thresh=args.signif_thresh,
        burst_sign=args.burst_sign,
    )
    chrom_order = resolve_chrom_order(list(loci["chrom"].unique()), args.chrom_order)
    loci = compute_genome_offsets(loci, chrom_order, gap_frac=args.gap_frac)

    width, height = (float(v) for v in args.figsize.split(","))
    plot_manhattan(
        loci=loci,
        out_base=os.path.join(args.outdir, args.out_prefix),
        image_formats=args.formats if args.formats else ["png"],
        figsize=(width, height),
        dpi=args.dpi,
        point_size=args.point_size,
        point_alpha=args.point_alpha,
        highlight_size=args.highlight_size,
        max_points=args.max_points,
        seed=args.seed,
        ylabel=args.ylabel,
        title=args.title,
        gap_frac=args.gap_frac,
        label_top=args.label_top,
        hlines=args.hlines,
        chrom_order=args.chrom_order,
    )


if __name__ == "__main__":
    main()
