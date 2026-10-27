#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genome-wide Manhattan-style scatter of per-locus TE activity.

Visualizes a "TE burst" as a genome-wide elevation of per-locus TE
values (e.g. DESeq2 log2FoldChange of locus-level TE counts): the x
axis lays out all TE loci by genomic coordinate (chromosomes placed
end-to-end), the y axis shows the chosen statistic.

Input data
----------
Two input modes:

* ``--deseq`` mode (recommended): a TElocal DESeq2 result TSV
  (``*_TE.tsv``, locus-id index + statistics + per-sample count
  columns), a ``group.tsv`` (``sample``/``group``), group names
  (``--control``/``--experiment``) and the TElocal locus index
  (``--loci-tsv <ref>_GENCODE_rmsk_TE.gtf.locInd``). Per-locus TPM,
  delta TPM (y axis) and mean TPM (point size) are computed
  internally, all in the same TPM scale.

* ``--values`` mode: a TSV, one row per TE locus, carrying locus
  coordinates (``--chrom-col``/``--start-col``/``--end-col``,
  defaults ``chrom``/``start``/``end``), optional
  ``--name-col``/``--family-col``/``--class-col`` label columns, a
  numeric ``--value-col`` for the plotted statistic, an optional
  ``--signif-col`` (e.g. ``padj``) and ``--size-col`` for point
  sizing.

Notes
-----
* TEcount/TEtranscripts output is FAMILY-level and has no genomic
  coordinates -- it cannot feed this plot. Locus-level counting
  (TElocal) is required.
* Chromosome order: ``--chrom-order`` file (one name per line) or
  comma list; default natural order chr1..chrN, chrX, chrY, chrM,
  then remaining names sorted.

Examples
--------
python te_manhattan.py --deseq mESC_vs_TLSC.TElocal_TE.tsv \\
    --group group.tsv --control mESC --experiment TLSC \\
    --loci-tsv GRCm39_GENCODE_rmsk_TE.gtf.locInd \\
    --signif-col padj --signif-thresh 0.05 --value-thresh 0 \\
    --highlight-family MERVL -o outdir -f png
"""

import argparse
import math
import os
import pickle
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


CHROM_COLOR_ODD = "#4C72B0"
CHROM_COLOR_EVEN = "#8DA0CB"
# Burst/down loci are colored by direction: red = up, blue = down.
BURST_COLOR_UP = "#C44E52"
BURST_COLOR_DOWN = "#2166AC"
DOWN_POINT_COLOR = "#1A4E8A"  # darker blue for background down points
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


def load_locind(locind_path: str) -> pd.DataFrame:
    """Dump per-locus coordinates from a TElocal ``.locInd`` pickle.

    The file is a pickled ``TElocal_Toolkit.TEindex.TEfeatures`` index:
    per chromosome an AVL tree of bin nodes, where each node's
    ``__namelist`` maps start -> [(name_idx, end)]. Names resolve via
    ``getEleName`` to the strand-free ``name:family:class`` locus id
    matching the DESeq2 result table's keys.

    Parameters
    ----------
    locind_path : str
        Path to the ``<ref>_GENCODE_rmsk_TE.gtf.locInd`` file.

    Returns
    -------
    pd.DataFrame
        Columns ``locus_id``, ``chrom``, ``start``, ``end``.
    """
    with open(locind_path, "rb") as handle:
        index = pickle.load(handle)
    rows: List[Tuple[str, str, int, int]] = []
    get_ele_name = index.getEleName
    for chrom, tree in index.indexlist.items():
        stack = [tree.root]
        while stack:
            node = stack.pop()
            if node is None:
                continue
            namelist = getattr(node, "_Node__namelist", {})
            for start, entries in namelist.items():
                for name_idx, end in entries:
                    rows.append((get_ele_name(name_idx), chrom, start, end))
            stack.extend([node.left, node.right])
    loci = pd.DataFrame(rows, columns=["locus_id", "chrom", "start", "end"])
    logger.info(f"Parsed {len(loci)} TE loci from {locind_path}")
    # cache the parsed table next to the pickle: unpickling the TElocal
    # index takes ~20-30 s per call, the parquet read ~0.3 s
    cache_path = locind_path + ".parquet"
    try:
        loci.to_parquet(cache_path)
        logger.info(f"Cached locus table to {cache_path}")
    except Exception as exc:  # pragma: no cover - cache is best effort
        logger.warning(f"Cannot write locus cache {cache_path}: {exc}")
    return loci


def load_loci_table(locind_path: str) -> pd.DataFrame:
    """Load the locus coordinate table, via parquet cache when possible.

    Parameters
    ----------
    locind_path : str
        Path to a TElocal ``.locInd`` pickle, its ``.parquet`` cache, or
        a plain TSV.

    Returns
    -------
    pd.DataFrame
        Columns ``locus_id``, ``chrom``, ``start``, ``end``.
    """
    if locind_path.endswith(".parquet"):
        return pd.read_parquet(locind_path)
    if locind_path.endswith(".locInd"):
        cache_path = locind_path + ".parquet"
        if os.path.isfile(cache_path):
            try:
                loci = pd.read_parquet(cache_path)
                logger.info(
                    f"Loaded {len(loci)} cached loci from {cache_path} "
                    f"(delete it to re-parse {locind_path})"
                )
                return loci
            except Exception as exc:
                logger.warning(f"Cache {cache_path} unreadable ({exc}); re-parsing")
        return load_locind(locind_path)
    return pd.read_csv(locind_path, sep="\t")


def load_locus_values(
    values_tsv: str,
    value_col: str,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    name_col: str = "name",
    family_col: str = "family",
    class_col: str = "class",
    signif_col: Optional[str] = None,
    size_col: Optional[str] = None,
) -> pd.DataFrame:
    """Load the per-locus value table and standardize its columns.

    Parameters
    ----------
    values_tsv : str
        Path to the value TSV (see module docstring for the contract).
    value_col : str
        Column holding the plotted numeric statistic.
    chrom_col, start_col, end_col, name_col, family_col, class_col : str
        Input column names (coordinate mode).
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

    needed = [chrom_col, start_col, end_col]
    missing = [c for c in needed if c not in values.columns]
    if missing:
        raise ValueError(
            f"Coordinate columns {missing} not in {values_tsv}; "
            f"fix --chrom-col/--start-col/--end-col"
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
    if size_col is not None:
        if size_col not in values.columns:
            raise ValueError(f"--size-col '{size_col}' not in {values_tsv}")
        merged["size_value"] = pd.to_numeric(merged[size_col], errors="coerce").fillna(0.0)
    if merged.empty:
        raise ValueError(f"No usable rows in {values_tsv}")
    logger.info(f"Loaded {len(merged)} loci from {values_tsv}")
    return merged


def load_deseq_tpm(
    deseq_tsv: str,
    group_tsv: str,
    control: str,
    experiment: str,
    loci_tsv: Optional[str] = None,
    signif_col: str = "padj",
) -> pd.DataFrame:
    """Build the plotting table from a DESeq2 result TSV plus group file.

    The DESeq2 result table (e.g. TElocal ``*_TE.tsv``) carries per-sample
    normalized count columns named as in ``group_tsv``. Per-locus TPM is
    computed inside this function: counts are divided by locus lengths
    (from the coordinates), scaled to 1e6 per sample, then averaged per
    group; ``delta_tpm`` = experiment mean - control mean and
    ``mean_tpm`` = grand mean, both in the same TPM scale.

    Parameters
    ----------
    deseq_tsv : str
        DESeq2 result TSV; first column is the locus id index
        (e.g. ``AmnL2-1_dup1083:AmnL2-1:L2:LINE``), plus statistic
        columns and per-sample count columns.
    group_tsv : str
        Two-column TSV (``sample``/``group``) mapping the count columns
        to groups.
    control, experiment : str
        Group names in ``group_tsv`` defining the two conditions.
    loci_tsv : str, optional
        Locus coordinate source: a TElocal ``.locInd`` pickle (parsed
        via :func:`load_locind`) or a plain TSV with
        ``locus_id``/``chrom``/``start``/``end`` columns. Required.
    signif_col : str, optional
        Significance column carried through. Default ``padj``.

    Returns
    -------
    pd.DataFrame
        Standardized columns: ``chrom``, ``start``, ``end``, ``name``,
        ``family``, ``class``, ``locus_id``, ``value`` (= delta TPM),
        ``mean_tpm``, ``signif``.

    Raises
    ------
    ValueError
        If inputs are missing columns or no rows survive the joins.
    """
    deseq = pd.read_csv(deseq_tsv, sep="\t", index_col=0)
    deseq = deseq.reset_index()
    deseq = deseq.rename(columns={deseq.columns[0]: "locus_id"})
    deseq["locus_id"] = deseq["locus_id"].astype(str)
    parts = deseq["locus_id"].str.split(":", expand=True)
    if parts.shape[1] >= 3:
        deseq["name"] = parts[0].str.rsplit("_dup", n=1).str[0]
        deseq["family"] = parts[1]
        deseq["class"] = parts[2]
    else:
        deseq["name"], deseq["family"], deseq["class"] = "NA", "NA", "NA"

    groups = pd.read_csv(group_tsv, sep="\t")
    if not {"sample", "group"}.issubset(groups.columns):
        raise ValueError(f"{group_tsv} must have 'sample' and 'group' columns")
    ctrl_samples = groups.loc[groups["group"] == control, "sample"].astype(str).tolist()
    exp_samples = groups.loc[groups["group"] == experiment, "sample"].astype(str).tolist()
    missing = [s for s in ctrl_samples + exp_samples if s not in deseq.columns]
    if missing:
        raise ValueError(f"samples {missing} from {group_tsv} not in {deseq_tsv}")
    if not ctrl_samples or not exp_samples:
        raise ValueError(f"groups '{control}'/'{experiment}' not found in {group_tsv}")

    # coordinates
    if loci_tsv is None:
        raise ValueError("need --loci-tsv (a TElocal .locInd file or a TSV) for locus coordinates")
    if loci_tsv.endswith(".locInd") or loci_tsv.endswith(".parquet"):
        loci = load_loci_table(loci_tsv)
    else:
        loci = pd.read_csv(loci_tsv, sep="\t")
    agg = (
        loci.sort_values(["locus_id", "start", "end"])
        .groupby("locus_id", as_index=False)
        .agg(chrom=("chrom", "first"), start=("start", "min"), end=("end", "max"))
    )
    merged = deseq.merge(agg, on="locus_id", how="inner")
    merged["length"] = (merged["end"] - merged["start"]).clip(lower=1)

    counts = merged[ctrl_samples + exp_samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    rpk = counts.div(merged["length"].values, axis=0)
    tpm = rpk.div(rpk.sum(axis=0).values, axis=1) * 1e6
    merged["ctrl_tpm"] = tpm[ctrl_samples].mean(axis=1).values
    merged["exp_tpm"] = tpm[exp_samples].mean(axis=1).values
    merged["value"] = merged["exp_tpm"] - merged["ctrl_tpm"]
    merged["mean_tpm"] = (merged["exp_tpm"] + merged["ctrl_tpm"]) / 2.0
    if signif_col in merged.columns:
        merged["signif"] = pd.to_numeric(merged[signif_col], errors="coerce")
    else:
        logger.warning(f"'{signif_col}' not in {deseq_tsv}; significance disabled")
    keep = ["chrom", "start", "end", "name", "family", "class", "locus_id",
            "length", "ctrl_tpm", "exp_tpm", "value", "mean_tpm", "signif"]
    merged = merged[[c for c in keep if c in merged.columns]]
    if merged.empty:
        raise ValueError(f"No usable rows joined from {deseq_tsv}")
    logger.info(
        f"Loaded {len(merged)} loci from {deseq_tsv} "
        f"(ctrl={ctrl_samples}, exp={exp_samples})"
    )
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


def transform_value(
    loci: pd.DataFrame,
    y_transform: str = "log2",
) -> pd.DataFrame:
    """Transform the ``value`` column for plotting.

    Parameters
    ----------
    loci : pd.DataFrame
        Standardized locus table with a numeric ``value`` column.
    y_transform : str, optional
        ``"log2"`` (default), ``"log10"`` or ``"none"`` (identity).
        Signed values use the symmetric log transform
        ``sign(x) * log(1 + |x|)`` so negative deltas stay plottable.

    Returns
    -------
    pd.DataFrame
        ``loci`` with ``value`` replaced by the transformed values.

    Raises
    ------
    ValueError
        If ``y_transform`` is not one of ``log2``/``log10``/``none``.
    """
    if y_transform == "none":
        return loci
    if y_transform not in ("log2", "log10"):
        raise ValueError(f"y_transform must be log2/log10/none, got {y_transform!r}")
    base = 2.0 if y_transform == "log2" else 10.0
    out = loci.copy()
    out["value"] = np.sign(out["value"]) * np.log1p(out["value"].abs()) / np.log(base)
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


def size_scale(
    size_value: "pd.Series",
    s_min: float,
    s_max: float,
) -> "pd.Series":
    """Map an abundance column (e.g. baseMean) to marker areas.

    Uses ``log10(1 + x)`` scaling so the heavy right tail of count
    data does not crowd everything at the maximum.

    Parameters
    ----------
    size_value : pd.Series
        Non-negative numeric values (NaN filled with 0 upstream).
    s_min, s_max : float
        Output marker area range (matplotlib ``s`` units).

    Returns
    -------
    pd.Series
        Marker areas in ``[s_min, s_max]``.
    """
    v = size_value.clip(lower=0.0)
    lv = np.log10(1.0 + v)
    lo, hi = float(lv.min()), float(lv.max())
    if hi <= lo:
        return pd.Series(s_min, index=size_value.index)
    return s_min + (lv - lo) / (hi - lo) * (s_max - s_min)


def detect_y_break(values: "np.ndarray", min_gap_frac: float = 0.10) -> Optional[Tuple[float, float]]:
    """Detect a y-axis break band worth collapsing (auto broken axis).

    A break is justified only when a contiguous empty band splits the
    data into a dense body and a sparse tail, and collapsing it saves a
    meaningful fraction of the axis height. The scan runs on sorted
    values: the best candidate gap must (a) isolate <= 10% of points as
    the outer tail beyond it, (b) be >= ``min_gap_frac`` of the full
    value range wide, (c) save >= 25% of the axis height after collapse
    (only the excess gap beyond ``min_gap_frac`` counts as savings).

    Parameters
    ----------
    values : np.ndarray
        Finite plotted y values (post transform).
    min_gap_frac : float, optional
        Minimum gap width as a fraction of the value range for the gap
        to even be considered a break candidate. Default ``0.10``.

    Returns
    -------
    Optional[Tuple[float, float]]
        ``(lo, hi)`` bounds of the empty band to collapse, or ``None``
        when no gap qualifies (then a plain continuous axis is better).
    """
    v = np.sort(np.asarray(values, dtype=float))
    v = v[np.isfinite(v)]
    if len(v) < 20:
        return None
    span = v[-1] - v[0]
    if span <= 0:
        return None

    tail_max_frac = 0.10
    min_saving_frac = 0.25
    n = len(v)
    best: Optional[Tuple[float, float]] = None
    best_saving = 0.0
    # candidate gaps: each consecutive-pair gap that leaves <=10% of
    # points on one side
    for i in range(1, n):
        gap = v[i] - v[i - 1]
        if gap < min_gap_frac * span:
            continue
        # smaller side size (points isolated beyond the gap)
        smaller = min(i, n - i)
        if smaller / n > tail_max_frac:
            continue
        saving = (gap - min_gap_frac * span) / span  # collapse only the excess
        if saving > best_saving and saving >= min_saving_frac:
            best_saving = saving
            best = (float(v[i - 1]), float(v[i]))
    if best is not None:
        logger.info(
            f"Auto y-axis break detected: collapsing empty band "
            f"[{best[0]:.3g}, {best[1]:.3g}] (saves {best_saving:.0%} of axis height)"
        )
    return best


def dodge_labels(
    ax: "plt.Axes",
    items: List[Tuple[float, float, str]],
    fontsize: float = 8.0,
    push_pt: float = 2.0,
    frame_pad_pt: float = 4.0,
) -> List["plt.Text"]:
    """Draw point labels with a short leader line and global collision push.

    Each label is a plain text artist placed above its anchor point,
    plus a separate thin leader line segment (point -> label base) in
    data coordinates. Using plain texts (not arrow annotations) is
    required because ``Annotation.set_position`` is a no-op for drawn
    arrow annotations, which would break collision dodging. Colliding
    labels are pushed vertically in a single bottom-up sweep over
    x-overlapping clusters; a final global pairwise check pushes any
    remaining overlaps (labels from different clusters can collide
    after their own adjustments).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes (data coordinates).
    items : List[Tuple[float, float, str]]
        ``(x, y, text)`` entries in data coordinates.
    fontsize : float, optional
        Label font size in points. Default is ``8.0``.
    push_pt : float, optional
        Vertical push step in points. Default is ``2.0``.

    Returns
    -------
    List[plt.Text]
        The created label artists (leaders are redrawn by
        :func:`sync_leaders`; callers can remove texts to rebuild).
    """
    if not items:
        return []
    fig = ax.figure
    fig.canvas.draw()  # resolve positions before measuring extents
    renderer = fig.canvas.get_renderer()
    leader_gap_pt = 12.0  # label base sits this far above the anchor
    texts: List["plt.Text"] = []
    leaders: List["plt.Line2D"] = []
    # points->data scale for the vertical offset
    y0, y1 = ax.get_ylim()
    ax_bb = ax.get_window_extent(renderer=renderer)
    data_per_pt_y = (y1 - y0) / max(ax_bb.y1 - ax_bb.y0, 1.0)
    tick_pt = 3.0  # short leader tick length below the label
    # keep labels inside the axes horizontally: labels near the right
    # edge switch to right-align (text extends leftward), near the left
    # edge to left-align (text extends rightward)
    x_lo, x_hi = ax.get_xlim()
    x_span = max(x_hi - x_lo, 1.0)
    fig.canvas.draw()
    renderer0 = fig.canvas.get_renderer()
    ax_bb0 = ax.get_window_extent(renderer=renderer0)
    data_per_px_x = x_span / max(ax_bb0.x1 - ax_bb0.x0, 1.0)
    for x, y, label in items:
        ty = y + leader_gap_pt * data_per_pt_y
        # approximate label half-width in data units for alignment choice
        half_w = len(label) * fontsize * 0.3 * data_per_px_x
        ha = "center"
        if x + half_w > x_hi:
            ha = "right"
        elif x - half_w < x_lo:
            ha = "left"
        text = ax.text(
            x, ty, label,
            ha=ha, va="bottom", fontsize=fontsize,
            color="#1A1A1A", clip_on=True, zorder=6,
        )
        # short tick directly under the label (moves with it), NOT a long
        # line to the anchor -- long verticals overlap each other and
        # cross other labels when labels get pushed apart
        leader, = ax.plot(
            [x, x], [ty - tick_pt * data_per_pt_y, ty],
            lw=0.6, color="#555555", zorder=5, clip_on=True,
        )
        texts.append(text)
        leaders.append(leader)

    px_per_pt = fig.dpi / 72.0
    gap_px = push_pt * px_per_pt

    def sync_leaders() -> None:
        # keep each tick glued to its (possibly pushed) label base
        for text, leader in zip(texts, leaders):
            lx, ly = text.get_position()
            leader.set_data([lx, lx], [ly - tick_pt * data_per_pt_y, ly])

    def shift(idx: int, delta_px: float) -> None:
        x_pos, y_pos = texts[idx].get_position()
        texts[idx].set_position((x_pos, y_pos + delta_px * data_per_pt_y))
        sync_leaders()

    def boxes() -> List["matplotlib.transforms.Bbox"]:
        fig.canvas.draw()
        return [t.get_window_extent(renderer=renderer) for t in texts]

    base_boxes = boxes()
    offsets = [0.0] * len(texts)

    # Pass 1: deterministic stacking inside x-overlap clusters (sweep by
    # x0, split when the next label starts beyond the current x span).
    order = sorted(range(len(texts)), key=lambda i: base_boxes[i].x0)
    x_margin = 10.0
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

    frame_pad_px = frame_pad_pt * px_per_pt
    axes_top = ax.get_window_extent(renderer=renderer).y1 - frame_pad_px
    for cluster in clusters:
        cluster.sort(key=lambda i: base_boxes[i].y0)  # bottom-most first
        prev_top = None
        for i in cluster:
            b = base_boxes[i]
            y0, y1 = b.y0 + offsets[i], b.y1 + offsets[i]
            if prev_top is not None and y0 < prev_top + gap_px:
                offsets[i] += (prev_top + gap_px) - y0
                y0, y1 = b.y0 + offsets[i], b.y1 + offsets[i]
            prev_top = y1
        if prev_top > axes_top:
            # stack would overflow the axes top: rebuild downward,
            # clamping the topmost label inside the axes first
            cluster.sort(key=lambda i: base_boxes[i].y0, reverse=True)
            prev_bottom = None
            for i in cluster:
                b = base_boxes[i]
                y0, y1 = b.y0 + offsets[i], b.y1 + offsets[i]
                if y1 > axes_top:
                    offsets[i] -= y1 - axes_top
                    y0, y1 = b.y0 + offsets[i], b.y1 + offsets[i]
                if prev_bottom is not None and y1 > prev_bottom - gap_px:
                    offsets[i] -= y1 - (prev_bottom - gap_px)
                    y0, y1 = b.y0 + offsets[i], b.y1 + offsets[i]
                prev_bottom = y0

    for idx, text in enumerate(texts):
        base_x, base_y = text.get_position()
        text.set_position((base_x, base_y + offsets[idx] * data_per_pt_y))
    sync_leaders()

    # Pass 2: global pairwise collision sweep (catches cross-cluster
    # overlaps created when different clusters adjusted against each
    # other or the y ranks differ from x ranks). Pushes are boundary
    # aware: when pushing the upper label up would cross the axes top,
    # the lower label is pushed down instead (proportionally when both
    # are tight) -- a naive always-up push plus a clamp step oscillates
    # and never converges.
    axes_bb = ax.get_window_extent(renderer=renderer)
    pad_px = frame_pad_pt * px_per_pt
    axes_top, axes_bot = axes_bb.y1 - pad_px, axes_bb.y0 + pad_px

    for _ in range(60):
        cur_boxes = boxes()
        # fixed obstacles drawn earlier (boundary count notes, hline
        # labels) must not be overlapped by pushed labels
        obstacles = [
            t.get_window_extent(renderer=renderer)
            for t in ax.texts
            if t not in texts and t.get_visible()
        ]
        moved = False
        for a in range(len(texts)):
            for b_idx in range(a + 1, len(texts)):
                ba, bb = cur_boxes[a], cur_boxes[b_idx]
                if ba.x1 < bb.x0 or bb.x1 < ba.x0:
                    continue  # no x overlap -> cannot collide
                if ba.y1 < bb.y0 or bb.y1 < ba.y0:
                    continue  # already separated vertically
                upper, lower = (a, b_idx) if ba.y0 > bb.y0 else (b_idx, a)
                need = cur_boxes[lower].y1 + gap_px - cur_boxes[upper].y0
                if need <= 0:
                    continue
                room_up = axes_top - cur_boxes[upper].y1
                room_dn = cur_boxes[lower].y0 - axes_bot
                if need <= room_up:
                    shift(upper, need)
                elif need <= room_dn:
                    shift(lower, -need)
                elif room_up > 0 and room_dn > 0:
                    # split the push between both sides, never past a bound
                    shift(upper, room_up)
                    shift(lower, -(need - room_up))
                elif room_dn > 0:
                    shift(lower, -min(need, room_dn))
                elif room_up > 0:
                    shift(upper, min(need, room_up))
                else:
                    continue  # no room anywhere: leave as-is (best effort)
                moved = True
                cur_boxes = boxes()
        # label vs fixed obstacle: push the label away, respecting the
        # axes bounds (choose the side with room, else the roomier side)
        for a in range(len(texts)):
            ba = cur_boxes[a]
            for ob in obstacles:
                if ba.x1 < ob.x0 or ob.x1 < ba.x0:
                    continue
                if ba.y1 < ob.y0 or ob.y1 < ba.y0:
                    continue
                need_up = ob.y1 + gap_px - ba.y0
                need_dn = ba.y1 - (ob.y0 - gap_px)
                room_up = axes_top - ba.y1
                room_dn = ba.y0 - axes_bot
                if need_dn <= room_dn:
                    shift(a, -need_dn)
                elif need_up <= room_up:
                    shift(a, need_up)
                elif room_dn >= room_up:
                    shift(a, -room_dn)
                else:
                    shift(a, room_up)
                moved = True
                cur_boxes = boxes()
        if not moved:
            break
    logger.debug(f"Label dodge finished for {len(items)} labels in {len(clusters)} clusters")
    return texts


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
    ylim: Optional[Tuple[float, float]] = None,
    y_break: Optional[Tuple[float, float]] = None,
    y_break_auto: bool = False,
    signif_thresh: Optional[float] = None,
    size_smin: Optional[float] = None,
    size_smax: Optional[float] = None,
    size_name: Optional[str] = None,
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
    ylim : Tuple[float, float], optional
        Fixed y axis limits ``(ymin, ymax)``. Loci beyond the limits
        are drawn as hollow markers pinned to the boundary (clipped,
        not dropped) with a count note at each boundary. Default is
        ``None`` (matplotlib autoscale).

    Returns
    -------
    None
    """
    if image_formats is None:
        image_formats = ["png"]

    # --- broken y axis: detect and collapse an empty band ------------
    break_lo = break_hi = None
    if y_break_auto:
        detected = detect_y_break(loci["value"].to_numpy())
        if detected is not None:
            y_break = detected
    if y_break is not None:
        break_lo, break_hi = float(y_break[0]), float(y_break[1])
        if break_lo >= break_hi:
            raise ValueError(f"y_break must be (lo, hi) with lo < hi, got {y_break}")
        n_above = int((loci["value"] > break_hi).sum())
        n_below = int((loci["value"] < break_lo).sum())
        if n_above == 0 or n_below == 0:
            logger.warning(
                f"y_break [{break_lo:g}, {break_hi:g}] does not separate two "
                "non-empty groups; ignoring it"
            )
            break_lo = break_hi = None
        else:
            collapse = break_hi - break_lo
            def fold(v: float) -> float:
                if v > break_hi:
                    return v - collapse
                return v
            loci = loci.copy()
            loci["value"] = loci["value"].map(fold)
            if "signif" in loci.columns:
                pass  # significance already folded into is_burst earlier
            # fold the same band in fixed ylim / hlines if given
            if ylim is not None:
                ylim = (
                    float(ylim[0]) if ylim[0] <= break_lo else float(ylim[0]) - collapse,
                    float(ylim[1]) if ylim[1] <= break_lo else float(ylim[1]) - collapse,
                )
            if hlines:
                folded: List[float] = []
                for lv in hlines:
                    if lv <= break_lo or lv >= break_hi:
                        folded.append(lv if lv <= break_lo else lv - collapse)
                hlines = folded
            logger.info(
                f"y axis broken at [{break_lo:g}, {break_hi:g}]: "
                f"{n_below} loci below, {n_above} above"
            )

    if ylim is not None:
        y_min, y_max = float(ylim[0]), float(ylim[1])
        if y_min >= y_max:
            raise ValueError(f"ylim must be (ymin, ymax) with ymin < ymax, got {ylim}")
        out_of_range = loci[(loci["value"] < y_min) | (loci["value"] > y_max)]
        loci = loci.copy()
        if len(out_of_range):
            logger.info(
                f"y axis clipped to [{y_min:g}, {y_max:g}]: {len(out_of_range)} loci "
                "pinned to the boundary as hollow markers"
            )
            loci["value"] = loci["value"].clip(lower=y_min, upper=y_max)
            loci["is_pinned"] = loci.index.isin(out_of_range.index)
        else:
            loci["is_pinned"] = False
    else:
        y_min = y_max = 0.0  # sentinel; only referenced under `ylim is not None`
        loci = loci.copy()
        loci["is_pinned"] = False

    n_highlight_families = loci.loc[loci["is_highlight"], "family"].nunique() if "is_highlight" in loci else 0
    has_legend = bool(loci["is_burst"].any()) or n_highlight_families > 0
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax_leg = ax

    chrom_names = resolve_chrom_order(list(loci["chrom"].unique()), chrom_order)
    # marker area from the abundance column (e.g. baseMean) when present;
    # default range keeps the original point sizes: background points
    # scale up to the burst marker size (highlight_size + 1)
    if "size_value" in loci.columns:
        s_lo = point_size if size_smin is None else size_smin
        s_hi = highlight_size + 1.0 if size_smax is None else size_smax
        loci["marker_s"] = size_scale(loci["size_value"], s_lo, s_hi)
    else:
        loci["marker_s"] = point_size
    background = loci[~loci["is_highlight"] & ~loci["is_burst"]]
    # When significance is known, non-significant background points are
    # gray so the significant ones stand out by color, not just direction.
    if "signif" in background.columns and signif_thresh is not None:
        bg_sig_mask = background["signif"] <= signif_thresh
    else:
        bg_sig_mask = None
    if max_points > 0 and len(background) > max_points:
        logger.info(
            f"Background loci downsampled for rendering: {len(background)} -> {max_points} "
            "(burst/highlight points kept in full)"
        )
        background = background.sample(n=max_points, random_state=seed)
    for idx, chrom in enumerate(chrom_names):
        sub = background[background["chrom"] == chrom]
        # Down-regulated points (value < 0) get a solid blue so they stay
        # readable against the pale alternating chromosome colors;
        # non-significant points (when testable) are drawn gray.
        if bg_sig_mask is not None:
            sub_mask = bg_sig_mask.loc[sub.index]
            sig, nonsig = sub[sub_mask], sub[~sub_mask]
        else:
            sig, nonsig = sub, sub.iloc[0:0]
        if len(nonsig):
            ax.scatter(
                nonsig["xpos"], nonsig["value"], s=nonsig["marker_s"], c="#B0B0B0",
                alpha=min(point_alpha + 0.2, 1.0), linewidths=0,
                rasterized=True, zorder=2,
            )
        chrom_color = CHROM_COLOR_ODD if idx % 2 == 0 else CHROM_COLOR_EVEN
        down = sig[sig["value"] < 0]
        up = sig[sig["value"] >= 0]
        ax.scatter(
            down["xpos"], down["value"], s=down["marker_s"], c=DOWN_POINT_COLOR,
            alpha=min(point_alpha + 0.25, 1.0), linewidths=0, rasterized=True, zorder=2,
        )
        ax.scatter(
            up["xpos"], up["value"], s=up["marker_s"], c=chrom_color,
            alpha=point_alpha, linewidths=0, rasterized=True, zorder=2,
        )

    legend_handles: List[object] = []
    # background point legend (proxy markers -- the actual scatters span
    # many chromosome colors / sizes)
    if bg_sig_mask is not None:
        legend_handles.extend(
            [
                plt.Line2D([], [], marker="o", ls="", markerfacecolor=BURST_COLOR_UP,
                           markeredgecolor="none", markersize=6, label="up (significant)"),
                plt.Line2D([], [], marker="o", ls="", markerfacecolor=DOWN_POINT_COLOR,
                           markeredgecolor="none", markersize=6, label="down (significant)"),
                plt.Line2D([], [], marker="o", ls="", markerfacecolor="#B0B0B0",
                           markeredgecolor="none", markersize=6, label="not significant"),
            ]
        )
    highlight_family_names = sorted(loci.loc[loci["is_highlight"], "family"].astype(str).unique())
    # significance encoded within family color: significant loci keep full
    # opacity with a thin white edge, non-significant ones fade -- avoids
    # the ugly dark-ring look while keeping significance readable
    hl_sig = None
    if "signif" in loci.columns and signif_thresh is not None:
        hl_sig = (loci["signif"] <= signif_thresh).fillna(False)
    for idx, fam_name in enumerate(highlight_family_names):
        fam_mask = loci["is_highlight"] & (loci["family"].astype(str) == fam_name)
        sub = loci[fam_mask]
        if hl_sig is not None:
            sig_mask = hl_sig[fam_mask]
            sub_sig, sub_nos = sub[sig_mask], sub[~sig_mask]
        else:
            sub_sig, sub_nos = sub.iloc[0:0], sub
        fam_color = HIGHLIGHT_PALETTE[idx % len(HIGHLIGHT_PALETTE)]
        if len(sub_nos):
            ax.scatter(
                sub_nos["xpos"], sub_nos["value"], s=sub_nos["marker_s"],
                c=fam_color, alpha=0.25, linewidths=0, zorder=3,
            )
        if len(sub_sig):
            ax.scatter(
                sub_sig["xpos"], sub_sig["value"], s=sub_sig["marker_s"],
                c=fam_color, alpha=0.95, linewidths=0.6, edgecolors="white",
                zorder=4,
            )
        handle = plt.Line2D(
            [], [], marker="o", ls="", markerfacecolor=fam_color,
            markeredgecolor="white", markeredgewidth=0.6, markersize=6,
            label=f"{fam_name} (opaque = significant, pale = not)",
        )
        legend_handles.append(handle)

    if loci["is_burst"].any():
        burst = loci[loci["is_burst"] & ~loci["is_highlight"]]
        up = burst[burst["value"] >= 0]
        down = burst[burst["value"] < 0]
        if len(up):
            ax.scatter(
                up["xpos"], up["value"], s=up["marker_s"], c=BURST_COLOR_UP,
                alpha=0.95, linewidths=0, zorder=4,
            )
        if len(down):
            ax.scatter(
                down["xpos"], down["value"], s=down["marker_s"], c=BURST_COLOR_DOWN,
                alpha=0.95, linewidths=0, zorder=4,
            )

    if ylim is not None and loci["is_pinned"].any():
        pinned = loci[loci["is_pinned"]]
        ax.scatter(
            pinned["xpos"], pinned["value"], s=highlight_size, facecolors="none",
            edgecolors="#333333", linewidths=0.8, zorder=5,
        )
        n_below = int((loci["value"] <= y_min).sum())
        n_above = int((loci["value"] >= y_max).sum())
        for note_y, note_n, va in (
            (y_min, n_below, "bottom"),
            (y_max, n_above, "top"),
        ):
            if note_n > 0:
                ax.annotate(
                    f"{note_n} loci beyond limit",
                    xy=(0.0, note_y),
                    xycoords=("axes fraction", "data"),
                    ha="left", va=va, fontsize=8, color="#333333",
                    annotation_clip=False, zorder=6,
                )

    if ylim is not None:
        ax.set_ylim(y_min, y_max)

    centers = loci.groupby("chrom")["chrom_center"].first()
    ax.set_ylabel(ylabel)

    if break_lo is not None:
        # broken-axis marker: a diagonal slash pair at the fold position
        # on the y axis spine (data above break_hi was shifted down by
        # the collapse width, so the fold sits at break_lo)
        marker_kw = dict(
            transform=ax.get_yaxis_transform(), clip_on=False, zorder=10,
            linewidth=1.2, color="black",
        )
        ax.plot(
            [-0.012, 0.012], [break_lo - 0.006, break_lo + 0.006],
            solid_capstyle="butt", **marker_kw,
        )
        ax.plot(
            [-0.012, 0.012], [break_lo - 0.002, break_lo + 0.010],
            solid_capstyle="butt", **marker_kw,
        )
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

    label_items: List[Tuple[float, float, str]] = []
    if label_top > 0:
        # Label with the locus id's first ':' segment (the rmsk
        # transcript_id, e.g. AmnL2-1_dup1083), not `name` (repName,
        # shared by thousands of loci).
        label_col = "locus_id" if "locus_id" in loci.columns else "name"
        marked = loci[loci["is_burst"] | loci["is_highlight"]]
        top = marked.reindex(
            marked["value"].abs().sort_values(ascending=False).index
        ).head(label_top)
        for _, row in top.iterrows():
            label_items.append(
                (float(row["xpos"]), float(row["value"]), str(row[label_col]).split(":")[0])
            )
    if label_top > 0:
        dodge_labels(ax, label_items)

    if ax_leg is not None:
        # fit the legend within the axes width: try ncol columns and
        # reduce until the rendered legend box is no wider than the axes
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        ax_width_px = ax_leg.get_window_extent(renderer=renderer).width
        leg = None
        for ncol in range(max(1, len(legend_handles)), 0, -1):
            leg = ax_leg.legend(
                handles=legend_handles,
                loc="upper left",
                bbox_to_anchor=(0.0, -0.16),
                ncol=ncol,
                frameon=False,
                fontsize=9,
            )
            fig.canvas.draw()
            if leg.get_window_extent(renderer=renderer).width <= ax_width_px:
                break
        if leg is not None:
            logger.debug(
                f"Legend rendered at ncol={getattr(leg, '_ncols', '?')} within axes width"
            )

    fig.align_ylabels()
    for fmt in image_formats:
        out_path = f"{out_base}.{fmt}"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
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
    parser.add_argument("--values", default=None,
                        help="Per-locus value TSV (see module docstring). Required unless "
                             "--deseq mode is used.")
    parser.add_argument("--deseq", default=None,
                        help="DESeq2 result TSV (e.g. TElocal *_TE.tsv) with per-sample count "
                             "columns; with --group/--control/--experiment the per-locus TPM, "
                             "delta TPM and mean TPM are computed internally from it.")
    parser.add_argument("--group", default=None,
                        help="group.tsv (sample/group) mapping count columns to conditions.")
    parser.add_argument("--control", default=None, help="Control group name in --group.")
    parser.add_argument("--experiment", default=None, help="Experiment group name in --group.")
    parser.add_argument("--loci-tsv", default=None,
                        help="Locus coordinate table (locus_id/chrom/start/end), e.g. dumped "
                             "from a TElocal .locInd file.")
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
    parser.add_argument("--ylim", default=None, metavar="YMIN,YMAX",
                        help="Fixed y axis limits. Loci beyond the limits are drawn "
                             "pinned to the boundary as hollow markers.")
    parser.add_argument("--y-break", default=None, metavar="LO,HI",
                        help="Collapse the empty y band [LO,HI] (manual broken axis).")
    parser.add_argument("--y-break-auto", action="store_true",
                        help="Auto-detect a justified y-axis break: only when an empty "
                             "band isolates <=10%% of points as a tail and collapsing it "
                             "saves >=25%% of the axis height. Default: off.")
    parser.add_argument("--size-col", default=None,
                        help="Optional numeric column (e.g. baseMean) mapped to marker "
                             "size via log10(1+x) scaling.")
    parser.add_argument("--size-smin", type=float, default=None,
                        help="Minimum marker area for --size-col mapping "
                             "(default: the background point size).")
    parser.add_argument("--size-smax", type=float, default=None,
                        help="Maximum marker area for --size-col mapping "
                             "(default: the burst/highlight point size).")
    parser.add_argument("--only-significant", action="store_true",
                        help="Plot only loci passing --signif-col/--signif-thresh "
                             "(default: plot all loci).")
    parser.add_argument("--y-transform", choices=["log2", "log10", "none"], default="log2",
                        help="Y value transform. log2/log10 use sign(x)*log(1+|x|) "
                             "for signed values. Default: log2.")
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

    if args.deseq:
        # DESeq2-driven mode: TPM computed internally from the result
        # table's per-sample counts + group file + locus coordinates
        if not (args.group and args.control and args.experiment):
            raise ValueError("--deseq requires --group, --control and --experiment")
        loci = load_deseq_tpm(
            deseq_tsv=args.deseq,
            group_tsv=args.group,
            control=args.control,
            experiment=args.experiment,
            loci_tsv=args.loci_tsv,
            signif_col=args.signif_col or "padj",
        )
        # value column naming inside this mode
        if args.value_col != "value":
            logger.info(f"--deseq mode: --value-col '{args.value_col}' ignored, using delta TPM")
        if args.size_col is None:
            args.size_col = "mean_tpm"
            logger.info("--deseq mode: point size defaults to mean_tpm")
        if args.size_col == "mean_tpm" and "size_value" not in loci.columns:
            loci["size_value"] = loci["mean_tpm"]
    else:
        if not args.values:
            raise ValueError("--values is required unless --deseq mode is used")
        loci = load_locus_values(
            values_tsv=args.values,
            value_col=args.value_col,
            chrom_col=args.chrom_col,
            start_col=args.start_col,
            end_col=args.end_col,
            name_col=args.name_col,
            family_col=args.family_col,
            class_col=args.class_col,
            signif_col=args.signif_col,
            size_col=args.size_col,
        )
    loci = select_highlights(
        loci,
        highlight_families=args.highlight_families,
        burst_thresh=args.value_thresh,
        signif_thresh=args.signif_thresh,
        burst_sign=args.burst_sign,
    )
    if args.only_significant:
        if "signif" not in loci.columns or args.signif_thresh is None:
            raise ValueError(
                "--only-significant requires --signif-col and --signif-thresh"
            )
        n_before = len(loci)
        loci = loci[loci["signif"] <= args.signif_thresh].reset_index(drop=True)
        logger.info(f"--only-significant: kept {len(loci)}/{n_before} loci")
    loci = transform_value(loci, y_transform=args.y_transform)
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
        ylim=(float(args.ylim.split(",")[0]), float(args.ylim.split(",")[1]))
        if args.ylim else None,
        y_break=(float(args.y_break.split(",")[0]), float(args.y_break.split(",")[1]))
        if args.y_break else None,
        y_break_auto=args.y_break_auto,
        signif_thresh=args.signif_thresh,
        size_smin=args.size_smin,
        size_smax=args.size_smax,
        size_name=args.size_col,
    )


if __name__ == "__main__":
    main()
