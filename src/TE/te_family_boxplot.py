#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TE family expression distribution across sample groups (stages).

For each requested TE family pattern, plots the distribution of
normalized TE expression across the TE entries matching that family,
one box/violin per sample group (e.g. developmental stage). A burst
shows up as the whole distribution of a family shifting up in one
group.

Input data
----------
``--matrix`` : TSV with a TE id column and one column per sample.

* TEcount combined tables (``<ref>_TEcount.tsv``): first column
  ``gene/TE`` holds ids like ``L1MD1:L1:LINE`` (``name:family:class``)
  or gene ids (``ENSG...``); remaining columns are raw counts.
  Gene rows are ignored by family matching but DO count toward the
  CPM library size.
* DESeq2 result tables (``*.TEcount_TE.tsv``): first column is the TE
  id; the stat columns ``baseMean/log2FoldChange/lfcSE/stat/pvalue/
  padj`` are auto-excluded from the sample columns; the remaining
  columns are DESeq2-normalized counts (pass ``--normalize none``).

TE id parsing: ``name:family:class`` split on ``:`` -- field 0 =
name/subfamily (e.g. ``L1MD1``, ``MT2_Mm``), field 1 = family (e.g.
``L1``, ``ERVL``), field 2 = class (e.g. ``LINE``, ``LTR``).

``--annotation`` : TSV with header and columns ``sample`` and
``group`` (group = stage/condition). Optional column ``order`` gives
the group plotting order (smaller first). ``--group-order`` (comma
list) overrides it. Samples present in the matrix but missing from
the annotation are dropped with a warning; annotation rows without a
matrix column are ignored.

Examples
--------
python te_family_boxplot.py --matrix GRCm39_TEcount.tsv \\
    --annotation samples.tsv --family MT2_Mm --family L1Md --family IAP \\
    --normalize log2cpm --kind violin -o outdir -f png -f pdf
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

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

DESEQ2_STAT_COLS = ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
GROUP_PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#64B5CD", "#CCB974"]


def load_matrix(
    matrix_tsv: str,
    id_col: Optional[str] = None,
    exclude_cols: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """Load the TE expression matrix and identify sample columns.

    Parameters
    ----------
    matrix_tsv : str
        Path to the count/normalized matrix TSV (see module docstring).
    id_col : str, optional
        TE id column name. Default is the first column.
    exclude_cols : List[str], optional
        Additional non-sample columns to exclude.

    Returns
    -------
    Tuple[pd.DataFrame, List[str]]
        ``(matrix, sample_columns)`` where ``matrix`` is indexed by
        the raw TE id string and holds only numeric sample columns.

    Raises
    ------
    ValueError
        If no sample columns remain or values are not numeric.
    """
    raw = pd.read_csv(matrix_tsv, sep="\t")
    if id_col is None:
        id_col = raw.columns[0]
    if id_col not in raw.columns:
        raise ValueError(f"id column '{id_col}' not in {matrix_tsv}: {list(raw.columns)}")

    excluded = set(DESEQ2_STAT_COLS) | set(exclude_cols or [])
    dropped = sorted(excluded & set(raw.columns))
    if dropped:
        logger.info(f"Non-sample columns excluded from matrix: {dropped}")

    sample_cols = [c for c in raw.columns if c != id_col and c not in excluded]
    numeric_cols = [
        c for c in sample_cols
        if pd.api.types.is_numeric_dtype(pd.to_numeric(raw[c], errors="coerce"))
    ]
    non_numeric = sorted(set(sample_cols) - set(numeric_cols))
    if non_numeric:
        logger.warning(f"Non-numeric columns excluded from samples: {non_numeric}")
    if not numeric_cols:
        raise ValueError(f"No numeric sample columns found in {matrix_tsv}")

    matrix = raw.set_index(id_col)[numeric_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    logger.info(f"Loaded matrix {matrix.shape[0]} TE entries x {matrix.shape[1]} samples from {matrix_tsv}")
    return matrix, numeric_cols


def load_annotation(
    annotation_tsv: str,
    sample_cols: List[str],
    group_order: Optional[str] = None,
) -> "pd.DataFrame":
    """Load the sample-to-group annotation and the group order.

    Parameters
    ----------
    annotation_tsv : str
        TSV with header columns ``sample``, ``group`` and optionally
        ``order`` (int; smaller plots first).
    sample_cols : List[str]
        Sample columns of the matrix (used for validation).
    group_order : str, optional
        Comma-separated group order overriding the ``order`` column.

    Returns
    -------
    pd.DataFrame
        Columns ``sample``, ``group`` restricted to matrix samples,
        with categorical ``group`` in plotting order.

    Raises
    ------
    ValueError
        If required columns are missing or no annotated sample
        intersects the matrix.
    """
    annot = pd.read_csv(annotation_tsv, sep="\t")
    for col in ("sample", "group"):
        if col not in annot.columns:
            raise ValueError(f"annotation column '{col}' missing in {annotation_tsv}: {list(annot.columns)}")

    annot = annot[annot["sample"].isin(sample_cols)].copy()
    missing = sorted(set(sample_cols) - set(annot["sample"]))
    if missing:
        logger.warning(f"Matrix samples without annotation dropped: {missing}")
    if annot.empty:
        raise ValueError(f"No annotated samples overlap the matrix columns: {sample_cols}")

    if group_order:
        order = [g.strip() for g in group_order.split(",") if g.strip()]
    elif "order" in annot.columns:
        order = (
            annot[["group", "order"]]
            .drop_duplicates()
            .sort_values("order")["group"]
            .astype(str)
            .tolist()
        )
    else:
        order = sorted(annot["group"].astype(str).unique())

    unseen = sorted(set(annot["group"].astype(str)) - set(order))
    if unseen:
        logger.warning(f"Groups not in group order appended at end: {unseen}")
        order.extend(unseen)
    annot["group"] = pd.Categorical(annot["group"].astype(str), categories=order, ordered=True)
    logger.info(f"Group order: {order}")
    return annot[["sample", "group"]].reset_index(drop=True)


def parse_te_id(te_id: str) -> Dict[str, str]:
    """Split a ``name:family:class`` TE id into its fields.

    Parameters
    ----------
    te_id : str
        TE id such as ``L1MD1:L1:LINE``. Ids without ``:`` (gene ids)
        return all fields equal to the id.

    Returns
    -------
    Dict[str, str]
        Keys ``name``, ``family``, ``class``.
    """
    fields = str(te_id).split(":")
    if len(fields) >= 3:
        return {"name": fields[0], "family": fields[1], "class": fields[2]}
    return {"name": str(te_id), "family": str(te_id), "class": str(te_id)}


def select_family_entries(
    matrix: pd.DataFrame,
    family_pattern: str,
    match_field: str = "any",
    match_mode: str = "substring",
) -> pd.DataFrame:
    """Return the matrix rows matching one family pattern.

    Parameters
    ----------
    matrix : pd.DataFrame
        Matrix indexed by TE id (see :func:`load_matrix`).
    family_pattern : str
        Pattern to match.
    match_field : str, optional
        ``"name"``, ``"family"``, ``"class"`` or ``"any"`` (match
        against the full id string). Default is ``"any"``.
    match_mode : str, optional
        ``"substring"`` (case-insensitive), ``"regex"`` or
        ``"exact"`` (case-insensitive on the selected field).
        Default is ``"substring"``.

    Returns
    -------
    pd.DataFrame
        Matching rows of ``matrix`` (possibly empty).

    Raises
    ------
    ValueError
        If ``match_field`` or ``match_mode`` is invalid.
    """
    if match_field not in ("name", "family", "class", "any"):
        raise ValueError(f"match_field must be name/family/class/any, got {match_field!r}")
    if match_mode not in ("substring", "regex", "exact"):
        raise ValueError(f"match_mode must be substring/regex/exact, got {match_mode!r}")

    if match_field == "any":
        haystack = matrix.index.astype(str).str.lower()
    else:
        haystack = pd.Index(
            [parse_te_id(i)[match_field].lower() for i in matrix.index], name=matrix.index.name
        )

    pattern = family_pattern.lower() if match_mode != "regex" else family_pattern
    if match_mode == "exact":
        mask = haystack == pattern
    elif match_mode == "regex":
        mask = haystack.str.contains(pattern, regex=True, na=False)
    else:
        mask = haystack.str.contains(pattern, regex=False, na=False)
    return matrix[np.asarray(mask, dtype=bool)]


def normalize_matrix(
    matrix: pd.DataFrame,
    method: str = "log2cpm",
) -> pd.DataFrame:
    """Normalize raw counts to CPM / log2(CPM+1).

    Parameters
    ----------
    matrix : pd.DataFrame
        Raw count matrix (all gene + TE rows; library sizes are column
        sums over everything given).
    method : str, optional
        ``"log2cpm"`` (default), ``"cpm"`` or ``"none"``.

    Returns
    -------
    pd.DataFrame
        Normalized matrix.

    Raises
    ------
    ValueError
        If ``method`` is invalid or a library size is zero.
    """
    if method not in ("log2cpm", "cpm", "none"):
        raise ValueError(f"normalize method must be log2cpm/cpm/none, got {method!r}")
    if method == "none":
        return matrix.copy()

    lib_sizes = matrix.sum(axis=0)
    zero_libs = lib_sizes[lib_sizes <= 0]
    if not zero_libs.empty:
        raise ValueError(f"Zero library size for samples: {list(zero_libs.index)}")
    cpm = matrix.div(lib_sizes, axis=1) * 1e6
    if method == "cpm":
        return cpm
    return np.log2(cpm + 1.0)


def collect_group_values(
    family_matrix: pd.DataFrame,
    annot: "pd.DataFrame",
    value_scope: str = "group_mean_per_entry",
) -> Dict[str, np.ndarray]:
    """Aggregate per-group value distributions for one family.

    Parameters
    ----------
    family_matrix : pd.DataFrame
        Normalized values for the entries of one family (entries x
        samples).
    annot : pd.DataFrame
        Sample/group annotation (see :func:`load_annotation`).
    value_scope : str, optional
        ``"group_mean_per_entry"`` (default): each entry contributes
        its mean across the group's replicates -> distribution across
        entries. ``"all_samples"``: pool every (entry, sample) value.

    Returns
    -------
    Dict[str, np.ndarray]
        Group name -> 1-D array of values (empty array if the group
        has no samples).

    Raises
    ------
    ValueError
        If ``value_scope`` is invalid.
    """
    if value_scope not in ("group_mean_per_entry", "all_samples"):
        raise ValueError(f"value_scope must be group_mean_per_entry/all_samples, got {value_scope!r}")

    out: Dict[str, np.ndarray] = {}
    for group in annot["group"].cat.categories:
        samples = annot.loc[annot["group"] == group, "sample"].tolist()
        samples = [s for s in samples if s in family_matrix.columns]
        if not samples:
            out[str(group)] = np.array([], dtype=float)
            continue
        sub = family_matrix[samples]
        if value_scope == "all_samples":
            out[str(group)] = sub.to_numpy().ravel()
        else:
            out[str(group)] = sub.mean(axis=1).to_numpy()
    return out


def plot_family_distributions(
    family_data: List[Tuple[str, pd.DataFrame]],
    annot: "pd.DataFrame",
    out_base: str,
    image_formats: Optional[List[PlotFormat]] = None,
    kind: str = "box",
    value_scope: str = "group_mean_per_entry",
    ncols: int = 3,
    figsize_unit: Tuple[float, float] = (4.2, 3.6),
    dpi: int = 300,
    ylabel: str = "log2(CPM+1)",
    show_points: bool = True,
    point_jitter: float = 0.12,
    point_size: float = 2.5,
    point_alpha: float = 0.5,
    seed: int = 0,
    test: str = "none",
    test_ref_group: Optional[str] = None,
    x_rotation: float = 30.0,
    box_width: float = 0.6,
    violin_bw: float = 0.3,
    row_spacing: float = 0.6,
    stat_headroom: float = 2.4,
) -> None:
    """Render the per-group distribution panels, one panel per family.

    Parameters
    ----------
    family_data : List[Tuple[str, pd.DataFrame]]
        ``(family_pattern, family_matrix)`` in requested order;
        an empty matrix renders a labeled placeholder panel.
    annot : pd.DataFrame
        Sample/group annotation with ordered categorical ``group``.
    out_base : str
        Output path without extension.
    image_formats : List[str], optional
        Save formats. Default is ``["png"]``.
    kind : str, optional
        ``"box"``, ``"violin"`` or ``"boxen"``-like box with points.
        Default is ``"box"``.
    value_scope : str, optional
        See :func:`collect_group_values`. Default is
        ``"group_mean_per_entry"``.
    ncols : int, optional
        Panel columns. Default is ``3``.
    figsize_unit : Tuple[float, float], optional
        Size per panel in inches. Default is ``(4.2, 3.6)``.
    dpi : int, optional
        Output resolution. Default is ``300``.
    ylabel : str, optional
        Y axis label shown on the left panels.
    show_points : bool, optional
        Overlay jittered per-entry points. Default is ``True``.
    point_jitter : float, optional
        Half-width of the point jitter in x units. Default ``0.12``.
    point_size : float, optional
        Point size. Default is ``2.5``.
    point_alpha : float, optional
        Point alpha. Default is ``0.5``.
    seed : int, optional
        Jitter RNG seed. Default is ``0``.
    test : str, optional
        ``"none"`` (default) or ``"mannwhitney"``: annotate the
        one-sided-greater Mann-Whitney p of each group vs the
        reference group.
    test_ref_group : str, optional
        Reference group for ``test``. Default is the first group.
    x_rotation : float, optional
        Group tick label rotation in degrees. Default is ``30.0``.
    box_width : float, optional
        Box width. Default is ``0.6``.
    violin_bw : float, optional
        Violin kernel bandwidth. Default is ``0.3``.
    row_spacing : float, optional
        Vertical spacing between panel rows (matplotlib ``hspace``),
        reserved for the rotated x tick labels of the row above.
        Default is ``0.6``.
    stat_headroom : float, optional
        Top margin as a multiple of the data range, reserved for the
        ``test`` annotations inside the axes. Default is ``2.4``.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If ``kind`` or ``test`` is invalid.
    """
    if image_formats is None:
        image_formats = ["png"]
    if kind not in ("box", "violin"):
        raise ValueError(f"kind must be box/violin, got {kind!r}")
    if test not in ("none", "mannwhitney"):
        raise ValueError(f"test must be none/mannwhitney, got {test!r}")

    groups = [str(g) for g in annot["group"].cat.categories]
    n_panels = len(family_data)
    ncols = max(1, min(ncols, n_panels))
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_unit[0] * ncols, figsize_unit[1] * nrows),
        squeeze=False,
        gridspec_kw={"hspace": row_spacing},
    )

    rng = np.random.default_rng(seed)
    for idx, (pattern, fam_matrix) in enumerate(family_data):
        ax = axes[idx // ncols][idx % ncols]
        data = collect_group_values(fam_matrix, annot, value_scope=value_scope)
        non_empty = [g for g in groups if data[g].size > 0]

        if not non_empty:
            ax.text(
                0.5, 0.5, f"no matching TE\nfor '{pattern}'",
                ha="center", va="center", fontsize=10, color="#666666", transform=ax.transAxes,
            )
            ax.set_title(f"{pattern} (n=0 entries)", fontsize=11)
            ax.set_xticks(range(len(groups)))
            ax.set_xticklabels(groups, rotation=x_rotation, ha="right")
            continue

        positions = np.arange(len(groups), dtype=float)
        plot_data = [data[g] for g in groups]

        if kind == "violin":
            parts = ax.violinplot(
                plot_data, positions=positions, showmeans=False, showextrema=False,
                widths=box_width, bw_method=violin_bw,
            )
            for body in parts["bodies"]:
                body.set_facecolor("#A0C4E8")
                body.set_edgecolor("#4C72B0")
                body.set_alpha(0.75)
            for pos, values in zip(positions, plot_data):
                if values.size == 0:
                    continue
                ax.vlines(pos, np.percentile(values, 25), np.percentile(values, 75),
                          color="#2B4C7E", linewidth=3.0, zorder=3)
                ax.scatter([pos], [np.median(values)], s=18, color="white",
                           edgecolor="#2B4C7E", linewidth=1.2, zorder=4)
        else:
            boxes = ax.boxplot(
                plot_data, positions=positions, widths=box_width,
                patch_artist=True, showfliers=False,
                medianprops=dict(color="#2B4C7E", linewidth=1.4),
                whiskerprops=dict(color="#2B4C7E"),
                capprops=dict(color="#2B4C7E"),
            )
            for patch in boxes["boxes"]:
                patch.set_facecolor("#A0C4E8")
                patch.set_edgecolor("#2B4C7E")
                patch.set_alpha(0.85)

        if show_points:
            for pos, values in zip(positions, plot_data):
                if values.size == 0:
                    continue
                jitter = rng.uniform(-point_jitter, point_jitter, size=values.size)
                ax.scatter(
                    np.full(values.size, pos) + jitter, values,
                    s=point_size, color="#2B4C7E", alpha=point_alpha,
                    linewidths=0, zorder=5, rasterized=True,
                )

        if test == "mannwhitney":
            try:
                from scipy.stats import mannwhitneyu
            except ImportError as exc:
                logger.warning(f"scipy not available, skipping stats annotation: {exc}")
            else:
                ref = test_ref_group if test_ref_group in data else groups[0]
                ref_values = data[ref]
                y_top = max((v.max() for v in plot_data if v.size), default=1.0)
                y_bot = min((v.min() for v in plot_data if v.size), default=0.0)
                span = (y_top - y_bot) if y_top > y_bot else 1.0
                pad = span * 0.06
                for pos, group in zip(positions, groups):
                    if group == ref or data[group].size == 0 or ref_values.size == 0:
                        continue
                    _, p_val = mannwhitneyu(data[group], ref_values, alternative="greater")
                    ax.text(
                        pos, y_top + pad * (1.0 + 1.2 * (pos % 2)), f"p={p_val:.2g}",
                        ha="center", va="bottom", fontsize=7, color="#444444", clip_on=True,
                    )
                # reserve in-axes headroom so p texts never touch the title
                ax.set_ylim(y_bot - span * 0.08, y_top + span * stat_headroom)

        n_entries = fam_matrix.shape[0]
        ax.set_title(f"{pattern} (n={n_entries} entries)", fontsize=11)
        ax.set_xticks(positions)
        ax.set_xticklabels(groups, rotation=x_rotation, ha="right")
        ax.set_xlim(-0.6, len(groups) - 0.4)

    # leftover panels: turn off
    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    # ylabel ownership: left-column panels only
    for row in range(nrows):
        axes[row][0].set_ylabel(ylabel)
        for col in range(1, ncols):
            axes[row][col].tick_params(axis="y", labelleft=True)

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
        description="TE family expression distributions across sample groups.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--matrix", required=True, help="TE count/normalized matrix TSV (see module docstring).")
    parser.add_argument("--id-col", default=None, help="TE id column name (default: first column).")
    parser.add_argument("--annotation", required=True, help="sample/group annotation TSV.")
    parser.add_argument("--group-order", default=None, help="Comma-separated group order.")
    parser.add_argument("--family", action="append", dest="families", metavar="PATTERN", required=True,
                        help="TE family pattern. Repeatable.")
    parser.add_argument("--match-field", choices=["name", "family", "class", "any"], default="any",
                        help="TE id field to match the pattern against.")
    parser.add_argument("--match-mode", choices=["substring", "regex", "exact"], default="substring",
                        help="Pattern matching mode.")
    parser.add_argument("--normalize", choices=["log2cpm", "cpm", "none"], default="log2cpm",
                        help="Count normalization; use 'none' for DESeq2-normalized sample columns.")
    parser.add_argument("--value-scope", choices=["group_mean_per_entry", "all_samples"],
                        default="group_mean_per_entry",
                        help="Per-entry group means (distribution across entries) or pooled samples.")
    parser.add_argument("--exclude-col", action="append", dest="exclude_cols", metavar="COL",
                        help="Extra non-sample column to exclude. Repeatable.")
    parser.add_argument("--kind", choices=["box", "violin"], default="box", help="Distribution glyph.")
    parser.add_argument("--test", choices=["none", "mannwhitney"], default="none",
                        help="Annotate Mann-Whitney p per group vs reference.")
    parser.add_argument("--test-ref-group", default=None, help="Reference group for --test.")
    parser.add_argument("--ncols", type=int, default=3, help="Panel columns.")
    parser.add_argument("--figsize-unit", default="4.2,3.6", help="Panel size in inches, 'w,h'.")
    parser.add_argument("--ylabel", default="log2(CPM+1)", help="Y axis label.")
    parser.add_argument("--no-points", action="store_true", help="Hide the jittered points.")
    parser.add_argument("--point-jitter", type=float, default=0.12, help="Point jitter half-width.")
    parser.add_argument("--point-size", type=float, default=2.5, help="Point size.")
    parser.add_argument("--point-alpha", type=float, default=0.5, help="Point alpha.")
    parser.add_argument("--seed", type=int, default=0, help="Jitter RNG seed.")
    parser.add_argument("--x-rotation", type=float, default=30.0, help="Group label rotation (degrees).")
    parser.add_argument("--box-width", type=float, default=0.6, help="Box/violin width.")
    parser.add_argument("--violin-bw", type=float, default=0.3, help="Violin kernel bandwidth.")
    parser.add_argument("--row-spacing", type=float, default=0.6,
                        help="Vertical spacing between panel rows (hspace).")
    parser.add_argument("--stat-headroom", type=float, default=2.4,
                        help="Top margin (x data range) reserved for --test annotations.")
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution.")
    parser.add_argument("-f", "--format", action="append", dest="formats", metavar="FMT",
                        help="Image format (png/pdf/svg/...). Repeatable. Default: png.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory.")
    parser.add_argument("--out-prefix", default="te_family_box", help="Output file prefix.")
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

    matrix, sample_cols = load_matrix(args.matrix, id_col=args.id_col, exclude_cols=args.exclude_cols)
    annot = load_annotation(args.annotation, sample_cols, group_order=args.group_order)
    normalized = normalize_matrix(matrix, method=args.normalize)

    family_data: List[Tuple[str, pd.DataFrame]] = []
    for pattern in args.families:
        fam = select_family_entries(
            normalized, pattern, match_field=args.match_field, match_mode=args.match_mode
        )
        logger.info(f"Family pattern '{pattern}': {fam.shape[0]} entries")
        family_data.append((pattern, fam))

    w, h = (float(v) for v in args.figsize_unit.split(","))
    plot_family_distributions(
        family_data=family_data,
        annot=annot,
        out_base=os.path.join(args.outdir, args.out_prefix),
        image_formats=args.formats if args.formats else ["png"],
        kind=args.kind,
        value_scope=args.value_scope,
        ncols=args.ncols,
        figsize_unit=(w, h),
        dpi=args.dpi,
        ylabel=args.ylabel,
        show_points=not args.no_points,
        point_jitter=args.point_jitter,
        point_size=args.point_size,
        point_alpha=args.point_alpha,
        seed=args.seed,
        test=args.test,
        test_ref_group=args.test_ref_group,
        x_rotation=args.x_rotation,
        box_width=args.box_width,
        violin_bw=args.violin_bw,
        row_spacing=args.row_spacing,
        stat_headroom=args.stat_headroom,
    )


if __name__ == "__main__":
    main()
