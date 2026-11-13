#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quantitative bar plot of TE subfamily expression per sample group.

For each requested subfamily, aggregates the TEcount family-level
count matrix per sample group and draws grouped bars (mean ± sd
across samples). Two input flavors:

* TEcount matrix (default): rows are family/subfamily ids like
  ``MT2_Mm:MT2:SINE`` -- no lengths needed, values are used as-is
  (raw counts or DESeq2-normalized counts from a ``*_TE.tsv``).
* TElocal locus matrix (``--loci-tsv`` given): per-locus TPM is
  computed with locus lengths, then summed per subfamily.

Significance brackets use display-point (pt) offsets above the
tallest bar of each cluster, so spacing is identical on linear and
log axes; the y headroom is computed from the bracket stack before
rendering, so brackets never collide with error bars and are never
clipped.

Examples
--------
python te_subfamily_bar.py --matrix GRCh38_TEcount.tsv \
    --annotation samples.tsv --family HERVH-int --family SVA_D \
    -o outdir -f png -f pdf

Add ``--test`` to draw significance brackets (optional differential
analysis; test auto-selected via the statistics_plot submodule).
"""

import argparse
import os
import sys
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import logging
import matplotlib
import matplotlib.pyplot as plt
try:
    from .LogUtil import setup_logger
except ImportError:
    from LogUtil import setup_logger
logger = setup_logger(__name__, level=logging.INFO)
matplotlib.use("Agg")

# statistics_plot lives in the src/common submodule; make it importable
# (only used with --test).
_STATISTICS_PLOT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "src", "common", "statistics_plot", "Python",
)
if os.path.isdir(_STATISTICS_PLOT_DIR) and _STATISTICS_PLOT_DIR not in sys.path:
    sys.path.insert(0, _STATISTICS_PLOT_DIR)

PlotFormat = str

# matplotlib savefig's full format set (LogUtil.PlotFormat); anything else
# fails fast at parse time instead of dying inside savefig (modules.md #9)
VALID_FORMATS = ("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg",
                 "jpeg", "pgf", "raw", "rgba")

GROUP_PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#937860", "#64B5CD", "#CCB974",
]


def parse_te_id(te_id: str) -> dict:
    """Split a TE id into name/subfamily/family/class fields.

    Parameters
    ----------
    te_id : str
        3-segment (``name:family:class``) or 4-segment TElocal locus id
        (``name_dup:subfamily:family:class``).

    Returns
    -------
    dict
        Keys ``name``, ``subfamily``, ``family``, ``class``.
    """
    fields = str(te_id).split(":")
    if len(fields) >= 4:
        return {"name": fields[0], "subfamily": fields[1], "family": fields[2], "class": fields[3]}
    if len(fields) == 3:
        return {"name": fields[0], "subfamily": fields[0], "family": fields[1], "class": fields[2]}
    return {"name": str(te_id), "subfamily": str(te_id), "family": str(te_id), "class": str(te_id)}


def load_matrix(matrix_tsv: str) -> Tuple[pd.DataFrame, List[str]]:
    """Load the TElocus locus count matrix.

    Parameters
    ----------
    matrix_tsv : str
        Path to the count matrix TSV.

    Returns
    -------
    Tuple[pd.DataFrame, List[str]]
        ``(matrix, sample_columns)``.
    """
    raw = pd.read_csv(matrix_tsv, sep="\t")
    id_col = raw.columns[0]
    matrix = raw.set_index(id_col).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    logger.info(f"Loaded matrix {matrix.shape[0]} loci x {matrix.shape[1]} samples from {matrix_tsv}")
    return matrix, list(matrix.columns)


def load_annotation(annotation_tsv: str, sample_cols: List[str], group_order: Optional[str] = None) -> pd.DataFrame:
    """Load sample-to-group annotation with ordered categorical groups.

    Parameters
    ----------
    annotation_tsv : str
        TSV with ``sample``/``group`` (optional ``order``) columns.
    sample_cols : List[str]
        Matrix sample columns for validation.
    group_order : str, optional
        Comma-separated group order override.

    Returns
    -------
    pd.DataFrame
        Columns ``sample``, ``group`` (ordered categorical).
    """
    annot = pd.read_csv(annotation_tsv, sep="\t")
    annot = annot[annot["sample"].isin(sample_cols)].copy()
    if annot.empty:
        raise ValueError(f"No annotated samples overlap the matrix columns: {sample_cols}")
    if group_order:
        cats = [g.strip() for g in group_order.split(",") if g.strip()]
    elif "order" in annot.columns:
        cats = (annot.sort_values("order")["group"].astype(str).unique().tolist())
    else:
        cats = sorted(annot["group"].astype(str).unique())
    annot["group"] = pd.Categorical(annot["group"].astype(str), categories=cats, ordered=True)
    return annot.sort_values(["group", "sample"])


def load_lengths(loci_tsv: str, matrix_index: pd.Index) -> pd.Series:
    """Load per-locus lengths from a coordinate table.

    Parameters
    ----------
    loci_tsv : str
        ``.parquet`` cache, TElocal ``.locInd`` pickle, or TSV.
    matrix_index : pd.Index
        Matrix row ids to reindex to.

    Returns
    -------
    pd.Series
        Lengths indexed by locus id.
    """
    if loci_tsv.endswith(".parquet"):
        loci = pd.read_parquet(loci_tsv)
    elif loci_tsv.endswith(".locInd"):
        from te_locus_boxplot import load_loci_lengths  # noqa: PLC0415

        return load_loci_lengths(loci_tsv, matrix_index)
    else:
        loci = pd.read_csv(loci_tsv, sep="\t")
    lengths = (
        loci.assign(length=(loci["end"] - loci["start"]).clip(lower=1))
        .groupby("locus_id")["length"].max()
    )
    return lengths.reindex(matrix_index)


def compute_subfamily_value(
    matrix: pd.DataFrame,
    patterns: List[str],
    match_field: str,
    match_mode: str,
    value: str,
    lengths: Optional[pd.Series] = None,
    normalize: str = "cpm",
) -> pd.DataFrame:
    """Per-sample subfamily-level value (CPM or TPM sum).

    Parameters
    ----------
    matrix : pd.DataFrame
        Count matrix, rows = TE ids. With ``lengths`` given, treated as
        a TElocal locus matrix and converted to TPM; otherwise
        TEcount family counts normalized to CPM (``normalize="none"``
        uses raw values).
    lengths : pd.Series, optional
        Locus lengths indexed like ``matrix`` (TElocal mode).
    patterns : List[str]
        Subfamily patterns in plotting order.
    match_field : str
        TE id field to match (name/subfamily/family/class/any).
    match_mode : str
        substring/regex/exact matching mode.
    value : str
        ``"total"`` (sum over member rows, default) or ``"mean"``.
    normalize : str, optional
        ``"cpm"`` (default), ``"tpm"`` (requires ``lengths``) or
        ``"none"``.

    Returns
    -------
    pd.DataFrame
        Rows = patterns, columns = samples.

    Raises
    ------
    ValueError
        If ``normalize`` is invalid or TPM is requested without lengths.
    """
    if normalize not in ("cpm", "tpm", "none"):
        raise ValueError(f"normalize must be cpm/tpm/none, got {normalize!r}")
    if normalize == "tpm" and lengths is None:
        raise ValueError("normalize 'tpm' requires --loci-tsv")
    if lengths is not None:
        common = matrix.index.intersection(lengths.index)
        sub = matrix.loc[common]
        rpk = sub.div(lengths.loc[common].clip(lower=1), axis=0)
        data = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    elif normalize == "cpm":
        lib = matrix.sum(axis=0)
        data = matrix.div(lib.where(lib > 0), axis=1) * 1e6
    else:
        data = matrix

    fields = pd.DataFrame([parse_te_id(i) for i in data.index], index=data.index)
    out = {}
    for pattern in patterns:
        if match_field == "any":
            haystack = data.index.astype(str).str.lower()
        else:
            haystack = fields[match_field].astype(str).str.lower()
        needle = pattern.lower() if match_mode != "regex" else pattern
        if match_mode == "exact":
            mask = haystack == needle
        else:
            mask = haystack.str.contains(needle, regex=(match_mode == "regex"), na=False)
        member = data[np.asarray(mask, dtype=bool)]
        logger.info(f"Subfamily '{pattern}': {member.shape[0]} entries")
        out[pattern] = member.sum(axis=0) if value == "total" else member.mean(axis=0)
    return pd.DataFrame(out).T


def plot_subfamily_bars(
    subfamily_tpm: pd.DataFrame,
    annot: pd.DataFrame,
    out_base: str,
    image_formats: Optional[List[PlotFormat]] = None,
    figsize: Tuple[float, float] = (7.0, 4.2),
    dpi: int = 300,
    ylabel: str = "TPM",
    yscale: str = "linear",
    errbar_capsize: float = 3.0,
    x_rotation: float = 30.0,
    do_test: bool = False,
    test_ref_group: Optional[str] = None,
    test_alpha: float = 0.05,
) -> None:
    """Render grouped bars: one bar cluster per subfamily, one bar per group.

    Parameters
    ----------
    subfamily_tpm : pd.DataFrame
        Rows = subfamilies, columns = samples (TPM).
    annot : pd.DataFrame
        Sample/group annotation with ordered categorical ``group``.
    out_base : str
        Output path without extension.
    image_formats : List[str], optional
        Save formats. Default ``["png"]``.
    figsize : Tuple[float, float], optional
        Figure size in inches.
    dpi : int, optional
        Output resolution.
    ylabel : str, optional
        Y axis label.
    yscale : str, optional
        ``"linear"`` or ``"log"``.
    errbar_capsize : float, optional
        Error bar cap size in points.
    x_rotation : float, optional
        X tick label rotation.
    do_test : bool, optional
        Draw significance brackets (optional differential analysis). The
        statistical test per group pair is auto-selected by
        ``statistical_test.compare_groups`` (t-test / Welch / Mann-Whitney
        after normality and variance checks); stars follow
        ``bar_significance.p_to_star``. Brackets are laid out with
        display-point offsets above each cluster's tallest bar (incl.
        error bar), so they never collide with error bars and spacing is
        identical on linear and log axes. Default ``False`` (no test).
    test_ref_group : str, optional
        Reference group for the significance test (only used with
        ``do_test``; default: the first group, e.g. fibroblasts).
    test_alpha : float, optional
        Significance level for ``do_test``. Default ``0.05``.

    Returns
    -------
    None
    """
    if image_formats is None:
        image_formats = ["png"]
    bad = [f for f in image_formats if f not in VALID_FORMATS]
    if bad:
        raise ValueError(
            f"Unsupported image format(s): {', '.join(bad)}; "
            f"valid formats: {', '.join(VALID_FORMATS)}")
    groups = [str(g) for g in annot["group"].cat.categories]
    group_of = annot.set_index("sample")["group"].astype(str)

    # y scale first: subsequent display transforms must reflect it
    fig, ax = plt.subplots(figsize=figsize)
    if yscale == "log":
        ax.set_yscale("log")

    patterns = list(subfamily_tpm.index)
    n_g = len(groups)
    width = 0.8 / n_g
    stats = {}
    bar_tops = np.zeros(len(patterns), dtype=float)
    samples_by_group = {
        grp: [s for s in subfamily_tpm.columns if group_of.get(s) == grp] for grp in groups
    }
    for gi, grp in enumerate(groups):
        samples = samples_by_group[grp]
        means = subfamily_tpm[samples].mean(axis=1).reindex(patterns)
        sds = subfamily_tpm[samples].std(axis=1, ddof=1).reindex(patterns).fillna(0.0)
        xs = np.arange(len(patterns)) + (gi - n_g / 2 + 0.5) * width
        ax.bar(
            xs, means.values, width=width, yerr=sds.values, capsize=errbar_capsize,
            color=GROUP_PALETTE[gi % len(GROUP_PALETTE)], label=grp,
            error_kw=dict(lw=0.8, ecolor="#333333"),
        )
        stats[grp] = means
        bar_tops = np.maximum(bar_tops, (means + sds).to_numpy(dtype=float))

    ax.set_xticks(np.arange(len(patterns)))
    ax.set_xticklabels(patterns, rotation=x_rotation, ha="right")
    ax.set_ylabel(ylabel)
    # legend outside the axes (right side) so it never covers bars
    ax.legend(frameon=False, fontsize=9, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    # significance brackets (optional, do_test): each group vs the
    # reference group (default: first group, e.g. fibroblasts); test
    # auto-selected by statistical_test.compare_groups; bracket layout
    # and headroom come from bar_significance (statistics_plot
    # submodule) -- display-point based, so spacing is identical on
    # linear/log axes; headroom is computed after tight layout so
    # transforms are final.
    if do_test:
        sys.path.insert(
            0,
            os.path.join(os.path.dirname(os.path.abspath(__file__))),
        )
        # vendored copy in this directory (from statistics_plot submodule)
        from statistical_test import compare_groups  # noqa: PLC0415
        from bar_significance import draw_significance_bracket, fit_bracket_headroom, p_to_star  # type: ignore[import-not-found]  # noqa: PLC0415

        fig.canvas.draw()
        ref = test_ref_group if test_ref_group in groups else groups[0]
        x_of = {grp: np.arange(len(patterns)) + (gi - n_g / 2 + 0.5) * width
                for gi, grp in enumerate(groups)}

        # pre-pass: expand y headroom so every bracket stack fits
        fit_bracket_headroom(ax, bar_tops, max(1, n_g - 1))

        for pi, pattern in enumerate(patterns):
            base_data = max(float(bar_tops[pi]), 1e-9)
            level = 0
            for gi, grp in enumerate(groups):
                if grp == ref:
                    continue
                a = subfamily_tpm.loc[pattern, samples_by_group[ref]].to_numpy(dtype=float)
                b = subfamily_tpm.loc[pattern, samples_by_group[grp]].to_numpy(dtype=float)
                if len(a) < 2 or len(b) < 2:
                    logger.warning(f"{pattern}: {ref} vs {grp} skipped (n<2 in a group)")
                    continue
                result = compare_groups(a.tolist(), b.tolist(), alpha=test_alpha, verbose=False)
                p_val = result["p_value"]
                star = p_to_star(p_val)
                logger.info(
                    f"{pattern} {ref} vs {grp}: {result['test_name']}, "
                    f"p={p_val:.4g} ({star})"
                )
                draw_significance_bracket(
                    ax, x_of[ref][pi], x_of[grp][pi], base_data, star, level=level
                )
                level += 1

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
        description="Quantitative grouped-bar plot of TE subfamily TPM per sample group.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--matrix", required=True, help="TElocal locus count matrix TSV.")
    parser.add_argument("--annotation", required=True, help="sample/group annotation TSV.")
    parser.add_argument("--group-order", default=None, help="Comma-separated group order.")
    parser.add_argument("--family", action="append", dest="families", metavar="PATTERN", required=True,
                        help="Subfamily pattern. Repeatable; plotting order follows this list.")
    parser.add_argument("--match-field", choices=["name", "subfamily", "family", "class", "any"],
                        default="subfamily", help="TE id field to match the pattern against.")
    parser.add_argument("--match-mode", choices=["substring", "regex", "exact"], default="exact",
                        help="Pattern matching mode.")
    parser.add_argument("--value", choices=["total", "mean"], default="total",
                        help="Subfamily-level statistic: sum or mean over member rows.")
    parser.add_argument("--normalize", choices=["cpm", "tpm", "none"], default="cpm",
                        help="Normalization for TEcount mode: CPM (default), raw counts (none). "
                             "TPM applies only with --loci-tsv (TElocal mode).")
    parser.add_argument("--loci-tsv", default=None,
                        help="Optional locus coordinate table; when given the matrix is "
                             "treated as TElocal locus counts and converted to TPM. Without "
                             "it the matrix values are used as-is (TEcount).")
    parser.add_argument("--yscale", choices=["linear", "log"], default="linear",
                        help="Y axis scale. Use 'log' when subfamilies span orders of magnitude.")
    parser.add_argument("--test", dest="do_test", action="store_true",
                        help="Optional differential analysis: draw significance brackets of each "
                             "group vs the reference group. Test auto-selected by "
                             "statistics_plot.statistical_test.compare_groups (t-test/Welch/"
                             "Mann-Whitney after normality and variance checks).")
    parser.add_argument("--test-ref-group", default=None,
                        help="Reference group for --test brackets "
                             "(default: the first group in the annotation order).")
    parser.add_argument("--test-alpha", type=float, default=0.05,
                        help="Significance level for --test.")
    parser.add_argument("--ylabel", default="TPM", help="Y axis label.")
    parser.add_argument("--figsize", default="7.0,4.2", help="Figure size in inches, 'w,h'.")
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution.")
    parser.add_argument("-f", "--format", action="append", dest="formats", metavar="FMT",
                        help="Image format (png/pdf/svg/...). Repeatable. Default: png.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory.")
    parser.add_argument("--out-prefix", default="te_subfamily_bar", help="Output file prefix.")
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

    matrix, sample_cols = load_matrix(args.matrix)
    annot = load_annotation(args.annotation, sample_cols, group_order=args.group_order)
    lengths = load_lengths(args.loci_tsv, matrix.index) if args.loci_tsv else None
    subfamily_vals = compute_subfamily_value(
        matrix, args.families, args.match_field, args.match_mode, args.value,
        lengths=lengths, normalize=args.normalize,
    )
    if subfamily_vals.empty:
        raise ValueError("No entries matched any pattern")
    label = args.ylabel
    if label == "TPM" and args.value == "total":
        if lengths is not None:
            label = "sum TPM over loci"
        elif args.normalize == "cpm":
            label = "CPM"
        else:
            label = "counts"
    w, h = (float(v) for v in args.figsize.split(","))
    plot_subfamily_bars(
        subfamily_vals,
        annot,
        out_base=os.path.join(args.outdir, args.out_prefix),
        image_formats=args.formats if args.formats else ["png"],
        figsize=(w, h),
        dpi=args.dpi,
        ylabel=label,
        yscale=args.yscale,
        do_test=args.do_test,
        test_ref_group=args.test_ref_group,
        test_alpha=args.test_alpha,
    )


if __name__ == "__main__":
    main()
