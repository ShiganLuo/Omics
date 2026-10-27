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
import types
from typing import Any, List, Optional, Tuple

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

PlotFormat = str

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
    """Load the TElocal locus count matrix.

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


def load_statistics_modules() -> Tuple[Any, Any]:
    """Import ``compare_groups`` and ``p_to_star`` from the statistics_plot submodule.

    The statistics_plot modules live in ``src/common/statistics_plot/Python``
    and import ``common.LogUtil``; that package does not exist inside the
    submodule, so the workflow ``LogUtil`` module is aliased under that name
    before importing.

    Returns
    -------
    Tuple[Any, Any]
        ``(compare_groups, p_to_star)`` callables from ``statistical_test``
        and ``bar_significance``.

    Raises
    ------
    SystemExit
        If the submodule or its scipy dependency is unavailable.
    """
    stats_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "common", "statistics_plot", "Python"
    )
    sys.path.insert(0, stats_dir)
    logutil = sys.modules.setdefault("LogUtil", types.ModuleType("LogUtil"))
    common_pkg = sys.modules.setdefault("common", types.ModuleType("common"))
    setattr(common_pkg, "LogUtil", logutil)
    sys.modules.setdefault("common.LogUtil", logutil)
    try:
        from bar_significance import p_to_star  # type: ignore[import-not-found]  # noqa: PLC0415
        from statistical_test import compare_groups  # type: ignore[import-not-found]  # noqa: PLC0415
    except ImportError as exc:
        raise SystemExit(
            f"statistics_plot modules unavailable ({exc}); --test requires "
            "src/common/statistics_plot (git submodule) and scipy."
        )
    return compare_groups, p_to_star


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
    bracket_lift: float = 1.15,
    bracket_gap: float = 1.18,
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
        ``bar_significance.p_to_star``. Default ``False`` (no test).
    test_ref_group : str, optional
        Reference group for the significance test (only used with
        ``do_test``; default: the first group, e.g. fibroblasts).
    test_alpha : float, optional
        Significance level for ``do_test``. Default ``0.05``.
    bracket_lift : float, optional
        Height of the first significance bracket as a fraction of the
        panel data range. Default ``1.15`` (above the tallest bar).
    bracket_gap : float, optional
        Vertical step between stacked brackets per group. Default ``1.18``.

    Returns
    -------
    None
    """
    if image_formats is None:
        image_formats = ["png"]
    groups = [str(g) for g in annot["group"].cat.categories]
    group_of = annot.set_index("sample")["group"].astype(str)

    patterns = list(subfamily_tpm.index)
    fig, ax = plt.subplots(figsize=figsize)
    n_g = len(groups)
    width = 0.8 / n_g
    stats = {}
    for gi, grp in enumerate(groups):
        samples = [s for s in subfamily_tpm.columns if group_of.get(s) == grp]
        means = subfamily_tpm[samples].mean(axis=1).reindex(patterns)
        sds = subfamily_tpm[samples].std(axis=1, ddof=1).reindex(patterns).fillna(0.0)
        xs = np.arange(len(patterns)) + (gi - n_g / 2 + 0.5) * width
        ax.bar(
            xs, means.values, width=width, yerr=sds.values, capsize=errbar_capsize,
            color=GROUP_PALETTE[gi % len(GROUP_PALETTE)], label=grp,
            error_kw=dict(lw=0.8, ecolor="#333333"),
        )
        stats[grp] = means

    # significance brackets (optional, do_test): each group vs the
    # reference group (default: first group, e.g. fibroblasts); the test is
    # auto-selected by statistical_test.compare_groups and the star label
    # comes from bar_significance.p_to_star (statistics_plot submodule)
    if do_test:
        compare_groups, p_to_star = load_statistics_modules()
        ref = test_ref_group if test_ref_group in groups else groups[0]
        ref_samples = [s for s in subfamily_tpm.columns if group_of.get(s) == ref]
        data_range = max(
            (float(subfamily_tpm.loc[p].max()) - float(subfamily_tpm.loc[p].min()))
            for p in patterns
        ) if patterns else 1.0
        floor = max(float(subfamily_tpm.min().min()), 1e-3)
        span = max(data_range, float(subfamily_tpm.max().max()) - floor)
        for pi, pattern in enumerate(patterns):
            level = bracket_lift
            for gi, grp in enumerate(groups):
                if grp == ref:
                    continue
                samples = [s for s in subfamily_tpm.columns if group_of.get(s) == grp]
                a = subfamily_tpm.loc[pattern, ref_samples].to_numpy(dtype=float)
                b = subfamily_tpm.loc[pattern, samples].to_numpy(dtype=float)
                if len(a) < 2 or len(b) < 2:
                    logger.warning(f"{pattern}: {ref} vs {grp} skipped (n<2 in a group)")
                    continue
                result = compare_groups(a.tolist(), b.tolist(), alpha=test_alpha, verbose=False)
                p_val = result["p_value"]
                logger.info(
                    f"{pattern} {ref} vs {grp}: {result['test_name']}, "
                    f"p={p_val:.4g} ({p_to_star(p_val)})"
                )
                y = floor + span * (level - 1.0) if yscale == "log" else span * level
                x_ref = pi + (groups.index(ref) - n_g / 2 + 0.5) * width
                x_grp = pi + (gi - n_g / 2 + 0.5) * width
                ax.plot([x_ref, x_ref, x_grp, x_grp],
                        [y - span * 0.012, y, y, y - span * 0.012],
                        lw=0.7, color="#333333", clip_on=False)
                ax.text((x_ref + x_grp) / 2, y, p_to_star(p_val),
                        ha="center", va="bottom", fontsize=8, color="#333333")
                level *= bracket_gap
            # expand y headroom so the highest bracket fits
            top_needed = floor + span * (level / bracket_gap - 1.0) if yscale == "log" else span * (level / bracket_gap)
            cur_top = ax.get_ylim()[1]
            if top_needed > cur_top:
                ax.set_ylim(bottom=ax.get_ylim()[0], top=top_needed * 1.05)
    ax.set_xticks(np.arange(len(patterns)))
    ax.set_xticklabels(patterns, rotation=x_rotation, ha="right")
    ax.set_ylabel(ylabel)
    if yscale == "log":
        ax.set_yscale("log")
    # legend outside the axes (right side) so it never covers bars
    ax.legend(frameon=False, fontsize=9, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
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
