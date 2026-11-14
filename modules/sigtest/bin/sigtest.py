#!/usr/bin/env python3
# Author: ShiganLuo
# Email: 25303020102@qq.com
# GitHub: https://github.com/ShiganLuo
"""Two-group significance testing with star labels and bracket annotation.

This module unifies the two legacy scripts ``statistical_test.py``
(assumption checks + automatic sample-based test selection) and
``bar_significance.py`` (count-based contingency tests, p-to-star labels,
display-point significance brackets) into one CLI-driven tool.

Data modes
----------
sample mode (methods: auto, t-test, welch, mannwhitney, wilcoxon)
    The long table holds one row per replicate: feature, group, value.
    Per feature the two groups' value vectors are compared. ``auto``
    checks normality (Shapiro-Wilk) and equal variance (Levene) and then
    picks paired t / Wilcoxon signed-rank (``--paired``) or Student t /
    Welch t / Mann-Whitney U. Effect sizes (Cohen's d or rank-biserial r)
    are reported.

count mode (methods: chi2, fisher)
    The long table holds one summarized count per feature x group. Per
    feature a 2x2 contingency table is built against the group totals
    (same construction as the legacy ``compute_significance``). chi2
    falls back to Fisher's exact test when expected frequencies hit zero.

Multiple-testing correction: none (default), bh (Benjamini-Hochberg) or
bonferroni. When a correction is active, stars and the ``significant``
flag are derived from the adjusted p-value.

Significance levels: p < 1e-4 "****", < 1e-3 "***", < 1e-2 "**",
< 0.05 "*", otherwise "ns".

Outputs
-------
<outdir>/significance_test.tsv
    Per-feature table: n, mean, std per group, test used, statistic,
    p_value, p_adjusted, stars, significant, effect size, note.
<outdir>/significance_bars.png
    Grouped bars (SEM error bars in sample mode) with significance
    brackets drawn in display points, so spacing is identical on any
    axis scale; the y headroom is auto-fitted to the tallest stack.

Usage:
    # continuous replicates, auto test selection, BH correction
    python sigtest.py -i expr_long.tsv -o sig_out/ \
        -f gene -g Control Mutant -m auto -c bh -y "Expression (TPM)"

    # summarized counts per feature (chi-square, Fisher fallback)
    python sigtest.py -i sv_counts.tsv -o sig_out/ \
        -f svtype -g Control Experiment -m chi2 --no-errorbars

    # paired design, auto picks paired t / Wilcoxon
    python sigtest.py -i paired.tsv -o sig_out/ --paired

Library use:
    from sigtest import test_long_table, plot_grouped_bars

    res = test_long_table(df, feature_col="gene", groups=("ctl", "mut"))
    plot_grouped_bars(df, res, "out.png", ylabel="TPM")
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, TextIO, Tuple

import numpy as np
import pandas as pd

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy import stats  # noqa: E402


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SAMPLE_METHODS: Tuple[str, ...] = ("auto", "t-test", "welch", "mannwhitney", "wilcoxon")
COUNT_METHODS: Tuple[str, ...] = ("chi2", "fisher")
ALL_METHODS: Tuple[str, ...] = SAMPLE_METHODS + COUNT_METHODS
CORRECTIONS: Tuple[str, ...] = ("none", "bh", "bonferroni")

# Statistical test + bracket layout come from the vendored house scripts
# (verbatim copies of src/common/statistics_plot/Python/{statistical_test,
# bar_significance}.py — modules.md 脚本自包含规范 #0): compare_groups does
# the testing, draw_significance_bracket/fit_bracket_headroom the layout.
from statistical_test import compare_groups  # noqa: E402
from bar_significance import (  # noqa: E402
    LEG_PT, TEXT_PT, FONT_ASCENT_PT, BRACKET_STEP_PT,
    y_to_px, px_to_y, fit_bracket_headroom, draw_significance_bracket,
)
from boxplot import plot_comparison  # noqa: E402

BAR_COLORS: Tuple[str, str] = ("#4C72B0", "#DD8452")

# matplotlib savefig's full format set (LogUtil.PlotFormat); anything else
# fails fast at parse time instead of dying inside savefig (modules.md #9)
VALID_FORMATS: Tuple[str, ...] = ("png", "pdf", "svg", "ps", "eps", "tif",
                                  "tiff", "jpg", "jpeg", "pgf", "raw", "rgba")
BAR_WIDTH = 0.36


# Copied verbatim from src/common/util/LogUtil.py (module bin scripts are
# self-contained — see modules.md 脚本自包含规范).
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


logger = setup_logger("SigTest", stream=sys.stdout)


# ---------------------------------------------------------------------------
# Stars and multiple-testing correction
# ---------------------------------------------------------------------------
def p_to_star(p: float) -> str:
    """Convert p-value to significance star label.

    Thresholds: p < 1e-4 "****", < 1e-3 "***", < 1e-2 "**", < 0.05 "*",
    otherwise "ns".
    """
    if p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def bh_adjust(pvalues: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg step-up adjustment with monotonicity enforcement.

    Parameters
    ----------
    pvalues : sequence of float
        Raw p-values (NaN entries are preserved as NaN).

    Returns
    -------
    np.ndarray
        Adjusted p-values, clipped to [0, 1], in the input order.
    """
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = int(ok.sum())
    if n == 0:
        return out
    q = p[ok]
    order = np.argsort(q, kind="stable")
    ranked = q[order] * n / (np.arange(n) + 1)
    # Enforce monotonicity from the largest p downwards.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(ranked, 0.0, 1.0)
    out[ok] = adjusted
    return out


def p_adjust(pvalues: Sequence[float], method: str = "bh") -> np.ndarray:
    """Adjust raw p-values for multiple testing.

    Parameters
    ----------
    pvalues : sequence of float
        Raw p-values (NaN entries are preserved as NaN).
    method : str, default="bh"
        'bh' (Benjamini-Hochberg FDR) or 'bonferroni'.

    Returns
    -------
    np.ndarray
        Adjusted p-values in the input order.

    Raises
    ------
    ValueError
        If the correction method is unsupported.
    """
    if method == "bh":
        return bh_adjust(pvalues)
    if method == "bonferroni":
        p = np.asarray(pvalues, dtype=float)
        n = max(int((~np.isnan(p)).sum()), 1)  # NaN rows do not inflate the family
        return np.where(np.isnan(p), np.nan, np.minimum(p * n, 1.0))
    raise ValueError(f"Unsupported correction: {method}")


# ---------------------------------------------------------------------------
# Per-feature testing over one long table (the unifying core)
# ---------------------------------------------------------------------------
def test_long_table(
    df: pd.DataFrame,
    feature_col: str = "feature",
    group_col: str = "group",
    value_col: str = "value",
    groups: Optional[Sequence[str]] = None,
    method: str = "auto",
    paired: bool = False,
    alpha: float = 0.05,
    correction: str = "none",
) -> pd.DataFrame:
    """Test every feature of a long table between two groups.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-format table. Sample mode: one row per replicate
        (feature, group, value). Count mode: exactly one row per
        feature x group carrying the summarized count.
    feature_col : str, default="feature"
        Column name identifying features (genes, SV types, ...).
    group_col : str, default="group"
        Column name identifying the two groups.
    value_col : str, default="value"
        Column name holding values or counts.
    groups : sequence of str, optional
        The two group labels. Defaults to the first two distinct labels
        in order of appearance; an error if the column holds a different
        number of distinct labels.
    method : str, default="auto"
        'auto', 't-test', 'welch', 'mannwhitney', 'wilcoxon' (sample
        mode) or 'chi2', 'fisher' (count mode).
    paired : bool, default=False
        Paired design for sample mode (ignored in count mode).
    alpha : float, default=0.05
        Significance level.
    correction : str, default="none"
        'none', 'bh' or 'bonferroni'. When active, stars and the
        ``significant`` flag use the adjusted p-value.

    Returns
    -------
    pandas.DataFrame
        One row per feature with columns: feature, group1, group2, n1,
        n2, mean1, mean2, std1, std2, test, statistic, p_value,
        p_adjusted, stars, significant, effect_size, effect_size_name,
        note. Features whose data cannot feed the requested test are
        kept with test='skipped' and the reason in ``note``.

    Raises
    ------
    ValueError
        On missing columns, ambiguous or invalid group selection,
        unsupported method/correction, or duplicate cells in count mode.
    """
    if method not in ALL_METHODS:
        raise ValueError(f"Unsupported method: {method} (choose from {', '.join(ALL_METHODS)})")
    if method in COUNT_METHODS and paired:
        raise ValueError(f"Method '{method}' tests summarized counts; --paired applies to sample methods only")
    if correction not in CORRECTIONS:
        raise ValueError(f"Unsupported correction: {correction} (choose from {', '.join(CORRECTIONS)})")
    for col in (feature_col, group_col, value_col):
        if col not in df.columns:
            raise ValueError(f"Column not found in input table: {col}")

    work = df.copy()
    work[value_col] = pd.to_numeric(work[value_col], errors="coerce")
    n_bad = int(work[value_col].isna().sum())
    if n_bad:
        logger.warning(f"Dropping {n_bad} non-numeric value rows")
    work = work.dropna(subset=[value_col])

    # Resolve the two groups.
    if groups is not None:
        if len(groups) != 2:
            raise ValueError(f"Exactly two group labels required, got: {list(groups)}")
        g1_name, g2_name = str(groups[0]), str(groups[1])
        missing = [g for g in (g1_name, g2_name) if g not in set(work[group_col])]
        if missing:
            raise ValueError(f"Group label(s) not found in column {group_col}: {', '.join(missing)}")
    else:
        distinct = list(dict.fromkeys(work[group_col].astype(str)))
        if len(distinct) != 2:
            raise ValueError(
                f"Column {group_col} holds {len(distinct)} distinct labels "
                f"({', '.join(distinct)}); pass exactly two via groups="
            )
        g1_name, g2_name = distinct

    work = work[work[group_col].isin([g1_name, g2_name])]
    features = list(dict.fromkeys(work[feature_col].astype(str)))

    rows: List[Dict[str, Any]] = []
    if method in COUNT_METHODS:
        dup = work.groupby([feature_col, group_col]).size()
        dup = dup[dup > 1]
        if len(dup):
            examples = [f"{f} / {g}" for f, g in list(dup.index)[:5]]
            raise ValueError(
                "Count mode needs exactly one value per feature x group; duplicates: "
                + ", ".join(examples)
                + (f" ({len(dup)} duplicated cells total)" if len(dup) > 5 else "")
            )
        pivot = work.pivot(index=feature_col, columns=group_col, values=value_col)
        pivot = pivot.reindex(features).fillna(0.0)
        total1 = float(pivot[g1_name].sum())
        total2 = float(pivot[g2_name].sum())

        for feature in features:
            c1 = float(pivot.loc[feature, g1_name])
            c2 = float(pivot.loc[feature, g2_name])
            table = np.array([
                [c1, total1 - c1],
                [c2, total2 - c2],
            ])
            note = ""
            if method == "fisher":
                test_name = "Fisher's exact test"
                statistic, p_value = stats.fisher_exact(table)
            else:
                try:
                    statistic, p_value, _, _ = stats.chi2_contingency(table)
                    test_name = "Chi-square test"
                except ValueError:
                    # Zero expected frequencies: fall back to Fisher.
                    statistic, p_value = stats.fisher_exact(table)
                    test_name = "Fisher's exact test (chi2 fallback)"
                    note = "expected frequencies hit zero"
            rows.append({
                "feature": feature,
                "group1": g1_name, "group2": g2_name,
                "n1": total1, "n2": total2,
                "mean1": c1, "mean2": c2,          # in count mode the "means" are the counts
                "std1": np.nan, "std2": np.nan,
                "test": test_name,
                "statistic": float(statistic),
                "p_value": float(p_value),
                "effect_size": np.nan, "effect_size_name": "",
                "note": note,
            })
    else:
        for feature in features:
            sub = work[work[feature_col] == feature]
            g1_vals = sub.loc[sub[group_col] == g1_name, value_col].to_numpy(dtype=float)
            g2_vals = sub.loc[sub[group_col] == g2_name, value_col].to_numpy(dtype=float)
            try:
                # vendored statistical_test.compare_groups: auto test
                # selection (Shapiro/Levene -> t/welch/MW/wilcoxon); a
                # forced sample method goes through the same scipy calls
                res = compare_groups(
                    [float(x) for x in g1_vals], [float(x) for x in g2_vals],
                    paired=paired, alpha=alpha, verbose=False,
                )
                if method != "auto":
                    from scipy import stats as _st
                    if method == "wilcoxon":
                        stat, p_value = _st.wilcoxon(g1_vals, g2_vals)
                        res["test_name"] = "Wilcoxon signed-rank test"
                    elif method == "t-test":
                        stat, p_value = _st.ttest_ind(g1_vals, g2_vals, equal_var=True)
                        res["test_name"] = "Independent t-test (equal variance)"
                    elif method == "welch":
                        stat, p_value = _st.ttest_ind(g1_vals, g2_vals, equal_var=False)
                        res["test_name"] = "Welch's t-test"
                    else:  # mannwhitney
                        stat, p_value = _st.mannwhitneyu(g1_vals, g2_vals,
                                                         alternative="two-sided")
                        res["test_name"] = "Mann-Whitney U test"
                    res["statistic"] = float(stat)
                    res["p_value"] = float(p_value)
                    res["significant"] = bool(p_value <= alpha)
            except ValueError as exc:
                rows.append({
                    "feature": feature,
                    "group1": g1_name, "group2": g2_name,
                    "n1": len(g1_vals), "n2": len(g2_vals),
                    "mean1": float(np.mean(g1_vals)) if len(g1_vals) else np.nan,
                    "mean2": float(np.mean(g2_vals)) if len(g2_vals) else np.nan,
                    "std1": float(np.std(g1_vals, ddof=1)) if len(g1_vals) > 1 else np.nan,
                    "std2": float(np.std(g2_vals, ddof=1)) if len(g2_vals) > 1 else np.nan,
                    "test": "skipped",
                    "statistic": np.nan,
                    "p_value": np.nan,
                    "effect_size": np.nan, "effect_size_name": "",
                    "note": str(exc),
                })
                continue
            rows.append({
                "feature": feature,
                "group1": g1_name, "group2": g2_name,
                "n1": len(g1_vals), "n2": len(g2_vals),
                "mean1": res["group1_mean"], "mean2": res["group2_mean"],
                "std1": res["group1_std"], "std2": res["group2_std"],
                "test": res["test_name"],
                "statistic": res["statistic"],
                "p_value": res["p_value"],
                "effect_size": res["effect_size"],
                "effect_size_name": res["effect_size_name"],
                "note": "",
            })

    result = pd.DataFrame(rows)
    result["p_adjusted"] = (
        p_adjust(result["p_value"].to_numpy(), correction)
        if correction != "none" else result["p_value"]
    )
    p_used = result["p_adjusted"].to_numpy(dtype=float)
    result["stars"] = [p_to_star(p) if not np.isnan(p) else "" for p in p_used]
    result["significant"] = [(not np.isnan(p)) and p <= alpha for p in p_used]
    return result[
        ["feature", "group1", "group2", "n1", "n2", "mean1", "mean2", "std1", "std2",
         "test", "statistic", "p_value", "p_adjusted", "stars", "significant",
         "effect_size", "effect_size_name", "note"]
    ]


# ---------------------------------------------------------------------------
# Grouped bar plot with significance brackets
# ---------------------------------------------------------------------------
def _pivot_group_stats(
    sub: pd.DataFrame,
    feature_col: str,
    group_col: str,
    value_col: str,
    features: Sequence[str],
    g1: str,
    g2: str,
    aggfunc: str,
) -> pd.DataFrame:
    """Pivot to per-feature x per-group stats with both group columns guaranteed.

    ``pivot_table`` drops columns whose cells are all NaN (e.g. ``sem``
    on one-row-per-cell count tables), so the result is reindexed to the
    two group columns explicitly.
    """
    return (
        sub.pivot_table(index=feature_col, columns=group_col, values=value_col,
                        aggfunc=aggfunc)
        .reindex(index=list(features), columns=[g1, g2])
        .fillna(0.0)
        .astype(float)
    )


def plot_feature_distributions(
    df: pd.DataFrame,
    results: pd.DataFrame,
    out_dir: str,
    feature_col: str = "feature",
    group_col: str = "group",
    value_col: str = "value",
    groups: Optional[Sequence[str]] = None,
    features: Optional[Sequence[str]] = None,
    ylabel: str = "value",
    formats: Sequence[str] = ("png",),
    dpi: int = 300,
    processes: int = 1,
    alpha: float = 0.05,
) -> List[str]:
    """Draw ONE DISTRIBUTION FIGURE PER FEATURE (house boxplot style).

    Each figure is the house style from statistics_plot (boxplot.py::
    plot_comparison over statistical_test.compare_groups): box + jittered
    points showing the two groups' DISTRIBUTIONS, significance stars and
    p-value on a bracket, effect size (Cohen's d / rank-biserial with
    interpretation) top-right, test name bottom-left. The old mean-bar
    rendering dropped the distribution entirely — with thousands of loci
    per group the distribution IS the story.

    A single all-features figure exceeds matplotlib's 2^16-px canvas past
    ~380 features; per-feature figures are bounded. Figures are
    independent work units — ``processes`` > 1 shards them across a pool
    (a 1271-feature set renders in ~1 min at 8 processes).

    Parameters
    ----------
    df, results, feature_col, group_col, value_col, groups, features, ylabel
        Same meaning as in test_long_table / plot_grouped_bars.
    formats : sequence of str
        Output image formats (matplotlib savefig set); one file each.
    dpi : int
        Output resolution.
    processes : int, default 1
        Worker processes for figure rendering (1 = serial).
    alpha : float
        Significance level for the per-feature comparison.

    Returns
    -------
    list of str
        The written file paths.
    """
    work = df.dropna(subset=[value_col]).copy()
    if groups is None:
        distinct = list(dict.fromkeys(work[group_col].astype(str)))
        if len(distinct) < 2:
            raise ValueError("Need at least two groups to plot")
        groups = distinct[:2]
    g1_name, g2_name = str(groups[0]), str(groups[1])

    if features is None:
        features = results["feature"].astype(str).tolist()
    features = [str(f) for f in features]

    sub = work[work[group_col].isin([g1_name, g2_name])]
    by_feat_g1: Dict[str, List[float]] = {}
    by_feat_g2: Dict[str, List[float]] = {}
    g1_rows = sub[sub[group_col] == g1_name]
    g2_rows = sub[sub[group_col] == g2_name]
    for feat, vals in g1_rows.groupby(feature_col)[value_col]:
        by_feat_g1[str(feat)] = [float(v) for v in vals]
    for feat, vals in g2_rows.groupby(feature_col)[value_col]:
        by_feat_g2[str(feat)] = [float(v) for v in vals]

    # per-feature stats come from the RESULTS table (already computed by
    # test_long_table) — the renderer never re-tests; skipped features
    # (n<2 etc.) have no figure, same as they have no tested row
    res_rows = results.set_index(results["feature"].astype(str))
    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    jobs = []
    for feat in features:
        if feat not in by_feat_g1 or feat not in by_feat_g2:
            continue
        if feat not in res_rows.index or str(res_rows.loc[feat, "test"]) == "skipped":
            continue
        row = res_rows.loc[feat]
        jobs.append(
            (feat, str(out_root), g1_name, g2_name, ylabel,
             by_feat_g1[feat], by_feat_g2[feat],
             tuple(formats), dpi,
             {
                 "test_name": str(row["test"]),
                 "p_value": float(row["p_value"]),
                 "significant": str(row["significant"]).lower() in ("true", "1"),
                 "effect_size": float(row["effect_size"]),
                 "effect_size_name": str(row["effect_size_name"]) or "effect size",
                 # same thresholds as statistical_test.compare_groups
                 "effect_size_interpretation": (
                     "negligible" if abs(float(row["effect_size"])) < 0.2
                     else "small" if abs(float(row["effect_size"])) < 0.5
                     else "medium" if abs(float(row["effect_size"])) < 0.8
                     else "large"),
             })
        )
    if processes > 1 and len(jobs) > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(min(processes, len(jobs))) as pool:
            written_nested = pool.map(_render_feature_distribution, jobs)
        written = [p for sub_list in written_nested for p in sub_list]
    else:
        written = []
        for job in jobs:
            written.extend(_render_feature_distribution(job))
    return written


def _render_feature_distribution(job: tuple) -> List[str]:
    """Render one feature's distribution comparison; module-level for pickling.

    Stats arrive in the job (from the results table — never re-tested);
    boxplot.plot_comparison (vendored house style) does the rendering.
    """
    (feat, out_root, g1_name, g2_name, ylabel, g1_vals, g2_vals,
     formats, dpi, comparison) = job
    result = {
        "group1_name": g1_name, "group2_name": g2_name,
        "group1_data": g1_vals, "group2_data": g2_vals,
        "comparison": comparison,
    }
    safe = re.sub(r"[^\w.\-+=(),]", "_", feat)
    written: List[str] = []
    for fmt in formats:
        out_path = os.path.join(out_root, f"{safe}.{fmt}")
        fig = plot_comparison(
            result, title=feat, xlabel="", ylabel=ylabel,
            show_points=True, show_stats=True,
            save_path=out_path, dpi=dpi, verbose=False,
        )
        plt.close(fig)
        written.append(out_path)
    return written


def plot_grouped_bars(
    df: pd.DataFrame,
    results: pd.DataFrame,
    out_png: str,
    feature_col: str = "feature",
    group_col: str = "group",
    value_col: str = "value",
    groups: Optional[Sequence[str]] = None,
    feature_order: Optional[Sequence[str]] = None,
    ylabel: str = "value",
    legend_labels: Optional[Sequence[str]] = None,
    colors: Tuple[str, str] = BAR_COLORS,
    annotate_ns: bool = False,
    show_errorbars: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 300,
) -> None:
    """Draw grouped bars per feature with significance brackets on top.

    Uses the per-feature results from :func:`test_long_table` for the
    star labels. Bars show group means with SEM error bars (sample
    mode); the bracket stack headroom is fitted automatically.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-format input table (same shape as given to test_long_table).
    results : pandas.DataFrame
        Output of test_long_table for ``df``.
    out_png : str
        Output PNG path (parent directories are created).
    feature_col, group_col, value_col : str
        Column names in ``df``.
    groups : sequence of str, optional
        Two group labels; defaults to the first two distinct labels.
    feature_order : sequence of str, optional
        X-axis feature order; defaults to appearance order.
    ylabel : str, default="value"
        Y-axis label.
    legend_labels : sequence of str, optional
        Legend labels for the two groups; defaults to group names.
    colors : tuple of str
        Bar colors for the two groups.
    annotate_ns : bool, default=False
        Also draw brackets for non-significant features ("ns").
    show_errorbars : bool, default=True
        Draw SEM error bars (single-value cells get no error bar).
    figsize : tuple of float, optional
        Figure size; defaults to width max(9, 2.5 + 0.55 * n_features).
    dpi : int, default=300
        Output resolution.

    Returns
    -------
    None
        The figure is saved to ``out_png``.
    """
    work = df.dropna(subset=[value_col]).copy()
    if groups is None:
        distinct = list(dict.fromkeys(work[group_col].astype(str)))
        if len(distinct) < 2:
            raise ValueError("Need at least two groups to plot")
        groups = distinct[:2]
    g1_name, g2_name = str(groups[0]), str(groups[1])
    if legend_labels is None:
        legend_labels = (g1_name, g2_name)

    if feature_order is None:
        feature_order = results["feature"].tolist()
    features = [str(f) for f in feature_order]
    n = len(features)
    if figsize is None:
        figsize = (max(9.0, 2.5 + 0.55 * n), 5.0)

    sub = work[work[group_col].isin([g1_name, g2_name])]
    means = _pivot_group_stats(sub, feature_col, group_col, value_col, features,
                               g1_name, g2_name, "mean")
    sems = _pivot_group_stats(sub, feature_col, group_col, value_col, features,
                              g1_name, g2_name, "sem")
    counts = _pivot_group_stats(sub, feature_col, group_col, value_col, features,
                                g1_name, g2_name, "count")

    stars = dict(zip(results["feature"].astype(str), results["stars"]))

    x = np.arange(n)
    fig, ax = plt.subplots(figsize=figsize)

    for offset, name, color, label in (
        (-BAR_WIDTH / 2, g1_name, colors[0], legend_labels[0]),
        (+BAR_WIDTH / 2, g2_name, colors[1], legend_labels[1]),
    ):
        vals = means[name].to_numpy(dtype=float)
        err = sems[name].to_numpy(dtype=float)
        reps = counts[name].to_numpy(dtype=float)
        if not show_errorbars:
            err = np.zeros_like(err)
        else:
            err = np.where(reps > 1, err, 0.0)  # single value: no error bar
        ax.bar(x + offset, vals, BAR_WIDTH, color=color, label=str(label),
               yerr=err, capsize=2, error_kw={"lw": 0.8})

    # Formatting first: brackets need finalized transforms.
    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.canvas.draw()

    # Bracket pass: base = tallest bar top including its error bar.
    annotated = [
        (i, feat) for i, feat in enumerate(features)
        if stars.get(feat, "") and (annotate_ns or stars.get(feat, "") != "ns")
    ]
    bases = []
    for i, _feat in annotated:
        top1 = means.iloc[i][g1_name] + (sems.iloc[i][g1_name] if show_errorbars else 0.0)
        top2 = means.iloc[i][g2_name] + (sems.iloc[i][g2_name] if show_errorbars else 0.0)
        bases.append(float(max(top1, top2)))
    fit_bracket_headroom(ax, bases, n_brackets=1 if bases else 0)
    for (i, feat), base in zip(annotated, bases):
        draw_significance_bracket(
            ax, x[i] - BAR_WIDTH / 2, x[i] + BAR_WIDTH / 2, base,
            stars[feat], level=0,
        )

    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # savefig infers the format from the extension (multi-format support)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    logger.info(f"Plot saved: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def read_input_table(path: Path, sep: Optional[str] = None) -> pd.DataFrame:
    """Read a TSV/CSV table, auto-detecting tab vs comma when sep is None.

    Parameters
    ----------
    path : pathlib.Path
        Input table path.
    sep : str, optional
        Explicit delimiter; default sniffs the header line.

    Returns
    -------
    pandas.DataFrame
        The parsed table.

    Raises
    ------
    ValueError
        If the file is empty.
    """
    if sep is None:
        with open(path, "r", encoding="utf-8") as fh:
            header = fh.readline()
        sep = "\t" if header.count("\t") >= header.count(",") else ","
    df = pd.read_csv(path, sep=sep)
    if df.empty:
        raise ValueError(f"Input table is empty: {path}")
    return df


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments for sigtest.py."""
    p = argparse.ArgumentParser(
        prog="sigtest.py",
        description="Two-group per-feature significance testing with star "
                    "labels, multiple-testing correction and a bracket "
                    "bar plot.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", type=Path, required=True,
                   help="long-format input table (TSV or CSV)")
    p.add_argument("-o", "--outdir", type=Path, required=True,
                   help="output directory (results TSV and PNG)")
    p.add_argument("-f", "--feature-col", default="feature",
                   help="column naming the features (genes, SV types, ...)")
    p.add_argument("--group-col", default="group",
                   help="column naming the groups")
    p.add_argument("-v", "--value-col", default="value",
                   help="column holding values (replicates) or counts")
    p.add_argument("-g", "--groups", nargs=2, default=None, metavar=("G1", "G2"),
                   help="the two group labels; default: first two distinct labels")
    p.add_argument("-m", "--method", default="auto", choices=ALL_METHODS,
                   help="statistical test; sample-based (auto/t-test/welch/"
                        "mannwhitney/wilcoxon) or count-based (chi2/fisher)")
    p.add_argument("--paired", action="store_true",
                   help="paired design (sample methods only)")
    p.add_argument("-a", "--alpha", type=float, default=0.05,
                   help="significance level")
    p.add_argument("-c", "--correction", default="none", choices=CORRECTIONS,
                   help="multiple-testing correction across features")
    p.add_argument("-s", "--sep", default=None,
                   help="input delimiter; default auto-detects tab vs comma")
    p.add_argument("--formats", action="append", dest="formats", metavar="FMT",
                   default=None,
                   help="image format(s) for the bar figures (png/pdf/svg/...; "
                        "matplotlib savefig set). Repeatable, e.g. --formats png "
                        "--formats pdf. An explicit empty value (--formats '') "
                        "means no figure; default png when absent.")
    p.add_argument("--overview-threshold", type=int, default=60, metavar="N",
                   help="emit the all-features overview figure only when the "
                        "feature count is <= N (its canvas is unbounded; "
                        "per-feature figures are always written). 0 disables "
                        "the overview. Default 60.")
    p.add_argument("-t", "--threads", type=int, default=1, metavar="N",
                   help="worker processes for per-feature figure rendering "
                        "(1 = serial).")
    p.add_argument("--annotate-ns", action="store_true",
                   help="draw brackets also for non-significant features")
    p.add_argument("--no-errorbars", action="store_true",
                   help="hide SEM error bars on the plot")
    p.add_argument("-y", "--ylabel", default="value",
                   help="y-axis label on the plot")
    p.add_argument("--legend", nargs=2, default=None, metavar=("L1", "L2"),
                   help="legend labels for the two groups")
    p.add_argument("-d", "--dpi", type=int, default=300,
                   help="output PNG resolution")
    p.add_argument("-l", "--log", type=Path, default=None,
                   help="log file path (default: log to stdout only)")
    p.add_argument("-n", "--dry-run", action="store_true",
                   help="log the plan and exit without testing")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the per-feature significance workflow from the CLI."""
    global logger
    args = parse_args(argv)
    if args.log is not None:
        args.log.parent.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("SigTest", log_file=str(args.log) if args.log else None,
                          stream=sys.stdout)
    logger.info(f"Command: {' '.join(sys.argv)}")

    if not args.input.is_file():
        logger.error(f"Input table not found: {args.input}")
        return 2

    try:
        df = read_input_table(args.input, args.sep)
    except ValueError as exc:
        logger.error(str(exc))
        return 2

    feature_count = df[args.feature_col].nunique() if args.feature_col in df.columns else -1
    group_labels = (
        list(dict.fromkeys(df[args.group_col].astype(str)))
        if args.group_col in df.columns else []
    )
    logger.info(
        f"Input: {args.input} rows={len(df)} features={feature_count} "
        f"groups={group_labels}"
    )

    plan = (
        f"Plan: method={args.method} paired={args.paired} alpha={args.alpha} "
        f"correction={args.correction} groups={args.groups or 'auto'} "
        f"formats={args.formats if args.formats is not None else ['png']}"
    )
    if args.dry_run:
        logger.info(f"[DRY-RUN] {plan}")
        return 0
    logger.info(plan)

    try:
        results = test_long_table(
            df,
            feature_col=args.feature_col,
            group_col=args.group_col,
            value_col=args.value_col,
            groups=args.groups,
            method=args.method,
            paired=args.paired,
            alpha=args.alpha,
            correction=args.correction,
        )
    except ValueError as exc:
        logger.error(str(exc))
        return 2

    args.outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = args.outdir / "significance_test.tsv"
    results.to_csv(tsv_path, sep="\t", index=False)
    logger.info(f"Results table saved: {tsv_path}")

    skipped = results[results["test"] == "skipped"]
    if len(skipped):
        logger.warning(
            f"{len(skipped)} feature(s) skipped: "
            + ", ".join(f"{r.feature} ({r.note})" for r in skipped.itertuples())
        )
        if len(skipped) == len(results):
            logger.error("All features were skipped; nothing tested")
            return 1

    tested = results[results["test"] != "skipped"]
    n_sig = int(tested["significant"].sum())
    logger.info(
        f"Tested {len(tested)} features: {n_sig} significant at "
        f"alpha={args.alpha} ({args.correction} correction)"
    )
    for test_name, count in tested["test"].value_counts().items():
        logger.info(f"  {test_name}: {count}")

    # formats: None (flag absent) = default ["png"]; an explicit empty
    # value (--formats '') = no figure. No separate boolean switch.
    formats = ([f for f in (str(x) for x in args.formats) if f]
               if args.formats is not None else ["png"])
    bad = [f for f in formats if f not in VALID_FORMATS]
    if bad:
        raise ValueError(
            f"Unsupported image format(s): {', '.join(bad)}; "
            f"valid formats: {', '.join(VALID_FORMATS)}")
    if formats:
        # ONE FIGURE PER FEATURE (default; bounded canvas at any feature
        # count — a single all-in figure exceeds matplotlib's 2^16-px
        # limit once features pass ~380 at 300 dpi)
        written = plot_feature_distributions(
            df, results, str(args.outdir / "significance_bars"),
            feature_col=args.feature_col, group_col=args.group_col,
            value_col=args.value_col, groups=args.groups,
            ylabel=args.ylabel,
            formats=formats, dpi=args.dpi,
            processes=max(1, args.threads),
            alpha=args.alpha,
        )
        logger.info(f"Per-feature figures: {len(written)} file(s) under "
                    f"{args.outdir / 'significance_bars'}")
        # DAG marker for the figure subdir (modules.md: no directory()
        # outputs; the rule declares significance_bars.done)
        (args.outdir / "significance_bars.done").write_text(
            "\n".join(written) + "\n")
        # overview (all features in one) only when the canvas stays small;
        # --overview-threshold <= 0 disables it entirely
        n_feat = len(results)
        if 0 < args.overview_threshold and n_feat <= args.overview_threshold:
            for fmt in formats:
                out_png = str(args.outdir / f"significance_bars_overview.{fmt}")
                plot_grouped_bars(
                    df, results, out_png,
                    feature_col=args.feature_col, group_col=args.group_col,
                    value_col=args.value_col, groups=args.groups,
                    feature_order=results["feature"].tolist(),
                    ylabel=args.ylabel, legend_labels=args.legend,
                    annotate_ns=args.annotate_ns,
                    show_errorbars=not args.no_errorbars,
                    dpi=args.dpi,
                )
                logger.info(f"Overview plot saved: {out_png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
