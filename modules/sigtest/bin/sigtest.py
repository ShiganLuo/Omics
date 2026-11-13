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

# Bracket geometry in display points (identical on linear and log axes):
# the first bracket sits LEG_PT above the cluster's tallest bar top; each
# further stacked bracket adds BRACKET_STEP_PT; the star label sits TEXT_PT
# above its bracket line; FONT_ASCENT_PT approximates the label font's
# ascent for headroom reservation.
LEG_PT = 10.0
TEXT_PT = 4.0
FONT_ASCENT_PT = 8.0
BRACKET_STEP_PT = 14.0

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
# Assumption checks (from statistical_test.py)
# ---------------------------------------------------------------------------
def check_normality(data: List[float], alpha: float = 0.05) -> Tuple[bool, float, str]:
    """Check if data follows a normal distribution using Shapiro-Wilk test.

    Parameters
    ----------
    data : list of float
        The data to test for normality.
    alpha : float, optional
        Significance level for the test (default: 0.05).

    Returns
    -------
    is_normal : bool
        True if data appears to be normally distributed (p > alpha).
    p_value : float
        The p-value from the Shapiro-Wilk test.
    message : str
        Interpretation of the test result.

    Notes
    -----
    - H0: The data is normally distributed.
    - If p > alpha, we fail to reject H0 (data appears normal).
    - Requires at least 3 data points; fewer are treated as non-normal.
    """
    if len(data) < 3:
        return False, 0.0, "Sample size too small for normality test (n < 3)"

    stat, p_value = stats.shapiro(data)
    is_normal = p_value > alpha

    if is_normal:
        message = f"Data appears normally distributed (Shapiro-Wilk p={p_value:.4f} > {alpha})"
    else:
        message = f"Data does not appear normally distributed (Shapiro-Wilk p={p_value:.4f} <= {alpha})"

    return is_normal, p_value, message


def check_equal_variance(group1: List[float], group2: List[float], alpha: float = 0.05) -> Tuple[bool, float, str]:
    """Check if two groups have equal variances using Levene's test.

    Parameters
    ----------
    group1 : list of float
        First group of data.
    group2 : list of float
        Second group of data.
    alpha : float, optional
        Significance level for the test (default: 0.05).

    Returns
    -------
    equal_variance : bool
        True if groups appear to have equal variances (p > alpha).
    p_value : float
        The p-value from Levene's test.
    message : str
        Interpretation of the test result.
    """
    stat, p_value = stats.levene(group1, group2)
    equal_variance = p_value > alpha

    if equal_variance:
        message = f"Groups appear to have equal variances (Levene p={p_value:.4f} > {alpha})"
    else:
        message = f"Groups appear to have unequal variances (Levene p={p_value:.4f} <= {alpha})"

    return equal_variance, p_value, message


# ---------------------------------------------------------------------------
# Sample-based comparison (from statistical_test.py, with method override)
# ---------------------------------------------------------------------------
def compare_groups(
    group1: List[float],
    group2: List[float],
    paired: bool = False,
    alpha: float = 0.05,
    verbose: bool = True,
    method: str = "auto",
) -> Dict[str, Any]:
    """Compare two groups to determine if there is a statistically significant difference.

    ``method='auto'`` selects the test from the data characteristics
    (normality, equal variance); the other choices force one test.

    Parameters
    ----------
    group1 : list of float
        First group of numerical data.
    group2 : list of float
        Second group of numerical data.
    paired : bool, optional
        If True, perform a paired test (default: False).
    alpha : float, optional
        Significance level for hypothesis tests (default: 0.05).
    verbose : bool, optional
        If True, log detailed information about the tests (default: True).
    method : str, optional
        'auto' (default), 't-test', 'welch', 'mannwhitney' or 'wilcoxon'.

    Returns
    -------
    result : dict
        Dictionary containing:
        - 'test_name': Name of the statistical test used
        - 'statistic': Test statistic value
        - 'p_value': p-value from the test
        - 'significant': Boolean indicating if difference is significant (p <= alpha)
        - 'alpha': Significance level used
        - 'group1_mean', 'group2_mean', 'group1_std', 'group2_std'
        - 'effect_size': Cohen's d (parametric) or rank-biserial r
        - 'effect_size_name', 'effect_size_interpretation'
        - 'interpretation': Human-readable interpretation of results
        - 'assumptions': Dictionary with assumption check results

    Raises
    ------
    ValueError
        If either group has fewer than 2 data points, pairing is requested
        for an independent-only test (or vice versa), or paired groups have
        unequal lengths.
    """
    if method not in SAMPLE_METHODS:
        raise ValueError(f"Unsupported sample method: {method}")
    if method == "mannwhitney" and paired:
        raise ValueError("Mann-Whitney U tests independent samples; drop --paired or use 'wilcoxon'")
    if method == "welch" and paired:
        raise ValueError("Welch's t-test is for independent samples; use paired t-test instead")
    if len(group1) < 2 or len(group2) < 2:
        raise ValueError("Each group must have at least 2 data points")

    use_paired = paired or method == "wilcoxon"
    if use_paired and len(group1) != len(group2):
        raise ValueError("Paired samples must have the same number of observations")

    # Convert to numpy arrays
    g1 = np.array(group1, dtype=float)
    g2 = np.array(group2, dtype=float)

    # Calculate basic statistics
    g1_mean, g1_std = np.mean(g1), np.std(g1, ddof=1)
    g2_mean, g2_std = np.mean(g2), np.std(g2, ddof=1)

    # Check assumptions
    assumptions: Dict[str, Any] = {}

    # Normality tests
    g1_normal, g1_normal_p, g1_normal_msg = check_normality(g1, alpha)
    g2_normal, g2_normal_p, g2_normal_msg = check_normality(g2, alpha)

    assumptions['group1_normality'] = {
        'is_normal': g1_normal,
        'p_value': g1_normal_p,
        'message': g1_normal_msg
    }
    assumptions['group2_normality'] = {
        'is_normal': g2_normal,
        'p_value': g2_normal_p,
        'message': g2_normal_msg
    }

    # Both groups normal?
    both_normal = g1_normal and g2_normal

    # Equal variance test (only if both groups are normal)
    if both_normal:
        equal_var, equal_var_p, equal_var_msg = check_equal_variance(g1, g2, alpha)
        assumptions['equal_variance'] = {
            'equal': equal_var,
            'p_value': equal_var_p,
            'message': equal_var_msg
        }
    else:
        equal_var = False
        assumptions['equal_variance'] = {
            'equal': False,
            'p_value': None,
            'message': 'Skipped (not both groups are normal)'
        }

    # Select and perform the test. Forced methods decide purely from the
    # method name; 'auto' decides from the assumption checks.
    need_paired = use_paired
    need_normal = both_normal if method == "auto" else method == "t-test"

    if need_paired:
        if need_normal:
            # Paired t-test
            test_name = "Paired t-test"
            statistic, p_value = stats.ttest_rel(g1, g2)
            # Effect size: Cohen's d for paired samples
            diff = g1 - g2
            effect_size = np.mean(diff) / np.std(diff, ddof=1)
            effect_size_name = "Cohen's d"
        else:
            # Wilcoxon signed-rank test
            test_name = "Wilcoxon signed-rank test"
            statistic, p_value = stats.wilcoxon(g1, g2)
            # Effect size: rank-biserial correlation
            n = len(g1)
            effect_size = statistic / (n * (n + 1) / 2)
            effect_size_name = "rank-biserial r"
    elif need_normal:
        if equal_var:
            # Independent t-test (equal variance)
            test_name = "Independent t-test (equal variance)"
            statistic, p_value = stats.ttest_ind(g1, g2, equal_var=True)
        else:
            # Welch's t-test (unequal variance)
            test_name = "Welch's t-test (unequal variance)"
            statistic, p_value = stats.ttest_ind(g1, g2, equal_var=False)
        # Effect size: Cohen's d for independent samples
        pooled_std = np.sqrt(((len(g1) - 1) * g1_std**2 + (len(g2) - 1) * g2_std**2) /
                             (len(g1) + len(g2) - 2))
        effect_size = (g1_mean - g2_mean) / pooled_std
        effect_size_name = "Cohen's d"
    else:
        # Mann-Whitney U test
        test_name = "Mann-Whitney U test"
        statistic, p_value = stats.mannwhitneyu(g1, g2, alternative='two-sided')
        # Effect size: rank-biserial correlation
        u1 = statistic
        u2 = len(g1) * len(g2) - u1
        u_max = max(u1, u2)
        effect_size = 1 - (2 * u_max) / (len(g1) * len(g2))
        effect_size_name = "rank-biserial r"

    # Determine significance
    significant = p_value <= alpha

    # Effect size interpretation
    abs_effect = abs(effect_size)
    if abs_effect < 0.2:
        effect_interp = "negligible"
    elif abs_effect < 0.5:
        effect_interp = "small"
    elif abs_effect < 0.8:
        effect_interp = "medium"
    else:
        effect_interp = "large"

    # Build interpretation
    if significant:
        interpretation = (
            f"There IS a statistically significant difference between the two groups "
            f"(p = {p_value:.4f} <= {alpha}). "
            f"The effect size is {effect_interp} ({effect_size_name} = {effect_size:.4f})."
        )
    else:
        interpretation = (
            f"There is NO statistically significant difference between the two groups "
            f"(p = {p_value:.4f} > {alpha}). "
            f"The effect size is {effect_interp} ({effect_size_name} = {effect_size:.4f})."
        )

    # Prepare result
    result = {
        'test_name': test_name,
        'statistic': float(statistic),
        'p_value': float(p_value),
        'significant': significant,
        'alpha': alpha,
        'group1_mean': float(g1_mean),
        'group2_mean': float(g2_mean),
        'group1_std': float(g1_std),
        'group2_std': float(g2_std),
        'effect_size': float(effect_size),
        'effect_size_name': effect_size_name,
        'effect_size_interpretation': effect_interp,
        'interpretation': interpretation,
        'assumptions': assumptions
    }

    # Log detailed output if verbose
    if verbose:
        logger.info("=" * 70)
        logger.info("STATISTICAL COMPARISON OF TWO GROUPS")
        logger.info("=" * 70)
        logger.info(f"\nSample sizes: Group 1 = {len(g1)}, Group 2 = {len(g2)}")
        logger.info(f"Group 1: Mean = {g1_mean:.4f}, SD = {g1_std:.4f}")
        logger.info(f"Group 2: Mean = {g2_mean:.4f}, SD = {g2_std:.4f}")
        logger.info(f"\nSample type: {'Paired' if use_paired else 'Independent'}")
        logger.info(f"Significance level (alpha): {alpha}")

        logger.info("\n" + "-" * 70)
        logger.info("ASSUMPTION CHECKS")
        logger.info("-" * 70)
        logger.info(f"Group 1 normality: {g1_normal_msg}")
        logger.info(f"Group 2 normality: {g2_normal_msg}")
        if assumptions['equal_variance']['p_value'] is not None:
            logger.info(f"Equal variance: {assumptions['equal_variance']['message']}")

        logger.info("\n" + "-" * 70)
        logger.info("TEST RESULTS")
        logger.info("-" * 70)
        logger.info(f"Test used: {test_name}")
        logger.info(f"Test statistic: {statistic:.4f}")
        logger.info(f"p-value: {p_value:.4f}")
        logger.info(f"{effect_size_name}: {effect_size:.4f} ({effect_interp} effect)")

        logger.info("\n" + "-" * 70)
        logger.info("CONCLUSION")
        logger.info("-" * 70)
        logger.info(interpretation)
        logger.info("=" * 70)

    return result


def compare_groups_simple(
    group1: List[float],
    group2: List[float],
    paired: bool = False,
    alpha: float = 0.05
) -> Tuple[bool, float, str]:
    """Simplified interface to compare two groups.

    Parameters
    ----------
    group1 : list of float
        First group of numerical data.
    group2 : list of float
        Second group of numerical data.
    paired : bool, optional
        If True, perform a paired test (default: False).
    alpha : float, optional
        Significance level (default: 0.05).

    Returns
    -------
    significant : bool
        True if there is a significant difference.
    p_value : float
        The p-value from the test.
    message : str
        Brief interpretation of the result.
    """
    result = compare_groups(group1, group2, paired=paired, alpha=alpha, verbose=False)
    return result['significant'], result['p_value'], result['interpretation']


# ---------------------------------------------------------------------------
# Count-based tests (from bar_significance.py)
# ---------------------------------------------------------------------------
def compute_pvalue(
    table: np.ndarray,
    method: str = "chi2"
) -> float:
    """Compute p-value for a 2x2 contingency table or two-sample comparison.

    Parameters
    ----------
    table : np.ndarray
        Input data:
        - For "chi2" / "fisher": shape (2, 2)
        - For "t-test" / "mannwhitney": shape (2, n)
    method : str, default="chi2"
        Statistical test method. Supported:
        - "chi2" : Chi-square test
        - "fisher" : Fisher's exact test
        - "t-test" : Independent t-test
        - "mannwhitney" : Mann-Whitney U test

    Returns
    -------
    float
        P-value.

    Raises
    ------
    ValueError
        If the method is unsupported.
    """
    if method == "chi2":
        _, p, _, _ = stats.chi2_contingency(table)
    elif method == "fisher":
        _, p = stats.fisher_exact(table)
    elif method == "t-test":
        p = stats.ttest_ind(table[0], table[1], equal_var=False).pvalue
    elif method == "mannwhitney":
        p = stats.mannwhitneyu(table[0], table[1], alternative="two-sided").pvalue
    else:
        raise ValueError(f"Unsupported method: {method}")
    return p


def compute_significance(
    pivot: pd.DataFrame,
    group_order: Tuple[str, str],
    method: str = "chi2"
) -> Dict[str, str]:
    """Compute significance stars per feature row of a count pivot table.

    Parameters
    ----------
    pivot : pandas.DataFrame
        Pivot table with features as index and groups as columns
        (one summarized count per feature x group cell).
    group_order : tuple of str
        Two groups to compare.
    method : str, default="chi2"
        'chi2' or 'fisher'.

    Returns
    -------
    dict
        Mapping: feature -> significance stars.
    """
    stars: Dict[str, str] = {}
    g1, g2 = group_order

    for feature in pivot.index:
        table = np.array([
            [pivot.loc[feature, g1], pivot[g1].sum() - pivot.loc[feature, g1]],
            [pivot.loc[feature, g2], pivot[g2].sum() - pivot.loc[feature, g2]],
        ])
        p = compute_pvalue(table, method)
        stars[feature] = p_to_star(p)

    return stars


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
                res = compare_groups(
                    g1_vals, g2_vals, paired=paired, alpha=alpha,
                    verbose=False, method=method,
                )
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
# Bracket drawing (display-point layout, from bar_significance.py)
# ---------------------------------------------------------------------------
def y_to_px(ax: plt.Axes, y_data: float) -> float:
    """Map a data y value on ``ax`` to display coordinates (pixels).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes whose transform is used.
    y_data : float
        Data-coordinate y value.

    Returns
    -------
    float
        Display-coordinate y (pixels).
    """
    return float(ax.transData.transform((0.0, y_data))[1])


def px_to_y(ax: plt.Axes, y_disp: float) -> float:
    """Map a display-coordinate y (pixels) to a data value on ``ax``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes whose transform is used.
    y_disp : float
        Display-coordinate y (pixels).

    Returns
    -------
    float
        Data-coordinate y value.
    """
    return float(ax.transData.inverted().transform((0.0, y_disp))[1])


def fit_bracket_headroom(ax: plt.Axes, bases: Sequence[float], n_brackets: int) -> None:
    """Raise the y upper limit so the tallest bracket stack fits.

    Finalize the axes layout (ticks, legend, tight_layout, canvas.draw)
    BEFORE calling this: display transforms must be stable afterwards.
    The new top is bisected so that the tallest base lands exactly the
    stack height below the axes-top pixel. Only ``ylim[1]`` moves, so
    this works identically on linear and log axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes already containing the bars and error bars.
    bases : array-like of float
        Per-cluster y of the tallest bar top including its error bar.
    n_brackets : int
        Maximum number of stacked brackets in any cluster.

    Returns
    -------
    None
    """
    if not n_brackets:
        return
    fig = ax.get_figure()
    px_per_pt = fig.dpi / 72.0
    stack_px = (LEG_PT + (n_brackets - 1) * BRACKET_STEP_PT + TEXT_PT + FONT_ASCENT_PT + 2.0) * px_per_pt
    positives = [float(b) for b in bases if b > 0]
    if not positives:
        return
    base = max(positives)
    axes_top_px = y_to_px(ax, ax.get_ylim()[1])  # top of the axes box
    target_px = axes_top_px - stack_px
    if y_to_px(ax, base) <= target_px:
        return  # already fits
    bottom = ax.get_ylim()[0]
    is_log = ax.get_yscale() == "log"
    lo = ax.get_ylim()[1]
    hi = lo * 10.0 if is_log else lo + (lo - bottom)
    while y_to_px(ax, base) > target_px:  # expand until it fits
        lo = hi
        hi = hi * 10.0 if is_log else hi + (hi - bottom)
        ax.set_ylim(bottom=bottom, top=hi)
        if hi > 1e12 * lo:
            break
    for _ in range(60):  # bisect the top value
        mid = np.sqrt(lo * hi) if is_log else (lo + hi) / 2.0
        ax.set_ylim(bottom=bottom, top=mid)
        if y_to_px(ax, base) > target_px:
            lo = mid
        else:
            hi = mid
    ax.set_ylim(bottom=bottom, top=hi)


def draw_significance_bracket(
    ax: plt.Axes,
    x1: float,
    x2: float,
    y_base: float,
    label: str,
    level: int = 0,
    lw: float = 0.7,
    fontsize: float = 8.0,
    color: str = "#333333",
) -> float:
    """Draw one bracket-with-ticks significance annotation above a bar cluster.

    The bracket is placed ``LEG_PT + level * BRACKET_STEP_PT`` display
    points above ``y_base`` (screen units, so spacing is identical on
    linear and log axes). The horizontal line spans ``x1`` to ``x2``
    (bar centers) with small outward diagonal ticks at both ends, drawn
    as a single polyline to avoid line-cap overlap at junctions.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes (layout already finalized, see fit_bracket_headroom).
    x1, x2 : float
        Data x of the left/right bar centers.
    y_base : float
        Data y of the cluster's tallest bar top including its error bar.
    label : str
        Star label (from p_to_star or custom text).
    level : int, default=0
        Stack level (0 = closest to the bars).
    lw : float, default=0.7
        Bracket line width.
    fontsize : float, default=8.0
        Star label font size.
    color : str, default="#333333"
        Bracket and label color.

    Returns
    -------
    float
        Display y (pixels) of the star label top, for further stacking.
    """
    fig = ax.get_figure()
    px_per_pt = fig.dpi / 72.0
    y_hat_disp = y_to_px(ax, y_base) + (LEG_PT + level * BRACKET_STEP_PT) * px_per_pt
    y_hat = px_to_y(ax, y_hat_disp)
    y_text = px_to_y(ax, y_hat_disp + TEXT_PT * px_per_pt)
    y_tick = px_to_y(ax, y_hat_disp - 3.0 * px_per_pt)
    ax.plot([x1, x1, x2, x2], [y_tick, y_hat, y_hat, y_tick],
            lw=lw, color=color, clip_on=False)
    ax.text((x1 + x2) / 2, y_text, label,
            ha="center", va="bottom", fontsize=fontsize, color=color)
    return y_hat_disp + (TEXT_PT + FONT_ASCENT_PT) * px_per_pt


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
                   help="image format(s) for the bar figure (png/pdf/svg/...; "
                        "matplotlib savefig set). Repeatable, e.g. --formats png "
                        "--formats pdf. An explicit empty value (--formats '') "
                        "means no figure; default png when absent.")
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
    for fmt in formats:
        out_png = str(args.outdir / f"significance_bars.{fmt}")
        plot_grouped_bars(
            df,
            results,
            out_png,
            feature_col=args.feature_col,
            group_col=args.group_col,
            value_col=args.value_col,
            groups=args.groups,
            feature_order=results["feature"].tolist(),
            ylabel=args.ylabel,
            legend_labels=args.legend,
            annotate_ns=args.annotate_ns,
            show_errorbars=not args.no_errorbars,
            dpi=args.dpi,
        )
        logger.info(f"Plot saved: {out_png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
