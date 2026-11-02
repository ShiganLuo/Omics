#!/usr/bin/env python3
# Author: ShiganLuo
# Email: 25303020102@qq.com
# GitHub: https://github.com/ShiganLuo
"""
Statistical significance test for comparing two groups of data.

This module provides functions to determine if there is a statistically
significant difference between two groups of numerical data.

Features:
- Automatic normality test (Shapiro-Wilk)
- Automatic homogeneity of variance test (Levene)
- Automatic selection of appropriate statistical test
- Support for both paired and independent samples
- Detailed output with test statistics and interpretation

Usage:
    from statistical_test import compare_groups

    group1 = [1.2, 2.3, 3.4, 4.5, 5.6]
    group2 = [2.3, 3.4, 4.5, 5.6, 6.7]

    result = compare_groups(group1, group2)
    logger.info(result)
"""

from __future__ import annotations
import sys
import os
from typing import List, Tuple, Optional, Dict, Any
import numpy as np
from scipy import stats
import logging
logger = logging.getLogger(__name__)

def check_normality(data: List[float], alpha: float = 0.05) -> Tuple[bool, float, str]:
    """
    Check if data follows a normal distribution using Shapiro-Wilk test.

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
    - If p <= alpha, we reject H0 (data does not appear normal).
    - Requires at least 3 data points.
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
    """
    Check if two groups have equal variances using Levene's test.

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

    Notes
    -----
    - H0: The groups have equal variances.
    - If p > alpha, we fail to reject H0 (variances appear equal).
    - If p <= alpha, we reject H0 (variances appear unequal).
    """
    stat, p_value = stats.levene(group1, group2)
    equal_variance = p_value > alpha

    if equal_variance:
        message = f"Groups appear to have equal variances (Levene p={p_value:.4f} > {alpha})"
    else:
        message = f"Groups appear to have unequal variances (Levene p={p_value:.4f} <= {alpha})"

    return equal_variance, p_value, message


def compare_groups(
    group1: List[float],
    group2: List[float],
    paired: bool = False,
    alpha: float = 0.05,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Compare two groups to determine if there is a statistically significant difference.

    This function automatically selects the appropriate statistical test based on
    the data characteristics (normality, equal variance).

    Parameters
    ----------
    group1 : list of float
        First group of numerical data.
    group2 : list of float
        Second group of numerical data.
    paired : bool, optional
        If True, perform a paired test (default: False).
        Use this when observations in group1 and group2 are related/matched.
    alpha : float, optional
        Significance level for hypothesis tests (default: 0.05).
    verbose : bool, optional
        If True, logger.info detailed information about the tests (default: True).

    Returns
    -------
    result : dict
        Dictionary containing:
        - 'test_name': Name of the statistical test used
        - 'statistic': Test statistic value
        - 'p_value': p-value from the test
        - 'significant': Boolean indicating if difference is significant (p <= alpha)
        - 'alpha': Significance level used
        - 'group1_mean': Mean of group1
        - 'group2_mean': Mean of group2
        - 'group1_std': Standard deviation of group1
        - 'group2_std': Standard deviation of group2
        - 'effect_size': Cohen's d effect size (for parametric tests) or
                         rank-biserial correlation (for non-parametric tests)
        - 'interpretation': Human-readable interpretation of results
        - 'assumptions': Dictionary with assumption check results

    Raises
    ------
    ValueError
        If either group has fewer than 2 data points.

    Examples
    --------
    >>> group1 = [1.2, 2.3, 3.4, 4.5, 5.6]
    >>> group2 = [2.3, 3.4, 4.5, 5.6, 6.7]
    >>> result = compare_groups(group1, group2)
    >>> logger.info(f"Significant difference: {result['significant']}")

    >>> # Paired sample test
    >>> before = [10, 12, 14, 16, 18]
    >>> after = [12, 14, 16, 18, 20]
    >>> result = compare_groups(before, after, paired=True)
    """
    # Input validation
    if len(group1) < 2 or len(group2) < 2:
        raise ValueError("Each group must have at least 2 data points")

    # Convert to numpy arrays
    g1 = np.array(group1, dtype=float)
    g2 = np.array(group2, dtype=float)

    # Calculate basic statistics
    g1_mean, g1_std = np.mean(g1), np.std(g1, ddof=1)
    g2_mean, g2_std = np.mean(g2), np.std(g2, ddof=1)

    # Check assumptions
    assumptions = {}

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

    # Select and perform appropriate test
    if paired:
        if both_normal:
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
    else:
        if both_normal:
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

    # logger.info detailed output if verbose
    if verbose:
        logger.info("=" * 70)
        logger.info("STATISTICAL COMPARISON OF TWO GROUPS")
        logger.info("=" * 70)
        logger.info(f"\nSample sizes: Group 1 = {len(g1)}, Group 2 = {len(g2)}")
        logger.info(f"Group 1: Mean = {g1_mean:.4f}, SD = {g1_std:.4f}")
        logger.info(f"Group 2: Mean = {g2_mean:.4f}, SD = {g2_std:.4f}")
        logger.info(f"\nSample type: {'Paired' if paired else 'Independent'}")
        logger.info(f"Significance level (α): {alpha}")

        logger.info("\n" + "-" * 70)
        logger.info("ASSUMPTION CHECKS")
        logger.info("-" * 70)
        logger.info(f"Group 1 normality: {g1_normal_msg}")
        logger.info(f"Group 2 normality: {g2_normal_msg}")
        if 'equal_variance' in assumptions and assumptions['equal_variance']['p_value'] is not None:
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
    """
    Simplified interface to compare two groups.

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

    Examples
    --------
    >>> group1 = [1.2, 2.3, 3.4, 4.5, 5.6]
    >>> group2 = [2.3, 3.4, 4.5, 5.6, 6.7]
    >>> sig, p, msg = compare_groups_simple(group1, group2)
    >>> logger.info(f"Significant: {sig}, p={p:.4f}")
    """
    result = compare_groups(group1, group2, paired=paired, alpha=alpha, verbose=False)
    return result['significant'], result['p_value'], result['interpretation']




if __name__ == "__main__":
    # Example usage
    logger.info("Example 1: Independent samples (normal data)")
    logger.info("-" * 70)
    np.random.seed(42)
    group1 = np.random.normal(loc=10, scale=2, size=30).tolist()
    group2 = np.random.normal(loc=12, scale=2, size=30).tolist()
    result1 = compare_groups(group1, group2)

    logger.info("\n\nExample 2: Paired samples")
    logger.info("-" * 70)
    before = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
    after = [12, 14, 15, 18, 19, 22, 23, 25, 28, 30]
    result2 = compare_groups(before, after, paired=True)

    logger.info("\n\nExample 3: Non-normal data (Mann-Whitney U)")
    logger.info("-" * 70)
    # Skewed data
    group3 = np.random.exponential(scale=2, size=20).tolist()
    group4 = np.random.exponential(scale=3, size=20).tolist()
    result3 = compare_groups(group3, group4)
