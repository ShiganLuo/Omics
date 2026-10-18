#!/usr/bin/env python3
"""Summarize Fiber-seq 5mC BED12 calls over TE annotations and plot A/B/C.

Input 5mC BED12 columns:
chrom start end read_name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts

Important:
- For fibertools-rs v0.13.1, the first and last BED12 blocks are sentinels;
  only interior blocks are treated as 5mC calls.
- This script does NOT calculate methylation percentage: the supplied format
  has no explicit denominator for callable unmethylated CpGs.
- The normalized metric is 5mC calls per kb of read span overlapping a TE,
  not a methylation fraction.
- Summary tables include TE insertions overlapped by at least one read; TE
  insertions with zero read overlap are omitted to reduce runtime and output size.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import logging
import math
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from scipy import stats

logger = logging.getLogger("te_5mc")


def open_text(path: str):
    """Open plain-text or gzip-compressed text files."""
    path_obj = Path(path)
    if path_obj.suffix.lower() in {".gz", ".bgz", ".bgzip"}:
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


@dataclass
class TE:
    idx: int
    chrom: str
    start: int
    end: int
    name: str
    family: str
    cls: str
    strand: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Fiber-seq 5mC BED12 calls over TE annotations."
    )
    parser.add_argument("--te-bed", required=True, help="TE BED: chrom start end name score strand family_id class_id")
    parser.add_argument(
        "--sample", action="append", required=True, metavar="LABEL=BED12",
        help="Sample label and 5mC BED12 path; provide twice, e.g. --sample Pluri=a.bed --sample Toti=b.bed",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--min-overlap-reads", type=int, default=5,
                        help="Minimum read-TE overlaps to retain an insertion (default: 5)")
    parser.add_argument("--top-subfamilies", type=int, default=50,
                        help="Number of subfamilies in heatmap (default: 50)")
    parser.add_argument("--top-loci", type=int, default=1,
                        help="Number of representative loci for panel C (default: 1; currently plots top locus)")
    parser.add_argument("--max-reads-per-locus", type=int, default=250,
                        help="Maximum read rows per sample in panel C (default: 250)")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return parser.parse_args()

def configure_chinese_font() -> str:
    """Configure a CJK-capable font for matplotlib and return selected font name."""
    candidates = [
        "Noto Sans CJK SC",
        "Source Han Sans SC",
        "WenQuanYi Micro Hei",
        "Microsoft YaHei",
        "SimHei",
        "PingFang SC",
        "Heiti SC",
    ]
    available_fonts = {f.name for f in font_manager.fontManager.ttflist}
    selected = ""
    for font_name in candidates:
        if font_name in available_fonts:
            selected = font_name
            break

    if selected:
        plt.rcParams["font.sans-serif"] = [selected, "DejaVu Sans", "Arial Unicode MS"]
    else:
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial Unicode MS"]

    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.family"] = "sans-serif"
    return selected

def apply_font_after_style(selected_font: str) -> None:
    """Re-apply font settings after sns.set_style() which resets rcParams."""
    if selected_font:
        plt.rcParams["font.sans-serif"] = [selected_font, "DejaVu Sans", "Arial Unicode MS"]
        plt.rcParams["axes.unicode_minus"] = False
        plt.rcParams["font.family"] = "sans-serif"

def parse_samples(items: List[str]) -> List[Tuple[str, str]]:
    samples = []
    seen = set()
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --sample {item!r}; expected LABEL=BED12")
        label, path = item.split("=", 1)
        label, path = label.strip(), path.strip()
        if not label or not path:
            raise ValueError(f"Invalid --sample {item!r}; expected LABEL=BED12")
        if label in seen:
            raise ValueError(f"Duplicate sample label: {label}")
        seen.add(label)
        if not Path(path).is_file():
            raise FileNotFoundError(path)
        samples.append((label, path))
    if len(samples) != 2:
        raise ValueError("Provide exactly two --sample arguments (one pluripotent and one totipotent sample).")
    return samples


def load_te_bed(path: str) -> Tuple[List[TE], Dict[str, List[TE]], Dict[str, List[int]], Dict[str, List[int]]]:
    """Load TE intervals and index them by chromosome, sorted by start."""
    tes: List[TE] = []
    by_chrom: Dict[str, List[TE]] = defaultdict(list)
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 8:
                logger.warning("Skipping TE BED line %d with fewer than 8 columns", line_no)
                continue
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            except ValueError:
                logger.warning("Skipping malformed TE BED line %d", line_no)
                continue
            if end <= start:
                continue
            te = TE(
                idx=len(tes), chrom=chrom, start=start, end=end,
                name=fields[3], family=fields[6], cls=fields[7],
                strand=fields[5],
            )
            tes.append(te)
            by_chrom[chrom].append(te)
    starts: Dict[str, List[int]] = {}
    prefix_max_end: Dict[str, List[int]] = {}
    for chrom, intervals in by_chrom.items():
        intervals.sort(key=lambda x: (x.start, x.end))
        starts[chrom] = [x.start for x in intervals]
        running = 0
        prefix = []
        for te in intervals:
            running = max(running, te.end)
            prefix.append(running)
        prefix_max_end[chrom] = prefix
    if not tes:
        raise ValueError(f"No TE intervals found in {path}")
    logger.info("Loaded %d TE intervals on %d chromosomes", len(tes), len(by_chrom))
    return tes, by_chrom, starts, prefix_max_end


def overlapping_tes(
    chrom: str, start: int, end: int,
    by_chrom: Dict[str, List[TE]], starts: Dict[str, List[int]],
    prefix_max_end: Dict[str, List[int]]
) -> Iterable[TE]:
    """Yield TE intervals overlapping [start, end); intervals are BED 0-based."""
    intervals = by_chrom.get(chrom)
    if not intervals or end <= start:
        return
    stop = bisect.bisect_left(starts[chrom], end)
    i = stop - 1
    # Walk backwards only while an earlier interval could still overlap start.
    while i >= 0 and prefix_max_end[chrom][i] > start:
        te = intervals[i]
        if te.end > start:
            yield te
        i -= 1


def parse_bed12(line: str, line_no: int) -> Optional[Tuple[str, int, int, str, List[int]]]:
    """Parse v0.13.1 ft extract -c BED12 and return real 5mC call positions.

    In this output convention, the first and last blocks are boundary sentinels:
    the first has size 0 and the last has size 1. Neither is a 5mC call.
    Actual calls are represented by the interior blocks.
    """
    fields = line.rstrip("\n").split()
    if len(fields) < 12:
        return None
    try:
        chrom = fields[0]
        start, end = int(fields[1]), int(fields[2])
        read_name = fields[3]
        n_blocks = int(fields[9])
        block_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x != ""]
        block_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x != ""]
    except (ValueError, IndexError):
        return None
    if start < 0 or end <= start or n_blocks < 2:
        return None
    if len(block_sizes) != n_blocks or len(block_starts) != n_blocks:
        return None

    # ft extract v0.13.1 writes two sentinel blocks, one at each read boundary.
    # Exclude both explicitly; filtering only size==0 would retain the terminal sentinel.
    positions = []
    for offset, size in zip(block_starts[1:-1], block_sizes[1:-1]):
        if size != 1:
            continue
        pos = start + offset
        if start <= pos < end:
            positions.append(pos)
    return chrom, start, end, read_name, positions


def iter_bed12(path: str):
    """Stream BED12 records; no whole-file loading."""
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parsed = parse_bed12(line, line_no)
            if parsed is None:
                continue
            yield parsed


def analyze_sample(
        label: str, path: str, tes: List[TE],
        by_chrom: Dict[str, List[TE]], starts: Dict[str, List[int]],
        prefix_max_end: Dict[str, List[int]]
    ) -> pd.DataFrame:
    """Aggregate only TE insertions overlapped by reads.

    Sparse counters avoid allocating a Python dictionary for every TE interval.
    Call counts are computed from the read's sorted call positions, avoiding a
    second interval-index lookup for every individual 5mC call.
    """
    # idx -> [5mC calls, overlapping reads, read-overlap bp]
    stats: Dict[int, List[int]] = {}
    total_records = 0
    malformed = 0

    for chrom, read_start, read_end, read_name, positions in iter_bed12(path):
        total_records += 1
        touched = overlapping_tes(
            chrom, read_start, read_end, by_chrom, starts, prefix_max_end
        )
        any_touched = False
        for te in touched:
            any_touched = True
            item = stats.get(te.idx)
            if item is None:
                item = [0, 0, 0]
                stats[te.idx] = item
            item[1] += 1
            overlap_start = max(read_start, te.start)
            overlap_end = min(read_end, te.end)
            overlap_bp = overlap_end - overlap_start
            if overlap_bp <= 0:
                continue
            item[2] += overlap_bp
            # positions are sorted; count calls within this TE/read intersection
            # in O(log n) rather than querying the TE interval index per call.
            left = bisect.bisect_left(positions, overlap_start)
            right = bisect.bisect_left(positions, overlap_end, lo=left)
            item[0] += right - left

        if total_records % 5_000_000 == 0:
            logger.info("%s: processed %s BED12 records", label, f"{total_records:,}")

    logger.info(
        "%s: processed %s BED12 records; %s TE insertions overlapped",
        label, f"{total_records:,}", f"{len(stats):,}"
    )

    rows = []
    for te_idx, item in stats.items():
        te = tes[te_idx]
        n_calls, n_reads, overlap_bp = item
        rows.append({
            "sample": label,
            "te_idx": te.idx,
            "chrom": te.chrom,
            "start": te.start,
            "end": te.end,
            "locus": f"{te.chrom}:{te.start}-{te.end}",
            "subfamily": te.name,
            "family": te.family,
            "class": te.cls,
            "strand": te.strand,
            "length_bp": te.end - te.start,
            "n_5mc_calls": n_calls,
            "n_read_overlaps": n_reads,
            "read_overlap_bp": overlap_bp,
            "calls_per_kb_read_overlap": (
                n_calls / (overlap_bp / 1000.0) if overlap_bp > 0 else np.nan
            ),
            "calls_per_read_overlap": n_calls / n_reads if n_reads > 0 else np.nan,
        })
    return pd.DataFrame.from_records(rows)


def aggregate_level(df: pd.DataFrame, keys: List[str]) -> pd.DataFrame:
    """Aggregate calls and read-overlap bases, then calculate a pooled call density."""
    grouped = df.groupby(keys, dropna=False).agg(
        n_5mc_calls=("n_5mc_calls", "sum"),
        n_read_overlaps=("n_read_overlaps", "sum"),
        read_overlap_bp=("read_overlap_bp", "sum"),
        n_insertions=("te_idx", "nunique"),
    ).reset_index()
    grouped["calls_per_kb_read_overlap"] = np.where(
        grouped["read_overlap_bp"] > 0,
        grouped["n_5mc_calls"] / (grouped["read_overlap_bp"] / 1000.0),
        np.nan,
    )
    grouped["calls_per_read_overlap"] = np.where(
        grouped["n_read_overlaps"] > 0,
        grouped["n_5mc_calls"] / grouped["n_read_overlaps"],
        np.nan,
    )
    return grouped

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

def plot_comparison(
        ax: matplotlib.axes.Axes,
        result: Dict[str, Any],
        title: str = "Mapping Rate Comparison",
        xlabel: str = "Sample Type",
        ylabel: str = "Mapping Rate (%)",
        palette: Tuple[str, str] = ("#3498db", "#e74c3c"),
        show_points: bool = True,
        show_stats: bool = True,
        save_path: Optional[str] = None,
        dpi: int = 300,
        verbose: bool = True
    ) -> matplotlib.axes.Axes:
    """
    Create a publication-quality comparison plot for two groups.

    Parameters
    ----------
    result : dict
        Result dictionary from prepare_data() containing:
        - group1_name, group2_name: Names of the groups
        - group1_data, group2_data: Data values for each group
        - comparison: Statistical comparison results
    title : str, optional
        Plot title (default: "Mapping Rate Comparison").
    xlabel : str, optional
        X-axis label (default: "Sample Type").
    ylabel : str, optional
        Y-axis label (default: "Mapping Rate (%)").
    figsize : tuple of int, optional
        Figure size (default: (10, 7)).
    palette : tuple of str, optional
        Colors for the two groups (default: ("#3498db", "#e74c3c")).
    show_points : bool, optional
        Whether to show individual data points (default: True).
    show_stats : bool, optional
        Whether to show statistics annotation (default: True).
    save_path : str, optional
        Path to save the figure (default: None, does not save).
    dpi : int, optional
        DPI for saved figure (default: 300).
    verbose : bool, optional
        Whether to logger.info information about the plot (default: True).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure object.
    """
    selected_font = configure_chinese_font()
    if verbose and selected_font:
        logger.info(f"Using Chinese font: {selected_font}")

    # Extract data from result
    group1_name = result['group1_name']
    group2_name = result['group2_name']
    group1_data = result['group1_data']
    group2_data = result['group2_data']
    comparison = result['comparison']

    # Set style
    sns.set_style("whitegrid")
    apply_font_after_style(selected_font)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    
    # Prepare data for plotting
    data_to_plot = [group1_data, group2_data]
    positions = [1, 2]

    # Create box plot
    bp = ax.boxplot(
        data_to_plot,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showmeans=True,
        meanprops={"marker": "D", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": 8},
        medianprops={"color": "black", "linewidth": 2},
        whiskerprops={"color": "black", "linewidth": 1.5},
        capprops={"color": "black", "linewidth": 1.5},
        flierprops={"marker": "o", "markerfacecolor": "gray", "markersize": 6, "alpha": 0.5}
    )

    # Color the boxes
    for patch, color in zip(bp['boxes'], palette):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)

    # Add individual data points if requested
    if show_points:
        for i, (data, pos, color) in enumerate(zip(data_to_plot, positions, palette)):
            # Add jitter to x-coordinates
            jitter = np.random.uniform(-0.15, 0.15, size=len(data))
            ax.scatter(
                np.full_like(data, pos) + jitter,
                data,
                color=color,
                alpha=0.6,
                s=50,
                edgecolors='black',
                linewidths=0.5,
                zorder=5
            )

    # Set x-axis labels with explicit font properties
    if selected_font:
        font_prop = FontProperties(family=selected_font, size=12)
        font_prop_bold = FontProperties(family=selected_font, size=14, weight='bold')
        font_prop_title = FontProperties(family=selected_font, size=16, weight='bold')
    else:
        font_prop = FontProperties(size=12)
        font_prop_bold = FontProperties(size=14, weight='bold')
        font_prop_title = FontProperties(size=16, weight='bold')

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{group1_name}\n(n={len(group1_data)})", f"{group2_name}\n(n={len(group2_data)})"],
                       fontproperties=font_prop)

    # Set axis labels
    ax.set_xlabel(xlabel, fontproperties=font_prop_bold)
    ax.set_ylabel(ylabel, fontproperties=font_prop_bold)
    ax.set_title(title, fontproperties=font_prop_title, pad=20)

    # Add statistics annotation if requested
    if show_stats:
        # Get statistical results
        p_value = comparison['p_value']
        significant = comparison['significant']
        effect_size = comparison['effect_size']
        effect_interp = comparison['effect_size_interpretation']

        # Determine significance stars
        if p_value < 0.001:
            sig_stars = "***"
        elif p_value < 0.01:
            sig_stars = "**"
        elif p_value < 0.05:
            sig_stars = "*"
        else:
            sig_stars = "ns"

        # Calculate position for annotation
        y_max = max(max(group1_data), max(group2_data))
        y_min = min(min(group1_data), min(group2_data))
        y_range = y_max - y_min

        # Position the bracket
        bracket_y = y_max + y_range * 0.15
        star_y = bracket_y + y_range * 0.05

        # Draw bracket
        ax.plot([1, 1, 2, 2], [bracket_y, bracket_y + y_range * 0.02, bracket_y + y_range * 0.02, bracket_y],
                color='black', linewidth=1.5)

        # Add significance stars
        ax.text(1.5, star_y, sig_stars, ha='center', va='bottom', fontsize=16, fontweight='bold')

        # Add p-value text
        p_text = f"p = {p_value:.4f}" if p_value >= 0.001 else f"p < 0.001"
        ax.text(1.5, star_y - y_range * 0.08, p_text, ha='center', va='top', fontsize=11, style='italic')

        # Add effect size info - use axes transform for better positioning
        effect_text = f"{comparison['effect_size_name']}: {effect_size:.2f} ({effect_interp})"
        ax.text(0.98, 0.98, effect_text, transform=ax.transAxes, fontsize=10, 
                color='gray', ha='right', va='top', style='italic')

        # Add test name - use axes transform for better positioning
        test_text = f"Test: {comparison['test_name']}"
        ax.text(0.02, 0.02, test_text, transform=ax.transAxes, fontsize=9, color='gray', style='italic')

    # Add mean and SD annotations
    for i, (data, pos, name) in enumerate(zip(data_to_plot, [0,3], [group1_name, group2_name])):
        mean_val = np.mean(data)
        std_val = np.std(data, ddof=1)
        ax.annotate(
            f"Mean: {mean_val:.1f}%\nSD: {std_val:.1f}%",
            xy=(pos, np.median(data)),
            xytext=(pos + 0.4 if i == 0 else pos - 0.4, np.median(data)),
            fontsize=9,
            ha='left' if i == 0 else 'right',
            va='center',
            arrowprops=dict(arrowstyle='->', color='gray', lw=1),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8)
        )

    # Adjust y-axis limits to make room for annotation
    ax.set_ylim(y_min - y_range * 0.1, y_max + y_range * 0.35)
    
    # Add grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Tight layout
    plt.tight_layout()

    if verbose:
        logger.info(f"Plot created successfully")
        logger.info(f"  Groups: {group1_name} (n={len(group1_data)}) vs {group2_name} (n={len(group2_data)})")
        logger.info(f"  p-value: {p_value:.4f} ({'significant' if significant else 'not significant'})")

    return ax

def plot_a(all_df: pd.DataFrame, labels: List[str], outdir: str, dpi: int, level: Literal["class", "family", "subfamily"] = "family") -> None:
    """A: class/family distribution of called 5mC density."""
    # Calculate density per insertion first, then show distribution across insertions.
    valid = all_df[
        (all_df["n_read_overlaps"] > 0) &
        all_df["calls_per_kb_read_overlap"].notna()
    ].copy()

    if valid.empty:
        logger.warning("No TE insertions with read overlap for panel A")
        return
    classes = valid.groupby(level)["n_read_overlaps"].sum().sort_values(ascending=False).index.tolist()
    metrics = [
        ("calls_per_kb_read_overlap", "5mC calls per kb read-overlap"),
        ("calls_per_read_overlap", "5mC calls per overlapping read"),
    ]

    for cls in classes:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        cls_data = valid.loc[valid[level] == cls]

        for ax, (metric, ylabel) in zip(axes, metrics):
            comparison_data = defaultdict(list)
            pos = 1

            for label in labels:
                values = (
                    cls_data.loc[cls_data["sample"] == label, metric]
                    .replace([np.inf, -np.inf], np.nan)
                    .dropna()
                    .values
                )

                comparison_data[f"group{pos}_name"] = label
                comparison_data[f"group{pos}_data"] = list(values)
                pos += 1
            if len(comparison_data["group1_data"]) < 20 or len(comparison_data["group2_data"]) < 20:
                logger.warning(
                    "Insufficient data for statistical comparison in class %s, metric %s: "
                    "group1 n=%d, group2 n=%d",
                    cls, metric, len(comparison_data["group1_data"]), len(comparison_data["group2_data"])
                )
                continue
            
            stats_result = compare_groups(
                comparison_data["group1_data"],
                comparison_data["group2_data"],
                paired=False,
                alpha=0.05,
                verbose=True,
            )

            comparison_data["comparison"] = stats_result
            plot_comparison(
                ax,
                comparison_data,
                title="",
                xlabel="Sample Type",
                ylabel=ylabel,
                dpi=dpi,
                show_points=False
            )

        fig.suptitle(f"Fiber-seq 5mC call density — {cls}")
        fig.tight_layout()

        safe_cls = str(cls).replace("/", "_").replace(" ", "_")
        fig.savefig(
            os.path.join(outdir, f"A_te_class_5mc_density_{safe_cls}.png"),
            dpi=dpi,
            bbox_inches="tight",
        )
        plt.close(fig)

    fam = aggregate_level(all_df, ["sample", "family"])
    fam.to_csv(
        os.path.join(outdir, "A_family_summary.tsv"),
        sep="\t",
        index=False,
    )


def plot_b(all_df: pd.DataFrame, labels: List[str], outdir: str, top_n: int, dpi: int) -> None:
    """B: heatmap of top subfamilies; values are log2 call density, not methylation %."""
    valid = all_df[all_df["n_read_overlaps"] > 0].copy()
    if valid.empty:
        logger.warning("No TE insertions with read overlap for panel B")
        return
    sub = aggregate_level(valid, ["sample", "subfamily"])
    sub["log2_density"] = np.log2(sub["calls_per_kb_read_overlap"] + 1e-9)
    pivot = sub.pivot(index="subfamily", columns="sample", values="log2_density")
    pivot = pivot.reindex(columns=labels)
    # Select subfamilies with the greatest pooled read-overlap bases, avoiding low-coverage labels.
    exposure = valid.groupby("subfamily")["read_overlap_bp"].sum()
    eligible = exposure[exposure > 0].index
    pivot = pivot.loc[pivot.index.intersection(eligible)]
    rank = exposure.reindex(pivot.index).sort_values(ascending=False)
    selected = rank.head(max(1, top_n)).index
    pivot = pivot.reindex(selected)
    pivot.to_csv(os.path.join(outdir, "B_subfamily_heatmap_matrix.tsv"), sep="\t", na_rep="NA")
    if pivot.empty:
        return
    # With only two samples, row-wise z-scores collapse to nearly +/-1 and are misleading.
    # Plot the log2 normalized call-density values directly instead.
    values = pivot.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size:
        vmin, vmax = np.percentile(finite, [2, 98])
        if math.isclose(float(vmin), float(vmax)):
            vmin, vmax = float(np.min(finite)), float(np.max(finite)) + 1e-6
    else:
        vmin, vmax = -1.0, 1.0
    fig_h = max(5, 0.22 * len(pivot))
    fig, ax = plt.subplots(figsize=(7, fig_h))
    im = ax.imshow(values, aspect="auto", interpolation="nearest", cmap="viridis", vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_title("B. TE subfamilies (log2 5mC call density)")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("log2(calls per kb read-overlap + 1e-9)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "B_te_subfamily_heatmap.png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def collect_locus_reads(path: str, te: TE, max_reads: int) -> List[Tuple[str, List[int]]]:
    """Collect read-level call positions relative to a TE insertion."""
    rows = []
    for chrom, read_start, read_end, read_name, positions in iter_bed12(path):
        if chrom != te.chrom or read_end <= te.start or read_start >= te.end:
            continue
        relative = [p - te.start for p in positions if te.start <= p < te.end]
        # Include overlapping reads even when no called sites fall inside this TE.
        rows.append((read_name, relative))
    if len(rows) > max_reads:
        # Deterministic sampling for display only; summary statistics use all reads.
        idx = np.linspace(0, len(rows) - 1, max_reads, dtype=int)
        rows = [rows[i] for i in idx]
    return rows


def plot_c(
    all_df: pd.DataFrame, samples: List[Tuple[str, str]], tes: List[TE],
    outdir: str, max_reads: int, dpi: int
) -> None:
    """C: read-level 5mC positions across the highest-call TE insertion."""
    agg = all_df.groupby("te_idx").agg(
        calls=("n_5mc_calls", "sum"),
        reads=("n_read_overlaps", "sum"),
    )
    agg = agg[agg["reads"] > 0]
    if agg.empty:
        logger.warning("No eligible TE insertion for panel C")
        return
    best_idx = int(agg.sort_values(["calls", "reads"], ascending=False).index[0])
    te = tes[best_idx]
    width = te.end - te.start
    fig, axes = plt.subplots(len(samples), 1, figsize=(12, max(4, len(samples) * 3.2)), sharex=True)
    if len(samples) == 1:
        axes = [axes]
    for ax, (label, path) in zip(axes, samples):
        rows = collect_locus_reads(path, te, max_reads)
        if not rows:
            ax.text(0.5, 0.5, "No reads overlap this TE in this sample",
                    transform=ax.transAxes, ha="center", va="center")
            ax.set_ylabel(label)
            continue
        for y, (read_name, rel_positions) in enumerate(rows):
            ax.hlines(y, 0, width, linewidth=0.35, alpha=0.55)
            if rel_positions:
                ax.scatter(rel_positions, [y] * len(rel_positions), s=5, marker="|")
        ax.set_ylim(-1, len(rows))
        ax.invert_yaxis()
        ax.set_ylabel(f"{label}\nRead")
        ax.text(0.99, 0.96, f"shown reads: {len(rows)}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8)
        ax.grid(axis="x", alpha=0.2)
    axes[-1].set_xlabel("Position relative to TE start (bp)")
    fig.suptitle(
        f"C. Single-molecule 5mC calls at {te.chrom}:{te.start}-{te.end} "
        f"({te.name}; {te.family}; {te.cls}; strand {te.strand})"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "C_representative_TE_single_molecule.png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("Panel C locus: %s:%d-%d (%s)", te.chrom, te.start, te.end, te.name)


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    samples = parse_samples(args.sample)
    labels = [x[0] for x in samples]
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    tes, by_chrom, starts, prefix_max_end = load_te_bed(args.te_bed)
    dfs = []
    outfile_insertions = os.path.join(outdir, "all_samples.te_insertion_summary.tsv")
    if os.path.exists(outfile_insertions) and not os.path.getsize(outfile_insertions) == 0:
        logger.info("Skipping all samples: %s (output exists)", outfile_insertions)
        all_df = pd.read_csv(outfile_insertions, sep="\t")
        dfs.append(all_df)
    else:
        for label, path in samples:
            outfile = os.path.join(outdir, f"{label}.te_insertion_summary.tsv")
            if os.path.exists(outfile) and not os.path.getsize(outfile) == 0:
                logger.info("Skipping sample %s: %s (output exists)", label, path)
                df = pd.read_csv(outfile, sep="\t")
                dfs.append(df)
            else:
                logger.info("Analyzing sample %s: %s", label, path)
                df = analyze_sample(label, path, tes, by_chrom, starts, prefix_max_end)
                dfs.append(df)
                df.to_csv(outfile, sep="\t", index=False)
        all_df = pd.concat(dfs, ignore_index=True)
        all_df.to_csv(outfile_insertions, sep="\t", index=False)

    # Remove low-support insertions for display/summary tables; raw per-insertion table is retained.
    plot_df = all_df[all_df["n_read_overlaps"] >= args.min_overlap_reads].copy()
    outfile_filtered = os.path.join(outdir, "filtered.te_insertion_summary.tsv")
    plot_df.to_csv(outfile_filtered, sep="\t", index=False)
    if plot_df.empty:
        logger.warning("No insertions pass --min-overlap-reads=%d; using all insertions for plots",
                    args.min_overlap_reads)
        plot_df = all_df.copy()

    plot_a(plot_df, labels, outdir, args.dpi, level = "locus")
    plot_b(plot_df, labels, outdir, args.top_subfamilies, args.dpi)
    plot_c(all_df, samples, tes, outdir, args.max_reads_per_locus, args.dpi)

    logger.info("Finished. Output directory: %s", outdir)


if __name__ == "__main__":
    main()