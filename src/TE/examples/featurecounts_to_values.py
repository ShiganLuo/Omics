#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert a locus-level featureCounts table to te_manhattan.py values.

Computes per-locus log2FC (mean CPM, group A vs group B) and a
Mann-Whitney p (BH-adjusted -> padj) from a featureCounts output table
built with ``-g transcript_id`` on the rmsk TE GTF. The result feeds
``te_manhattan.py --te-gtf <TE.gtf> --locus-col transcript_id``.

Usage
-----
python featurecounts_to_values.py --counts te_locus_counts.txt \
    --group-a ciTotiSC1-1,ciTotiSC1-2,... --group-b mESC,mESC1,mESC2 \
    -o locus_values.tsv
"""

import argparse
import os
import sys
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "modules", "function", "bin",
    ),
)
from LogUtil import setup_logger  # noqa: E402

logger = setup_logger(__name__)


def load_featurecounts(counts_tsv: str) -> pd.DataFrame:
    """Load a featureCounts output table (skip the comment header line).

    Parameters
    ----------
    counts_tsv : str
        featureCounts output path (``# Program:`` comment line + header
        starting with ``Geneid``).

    Returns
    -------
    pd.DataFrame
        Indexed by ``Geneid`` (transcript_id) with one numeric column
        per sample (basename of the BAM path).

    Raises
    ------
    ValueError
        If the header line is missing or no sample columns remain.
    """
    raw = pd.read_csv(counts_tsv, sep="\t", comment="#")
    if raw.empty or raw.columns[0] != "Geneid":
        raise ValueError(f"Not a featureCounts table (no 'Geneid' header): {counts_tsv}")
    sample_cols = [c for c in raw.columns if c not in ("Geneid", "Chr", "Start", "End", "Strand", "Length")]
    if not sample_cols:
        raise ValueError(f"No sample columns in {counts_tsv}")
    renamed = {c: os.path.splitext(os.path.basename(c))[0] for c in sample_cols}
    table = raw.set_index("Geneid")[sample_cols].rename(columns=renamed)
    table = table.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    logger.info(f"Loaded featureCounts table: {table.shape[0]} loci x {table.shape[1]} samples")
    return table


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment.

    Parameters
    ----------
    p_values : np.ndarray
        1-D raw p values (NaN allowed -> stays NaN).

    Returns
    -------
    np.ndarray
        Adjusted q values, same shape.
    """
    out = np.full(p_values.shape, np.nan, dtype=float)
    valid = ~np.isnan(p_values)
    p = p_values[valid]
    n = p.size
    if n == 0:
        return out
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = q
    out[valid] = restored
    return out


def compute_values(
    counts: pd.DataFrame,
    group_a: List[str],
    group_b: List[str],
    min_total: int = 10,
) -> pd.DataFrame:
    """Compute per-locus log2FC and Mann-Whitney padj (A vs B).

    Parameters
    ----------
    counts : pd.DataFrame
        Raw counts (loci x samples), see :func:`load_featurecounts`.
    group_a : List[str]
        Sample names of group A (numerator of the fold change).
    group_b : List[str]
        Sample names of group B (denominator).
    min_total : int, optional
        Loci with total counts below this across all used samples are
        reported with NaN statistics (kept in the table, never
        silently dropped). Default is ``10``.

    Returns
    -------
    pd.DataFrame
        Columns: ``transcript_id``, ``value`` (log2FC of mean CPM,
        ``log2((meanA+1)/(meanB+1))``), ``padj``.

    Raises
    ------
    ValueError
        If any requested sample is missing from ``counts`` or a group
        is empty.
    """
    try:
        from scipy.stats import mannwhitneyu
    except ImportError as exc:
        raise SystemExit(
            f"scipy is required for Mann-Whitney p values: {exc}. "
            "Install it into the current environment (no auto-install)."
        )
    for name in list(group_a) + list(group_b):
        if name not in counts.columns:
            raise ValueError(f"Sample {name} not in counts table: {list(counts.columns)}")
    if not group_a or not group_b:
        raise ValueError("Both groups must be non-empty")

    used = counts[list(group_a) + list(group_b)]
    lib_sizes = used.sum(axis=0)
    if (lib_sizes <= 0).any():
        raise ValueError(f"Zero library size for samples: {list(lib_sizes[lib_sizes <= 0].index)}")
    cpm = used.div(lib_sizes, axis=1) * 1e6

    mean_a = cpm[list(group_a)].mean(axis=1)
    mean_b = cpm[list(group_b)].mean(axis=1)
    value = np.log2((mean_a + 1.0) / (mean_b + 1.0))

    total = used.sum(axis=1)
    p_raw = np.full(len(used), np.nan, dtype=float)
    testable = total >= min_total
    logger.info(f"Testing {int(testable.sum())}/{len(used)} loci (>= {min_total} total counts)")
    for idx in np.flatnonzero(testable.to_numpy()):
        _, p_raw[idx] = mannwhitneyu(
            cpm.iloc[idx][list(group_a)].to_numpy(),
            cpm.iloc[idx][list(group_b)].to_numpy(),
            alternative="two-sided",
        )
    padj = bh_adjust(p_raw)

    out = pd.DataFrame(
        {
            "transcript_id": used.index.astype(str),
            "value": value.to_numpy(),
            "padj": padj,
        }
    )
    n_sig = int((out["padj"] < 0.05).sum())
    logger.info(f"Computed {len(out)} loci; {n_sig} with padj < 0.05")
    return out


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser.
    """
    parser = argparse.ArgumentParser(
        description="featureCounts locus table -> te_manhattan.py values (log2FC + padj).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--counts", required=True, help="featureCounts output table (-g transcript_id).")
    parser.add_argument("--group-a", required=True,
                        help="Comma-separated sample names (numerator, e.g. treated).")
    parser.add_argument("--group-b", required=True,
                        help="Comma-separated sample names (denominator, e.g. control).")
    parser.add_argument("--min-total", type=int, default=10,
                        help="Min total counts across used samples to test a locus.")
    parser.add_argument("-o", "--out", required=True, help="Output values TSV.")
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
    counts = load_featurecounts(args.counts)
    group_a = [s.strip() for s in args.group_a.split(",") if s.strip()]
    group_b = [s.strip() for s in args.group_b.split(",") if s.strip()]
    values = compute_values(counts, group_a, group_b, min_total=args.min_total)
    values.to_csv(args.out, sep="\t", index=False)
    logger.info(f"Saved {args.out}")


if __name__ == "__main__":
    main()
