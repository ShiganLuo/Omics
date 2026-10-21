#!/usr/bin/env python3
"""GSVA (Gene Set Variation Analysis) — pure Python implementation.

Computes single-sample Gene Set Enrichment Scores using a modified
Kolmogorov-Smirnov statistic, following Hänzelmann et al. (2013).

Input:
    - Normalized expression matrix (TSV, gene_name × samples)
    - GMT gene set file (tab-separated: name, description, gene1, gene2, ...)

Output:
    - gsva_scores.tsv: KS enrichment score matrix (geneset × samples;
      positive = gene set concentrated at the high-expression end)
    - gsva_mean_z.tsv: mean per-gene z-score matrix (geneset × samples;
      per-gene z on log2(TPM+1) across samples, averaged within the set)
    - gsva_heatmap.png: heatmap of the mean z-score matrix

Usage:
    python gsva.py -i expression_tpm.tsv -g geneset.gmt -o outdir [--title GSVA]
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List
import numpy as np
import pandas as pd
from scipy.stats import rankdata
try:
    from .LogUtil import setup_logger
except ImportError:
    from LogUtil import setup_logger
logger = setup_logger(__name__)

def parse_gmt(gmt_path: str) -> Dict[str, List[str]]:
    """Parse a GMT (Gene Matrix Transposed) file.

    Format: each line is a gene set:
        <set_name>\\t<description>\\t<gene1>\\t<gene2>\\t...

    Parameters
    ----------
    gmt_path : str
        Path to the GMT file.

    Returns
    -------
    dict
        Mapping from gene set name to list of gene symbols.
    """
    gene_sets: Dict[str, List[str]] = {}
    with open(gmt_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g for g in parts[2:] if g]
            if genes:
                gene_sets[name] = genes
    logger.info(f"Parsed {len(gene_sets)} gene sets from {gmt_path}")
    return gene_sets


def compute_gsva(
    expr: pd.DataFrame,
    gene_sets: Dict[str, List[str]],
    abs_rank: bool = True,
) -> pd.DataFrame:
    """Compute KS-based enrichment scores for each sample and gene set.

    Walks the within-sample expression ranking from LOW to HIGH expression.
    The score is the max deviation of (miss_frac - hit_frac); a positive
    score means the gene set is concentrated at the HIGH-expression end,
    matching the GSVA/GSEA sign convention.

    Parameters
    ----------
    expr : pd.DataFrame
        Expression matrix (genes × samples). Index must be gene symbols.
    gene_sets : dict
        Gene set name → list of gene symbols.
    abs_rank : bool
        If True, use absolute ranking (two-tailed KS).

    Returns
    -------
    pd.DataFrame
        Enrichment score matrix (gene_sets × samples).
    """
    gene_names = expr.index.tolist()
    n_genes = len(gene_names)
    n_samples = expr.shape[1]
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Rank expression values per sample
    ranks = expr.apply(lambda col: rankdata(col.values, method="average"), axis=0)
    ranks_np = ranks.values  # (n_genes, n_samples)

    results = {}
    for gs_name, gs_genes in gene_sets.items():
        hit_indices = [gene_to_idx[g] for g in gs_genes if g in gene_to_idx]
        if not hit_indices:
            logger.warning(f"Gene set '{gs_name}': no genes found in expression data")
            results[gs_name] = np.zeros(n_samples)
            continue

        hit_mask = np.zeros(n_genes, dtype=bool)
        hit_mask[hit_indices] = True

        scores = np.zeros(n_samples)
        for s in range(n_samples):
            # Sort genes by rank for this sample
            order = np.argsort(ranks_np[:, s])
            hit_ordered = hit_mask[order]

            # Walkman statistic: cumulative sum of hits minus misses
            hit_cum = np.cumsum(hit_ordered)
            miss_cum = np.cumsum(~hit_ordered)

            hit_frac = hit_cum / len(hit_indices)
            miss_frac = miss_cum / (n_genes - len(hit_indices))

            # Positive = gene set concentrated at the HIGH-expression end
            diff = miss_frac - hit_frac

            if abs_rank:
                scores[s] = diff[np.argmax(np.abs(diff))]
            else:
                # Positive-only enrichment at the high-expression end
                scores[s] = np.max(diff)

        results[gs_name] = scores

    score_df = pd.DataFrame(
        results, index=expr.columns.tolist()
    ).T
    score_df.index.name = "gene_set"
    return score_df

def compute_mean_z(
    expr: pd.DataFrame,
    gene_sets: Dict[str, List[str]],
    log2_transform: bool = True,
) -> pd.DataFrame:
    """Compute mean per-gene z-score for each sample and gene set.

    Each gene is standardized across samples (optionally after log2(x + 1)
    transformation), then averaged over the genes of the set. Unlike the
    KS-based score, this reflects absolute expression changes between
    samples instead of within-sample rank positions.

    Parameters
    ----------
    expr : pd.DataFrame
        Expression matrix (genes × samples). Index must be gene symbols.
    gene_sets : dict
        Gene set name → list of gene symbols.
    log2_transform : bool
        Apply log2(x + 1) before standardization (for TPM/FPKM input).
        Set False if the input is already log-scale.

    Returns
    -------
    pd.DataFrame
        Mean z-score matrix (gene_sets × samples). Per-gene z-scores sum
        to zero across samples, so each row has mean ~0.
    """
    x = np.log2(expr + 1.0) if log2_transform else expr.astype(float)
    gene_mean = x.mean(axis=1)
    gene_std = x.std(axis=1).replace(0.0, 1.0)
    z = x.sub(gene_mean, axis=0).div(gene_std, axis=0)

    results = {}
    for gs_name, gs_genes in gene_sets.items():
        hit = [g for g in gs_genes if g in z.index]
        if not hit:
            logger.warning(f"Gene set '{gs_name}': no genes found in expression data")
            results[gs_name] = np.zeros(expr.shape[1])
            continue
        results[gs_name] = z.loc[hit].mean(axis=0).values

    score_df = pd.DataFrame(
        results, index=expr.columns.tolist()
    ).T
    score_df.index.name = "gene_set"
    return score_df


def plot_heatmap(
    scores: pd.DataFrame,
    output_path: str,
    title: str = "GSVA Enrichment Scores",
    value_label: str = "Row z-score",
    row_scale: bool = True,
    label_wrap: int = 40,
    sample_wrap: int = 12,
) -> None:
    """Plot heatmap of enrichment scores with a symmetric color scale.

    Parameters
    ----------
    scores : pd.DataFrame
        Score matrix (genesets × samples).
    output_path : str
        Path to save the heatmap PNG.
    title : str
        Plot title.
    value_label : str
        Colorbar label describing the plotted values.
    row_scale : bool
        If True, standardize each row (gene set) across samples to zero
        mean and unit variance before plotting. Restores contrast when
        the raw scores span a narrow range (e.g. averaged gene z-scores).
    label_wrap : int
        Max characters per line for wrapped gene-set (y) labels.
    sample_wrap : int
        Max characters per line for wrapped sample (x) labels.

    Notes
    -----
    Uses constrained_layout so the colorbar stays attached to the axes,
    and sizes the figure width from the longest wrapped label so long
    gene-set names do not squeeze the heatmap into a sliver.
    """
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_gs, n_samples = scores.shape

    values = scores
    if row_scale:
        row_mean = scores.mean(axis=1)
        row_std = scores.std(axis=1).replace(0.0, 1.0)
        values = scores.sub(row_mean, axis=0).div(row_std, axis=0)

    y_labels = [
        "\n".join(textwrap.wrap(str(name), label_wrap)) or str(name)
        for name in scores.index
    ]
    x_labels = [
        "\n".join(textwrap.wrap(str(name), sample_wrap)) or str(name)
        for name in scores.columns
    ]
    max_label_chars = max(
        max(len(line) for line in lbl.split("\n")) for lbl in y_labels
    )

    fig_h = max(4, n_gs * 0.35 + 1.5)
    fig_w = max(6, n_samples * 0.6 + 2 + max_label_chars * 0.09)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    vmax = float(np.nanmax(np.abs(values.values))) or 1.0
    im = ax.imshow(
        values.values,
        cmap="RdBu_r",
        aspect="auto",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )

    ax.set_xticks(range(n_samples))
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_gs))
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_title(title, fontsize=11)

    # Thin white cell separators
    ax.set_xticks(np.arange(-0.5, n_samples), minor=True)
    ax.set_yticks(np.arange(-0.5, n_gs), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.colorbar(im, ax=ax, shrink=0.8, label=value_label)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    logger.info(f"Heatmap saved to {output_path}")


def main() -> None:
    p = argparse.ArgumentParser(
        description="GSVA: Gene Set Variation Analysis (Python implementation)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("-i", "--input", required=True,
                   help="Normalized expression matrix (TSV, gene_name × samples)")
    p.add_argument("-g", "--gmt", required=True,
                   help="GMT gene set file")
    p.add_argument("-o", "--outdir", required=True,
                   help="Output directory")
    p.add_argument("--title", default="GSVA Enrichment Scores",
                   help="Plot title (default: 'GSVA Enrichment Scores')")
    p.add_argument("--no-log2", action="store_true",
                   help="Disable log2(x + 1) for mean z-score (input already log-scale)")
    p.add_argument("--no-abs-rank", action="store_true",
                   help="Disable absolute ranking (positive-only enrichment)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load expression matrix
    expr = pd.read_csv(args.input, sep="\t", index_col=0)
    logger.info(f"Loaded expression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")

    # Parse GMT
    gene_sets = parse_gmt(args.gmt)

    # Compute KS-based scores (sign convention: positive = high-expression end)
    scores = compute_gsva(
        expr,
        gene_sets,
        abs_rank=not args.no_abs_rank,
    )

    # Save KS scores
    scores_path = os.path.join(args.outdir, "gsva_scores.tsv")
    scores.to_csv(scores_path, sep="\t")
    logger.info(f"GSVA KS scores saved to {scores_path}")

    # Compute mean per-gene z-score matrix (used for the heatmap)
    mean_z = compute_mean_z(
        expr,
        gene_sets,
        log2_transform=not args.no_log2,
    )
    mean_z_path = os.path.join(args.outdir, "gsva_mean_z.tsv")
    mean_z.to_csv(mean_z_path, sep="\t")
    logger.info(f"Mean z-score matrix saved to {mean_z_path}")

    # Plot heatmap of row-scaled mean z-scores
    heatmap_path = os.path.join(args.outdir, "gsva_heatmap.png")
    transform_note = "" if args.no_log2 else " (log2(x+1) of TPM before z)"
    value_label = f"Row z-score of mean gene z{transform_note}"
    plot_heatmap(mean_z, heatmap_path, title=args.title, value_label=value_label)

    logger.info("GSVA analysis completed successfully")


if __name__ == "__main__":
    main()
