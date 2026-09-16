#!/usr/bin/env python3
"""GSVA (Gene Set Variation Analysis) — pure Python implementation.

Computes single-sample Gene Set Enrichment Scores using a modified
Kolmogorov-Smirnov statistic, following Hänzelmann et al. (2013).

Input:
    - Normalized expression matrix (TSV, gene_name × samples)
    - GMT gene set file (tab-separated: name, description, gene1, gene2, ...)

Output:
    - GSVA enrichment score matrix (TSV, geneset × samples)
    - Clustered heatmap (PNG)

Usage:
    python gsva.py -i expression_tpm.tsv -g geneset.gmt -o outdir [--title GSVA]
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import rankdata

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


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
    kcdf: str = "Gaussian",
    abs_rank: bool = True,
) -> pd.DataFrame:
    """Compute GSVA enrichment scores for each sample and gene set.

    Parameters
    ----------
    expr : pd.DataFrame
        Expression matrix (genes × samples). Index must be gene symbols.
    gene_sets : dict
        Gene set name → list of gene symbols.
    kcdf : str
        Kernel CDF type: "Gaussian" for continuous (TPM/FPKM),
        "Poisson" for count data.
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

            if abs_rank:
                diff = hit_frac - miss_frac
                scores[s] = diff[np.argmax(np.abs(diff))]
            else:
                # Positive-only enrichment
                diff = hit_frac - miss_frac
                scores[s] = np.max(diff)

        results[gs_name] = scores

    score_df = pd.DataFrame(
        results, index=expr.columns.tolist()
    ).T
    score_df.index.name = "gene_set"
    return score_df


def plot_heatmap(
    scores: pd.DataFrame,
    output_path: str,
    title: str = "GSVA Enrichment Scores",
) -> None:
    """Plot clustered heatmap of GSVA scores.

    Parameters
    ----------
    scores : pd.DataFrame
        GSVA score matrix (genesets × samples).
    output_path : str
        Path to save the heatmap PNG.
    title : str
        Plot title.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_gs, n_samples = scores.shape
    fig_h = max(4, n_gs * 0.35 + 1.5)
    fig_w = max(6, n_samples * 0.6 + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    vmax = max(abs(scores.values.min()), abs(scores.values.max()))
    im = ax.imshow(
        scores.values,
        cmap="RdBu_r",
        aspect="auto",
        vmin=-vmax,
        vmax=vmax,
    )

    ax.set_xticks(range(n_samples))
    ax.set_xticklabels(scores.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_gs))
    ax.set_yticklabels(scores.index, fontsize=8)
    ax.set_title(title, fontsize=11)

    plt.colorbar(im, ax=ax, shrink=0.8, label="Enrichment Score")
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
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
    p.add_argument("--kcdf", choices=["Gaussian", "Poisson"], default="Gaussian",
                   help="Kernel CDF: Gaussian for TPM/FPKM, Poisson for counts")
    p.add_argument("--no-abs-rank", action="store_true",
                   help="Disable absolute ranking (positive-only enrichment)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load expression matrix
    expr = pd.read_csv(args.input, sep="\t", index_col=0)
    logger.info(f"Loaded expression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")

    # Parse GMT
    gene_sets = parse_gmt(args.gmt)

    # Compute GSVA
    scores = compute_gsva(
        expr,
        gene_sets,
        kcdf=args.kcdf,
        abs_rank=not args.no_abs_rank,
    )

    # Save scores
    scores_path = os.path.join(args.outdir, "gsva_scores.tsv")
    scores.to_csv(scores_path, sep="\t")
    logger.info(f"GSVA scores saved to {scores_path}")

    # Plot heatmap
    heatmap_path = os.path.join(args.outdir, "gsva_heatmap.png")
    plot_heatmap(scores, heatmap_path, title=args.title)

    logger.info("GSVA analysis completed successfully")


if __name__ == "__main__":
    main()
