#!/usr/bin/env python3
"""Compare HLA Class I vs Class II read coverage with stacked bar chart."""

import argparse
import csv
import math
import os

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLASS_I_GENES = ["A", "B", "C", "E", "F", "G", "H", "J", "K", "L", "V"]
CLASS_II_GENES = [
    "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB6", "DRB7", "DRB8", "DRB9",
    "DQA1", "DQB1", "DPA1", "DPB1", "DMA", "DMB", "DOA", "DOB", "DRA",
]


def load_gene_means(tsv_path):
    """Load per-gene mean read counts from the coverage TSV.

    Parameters
    ----------
    tsv_path : str
        Path to the TSV produced by summarize_read_coverage.py.

    Returns
    -------
    dict[str, float]
        Mapping of gene name to mean total_reads across all alleles.
    """
    sums = {}
    counts = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["allele"] == "Not_typed":
                continue
            gene = row["gene"]
            n = int(row["total_reads"])
            sums[gene] = sums.get(gene, 0) + n
            counts[gene] = counts.get(gene, 0) + 1
    return {g: sums[g] / counts[g] for g in sums}


def auto_ncol(n_items, fig_height):
    """Determine the number of legend columns based on item count and figure height.

    Each legend row occupies approximately 0.3 inches.  The available height
    for the legend is roughly 80% of the figure height.  Columns are added
    until all items fit within the available height, capped at 6.

    Parameters
    ----------
    n_items : int
        Total number of legend entries.
    fig_height : float
        Figure height in inches.

    Returns
    -------
    int
        Optimal number of columns for the legend.
    """
    row_height = 0.3
    available = fig_height * 0.8
    max_rows = max(1, int(available / row_height))
    ncol = min(6, max(1, math.ceil(n_items / max_rows)))
    return ncol


def plot_stacked_coverage(gene_means, c1_mean_sum, c2_mean_sum, output_path):
    """Plot stacked bar chart comparing Class I vs Class II read coverage.

    Parameters
    ----------
    gene_means : dict[str, float]
        Mapping of gene name to mean total_reads.
    c1_mean_sum : float
        Sum of per-gene mean reads for Class I (= bar height).
    c2_mean_sum : float
        Sum of per-gene mean reads for Class II (= bar height).
    output_path : str
        Path to save the output figure.
    """
    c1_genes = [g for g in CLASS_I_GENES if g in gene_means]
    c2_genes = [g for g in CLASS_II_GENES if g in gene_means]
    c1_vals = [gene_means[g] for g in c1_genes]
    c2_vals = [gene_means[g] for g in c2_genes]

    cmap1 = plt.cm.Blues(np.linspace(0.35, 0.85, len(c1_genes)))
    cmap2 = plt.cm.Oranges(np.linspace(0.35, 0.85, len(c2_genes)))

    fig_height = 6
    fig, ax = plt.subplots(figsize=(10, fig_height))
    x = np.array([0, 1])
    bar_w = 0.5

    # Class I stack
    bottom1 = 0
    for i, (g, v) in enumerate(zip(c1_genes, c1_vals)):
        ax.bar(x[0], v, bar_w, bottom=bottom1, color=cmap1[i],
               edgecolor="white", linewidth=0.5, label=g)
        bottom1 += v

    # Class II stack
    bottom2 = 0
    for i, (g, v) in enumerate(zip(c2_genes, c2_vals)):
        ax.bar(x[1], v, bar_w, bottom=bottom2, color=cmap2[i],
               edgecolor="white", linewidth=0.5, label=g)
        bottom2 += v

    # Total labels
    ax.text(x[0], bottom1, f"total={c1_mean_sum:,.0f}", ha="center", va="bottom", fontsize=10)
    ax.text(x[1], bottom2, f"total={c2_mean_sum:,.0f}", ha="center", va="bottom", fontsize=10)

    # Axes
    ax.set_xticks(x)
    ax.set_xticklabels(["Class I", "Class II"], fontsize=13, fontweight="bold")
    ax.set_ylabel("Mean total reads per allele", fontsize=12)
    ax.set_title("HLA Read Coverage: Class I vs Class II", fontsize=13, fontweight="bold")
    ax.set_ylim(0, max(bottom1, bottom2) * 1.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    all_handles = ax.get_legend_handles_labels()[0]
    all_labels = [f"I: {g}" for g in c1_genes] + [f"II: {g}" for g in c2_genes]
    ncol = auto_ncol(len(all_labels), fig_height)

    ax.legend(all_handles, all_labels, ncol=ncol, fontsize=8, frameon=False,
              loc="center left", bbox_to_anchor=(1.02, 0.45), borderaxespad=0.0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved to {output_path}")


def plot_per_gene_mean(gene_means, output_path):
    """Plot per-gene mean reads, colored by HLA class.

    Parameters
    ----------
    gene_means : dict[str, float]
        Mapping of gene name to mean total_reads.
    output_path : str
        Path to save the output figure.
    """
    c1_genes = [g for g in CLASS_I_GENES if g in gene_means]
    c2_genes = [g for g in CLASS_II_GENES if g in gene_means]
    all_genes = c1_genes + c2_genes
    all_vals = [gene_means[g] for g in all_genes]
    colors = [plt.cm.Blues(0.6)] * len(c1_genes) + [plt.cm.Oranges(0.6)] * len(c2_genes)

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(all_genes))
    ax.bar(x, all_vals, color=colors, edgecolor="white", linewidth=0.5)

    # Axes
    ax.set_xticks(x)
    ax.set_xticklabels(all_genes, fontsize=9, rotation=45, ha="right")
    ax.set_ylabel("Mean reads per allele", fontsize=12)
    ax.set_title("Per-Gene Mean Read Coverage", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Class divider and labels
    y_max = max(all_vals)
    ax.axvline(x=len(c1_genes) - 0.5, color="gray", linestyle="--", linewidth=0.8)
    ax.text(len(c1_genes) / 2 - 0.5, -y_max * 0.18, "Class I",
            ha="center", va="top", fontsize=11, fontweight="bold", color=plt.cm.Blues(0.6))
    ax.text(len(c1_genes) + len(c2_genes) / 2 - 0.5, -y_max * 0.18, "Class II",
            ha="center", va="top", fontsize=11, fontweight="bold", color=plt.cm.Oranges(0.6))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved to {output_path}")


def plot_mean_per_gene(gene_means, output_path):
    """Plot average read coverage per gene for Class I vs Class II.

    Each bar is the mean of per-gene mean reads within a class.
    Statistical test (Mann-Whitney U) is performed on per-gene means.

    Parameters
    ----------
    gene_means : dict[str, float]
        Mapping of gene name to mean total_reads.
    output_path : str
        Path to save the output figure.
    """
    c1_vals = [gene_means[g] for g in CLASS_I_GENES if g in gene_means]
    c2_vals = [gene_means[g] for g in CLASS_II_GENES if g in gene_means]

    c1_arr = np.array(c1_vals)
    c2_arr = np.array(c2_vals)

    c1_mean = c1_arr.mean()
    c2_mean = c2_arr.mean()
    c1_sem = c1_arr.std(ddof=1) / np.sqrt(len(c1_arr))
    c2_sem = c2_arr.std(ddof=1) / np.sqrt(len(c2_arr))

    u_stat, p_val = stats.mannwhitneyu(c1_arr, c2_arr, alternative="two-sided")
    if p_val < 0.001:
        sig_label = "***"
    elif p_val < 0.01:
        sig_label = "**"
    elif p_val < 0.05:
        sig_label = "*"
    else:
        sig_label = "n.s."

    fig, ax = plt.subplots(figsize=(5, 5))
    x = np.array([0, 1])
    bar_w = 0.4

    ax.bar(x, [c1_mean, c2_mean], bar_w, yerr=[c1_sem, c2_sem],
           color=[plt.cm.Blues(0.6), plt.cm.Oranges(0.6)],
           edgecolor="black", linewidth=0.8, capsize=6, error_kw={"linewidth": 1.2})

    # Significance bracket
    y_max = max(c1_mean + c1_sem, c2_mean + c2_sem)
    h = y_max * 0.06
    bracket_y = y_max + h * 1.5
    ax.plot([0, 0, 1, 1], [bracket_y, bracket_y + h, bracket_y + h, bracket_y],
            color="black", linewidth=1.2)
    ax.text(0.5, bracket_y + h, f"Mann-Whitney U test\n{sig_label} p={p_val:.2e}",
            ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Axes
    ax.set_xticks(x)
    ax.set_xticklabels(["Class I", "Class II"], fontsize=12, fontweight="bold")
    ax.set_ylabel("Mean reads per allele (per gene)", fontsize=12)
    ax.set_title("Average Read Coverage per Gene", fontsize=13, fontweight="bold")
    ax.set_ylim(0, bracket_y + h * 5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Class I vs Class II coverage plots.")
    parser.add_argument("tsv", help="Input TSV from summarize_read_coverage.py")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    gene_means = load_gene_means(args.tsv)

    c1_mean_sum = sum(v for g, v in gene_means.items() if g in CLASS_I_GENES)
    c2_mean_sum = sum(v for g, v in gene_means.items() if g in CLASS_II_GENES)

    plot_stacked_coverage(gene_means, c1_mean_sum, c2_mean_sum,
                          os.path.join(args.outdir, "class_I_vs_II_stacked.png"))

    plot_per_gene_mean(gene_means,
                       os.path.join(args.outdir, "per_gene_mean.png"))

    plot_mean_per_gene(gene_means,
                       os.path.join(args.outdir, "mean_per_gene.png"))


if __name__ == "__main__":
    main()
