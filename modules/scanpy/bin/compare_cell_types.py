#!/usr/bin/env python3
"""Compare cell type proportions between samples.

Each --group specifies an h5ad file followed by sample pairs.
Output subdirectories are named after the h5ad filename stem.

Outputs per group:
  - {prefix}_{id1}_vs_{id2}.tsv        proportion table
  - {prefix}_{id1}_vs_{id2}.png         proportion bar chart
  - {prefix}_bar_{id1}_vs_{id2}.png     DA bar comparison (proportions + FC)
  - {prefix}_volcano_{id1}_vs_{id2}.png volcano plot
  - {prefix}_heatmap.png                proportion change heatmap
  - {prefix}_heatmap_fc.png             log2FC heatmap
  - {prefix}_da_{id1}_vs_{id2}.tsv      DA results table
"""
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Proportion analysis
# ---------------------------------------------------------------------------
def get_cell_type_proportions(adata: ad.AnnData, sample1: str, sample2: str) -> Tuple[Optional[pd.Series], Optional[pd.Series]]:
    """Get cell type proportions for two samples."""
    mask1 = adata.obs["sample_id"] == sample1
    mask2 = adata.obs["sample_id"] == sample2

    if mask1.sum() == 0 or mask2.sum() == 0:
        log.warning("Sample not found: %s (%d) or %s (%d)", sample1, mask1.sum(), sample2, mask2.sum())
        return None, None

    ct1 = adata.obs[mask1]["cell_type"].value_counts(normalize=True) * 100
    ct2 = adata.obs[mask2]["cell_type"].value_counts(normalize=True) * 100

    return ct1, ct2


def plot_comparison(ct1: pd.Series, ct2: pd.Series, sample1: str, sample2: str, title: str, output_path: Path) -> None:
    """Plot cell type proportion comparison."""
    all_cts = sorted(set(ct1.index) | set(ct2.index))

    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(all_cts))
    width = 0.35

    vals1 = [ct1.get(ct, 0) for ct in all_cts]
    vals2 = [ct2.get(ct, 0) for ct in all_cts]

    bars1 = ax.bar(x - width/2, vals1, width, label=sample1, color="#4DBBD5", alpha=0.8)
    bars2 = ax.bar(x + width/2, vals2, width, label=sample2, color="#E64B35", alpha=0.8)

    ax.set_xlabel("Cell Type", fontsize=10)
    ax.set_ylabel("Proportion (%)", fontsize=10)
    ax.set_title(title, fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(all_cts, rotation=45, ha="right", fontsize=8)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    for bar in bars1:
        height = bar.get_height()
        if height > 1:
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}%', ha='center', va='bottom', fontsize=7)
    for bar in bars2:
        height = bar.get_height()
        if height > 1:
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}%', ha='center', va='bottom', fontsize=7)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", output_path)


def plot_proportion_heatmap(comparisons_data: List[Tuple[str, pd.Series, pd.Series]], title: str, output_path: Path) -> None:
    """Plot heatmap of cell type proportion changes (percentage point difference)."""
    all_cts: set = set()
    for comp_name, ct1, ct2 in comparisons_data:
        all_cts.update(ct1.index)
        all_cts.update(ct2.index)
    all_cts_list = sorted(all_cts)

    changes = []
    comp_names = []
    for comp_name, ct1, ct2 in comparisons_data:
        changes.append([ct2.get(ct, 0) - ct1.get(ct, 0) for ct in all_cts_list])
        comp_names.append(comp_name)

    changes_arr = np.array(changes)

    fig, ax = plt.subplots(figsize=(14, 5))
    im = ax.imshow(changes_arr, cmap="RdBu_r", aspect="auto", vmin=-30, vmax=30)

    ax.set_xticks(np.arange(len(all_cts_list)))
    ax.set_yticks(np.arange(len(comp_names)))
    ax.set_xticklabels(all_cts_list, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(comp_names, fontsize=9)

    for i in range(len(comp_names)):
        for j in range(len(all_cts_list)):
            val = changes_arr[i, j]
            if abs(val) > 0.5:
                ax.text(j, i, f'{val:+.1f}', ha='center', va='center', fontsize=7,
                       color='white' if abs(val) > 15 else 'black')

    ax.set_title(title, fontsize=12)
    plt.colorbar(im, ax=ax, label="Change (%)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", output_path)


def save_comparison_table(ct1: pd.Series, ct2: pd.Series, sample1: str, sample2: str, output_path: Path) -> None:
    """Save comparison table as TSV."""
    all_cts = sorted(set(ct1.index) | set(ct2.index))

    rows = []
    for ct in all_cts:
        pct1 = ct1.get(ct, 0)
        pct2 = ct2.get(ct, 0)
        change = pct2 - pct1
        rows.append({
            "cell_type": ct,
            f"pct_{sample1}": round(pct1, 2),
            f"pct_{sample2}": round(pct2, 2),
            "change": round(change, 2),
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep="\t", index=False)
    log.info("Saved: %s", output_path)


# ---------------------------------------------------------------------------
# DA analysis (Fisher's exact test)
# ---------------------------------------------------------------------------
def perform_da_analysis(adata: ad.AnnData, sample1: str, sample2: str) -> Tuple[pd.DataFrame, float, float]:
    """Perform differential abundance analysis via Fisher's exact test."""
    samples = [sample1, sample2]
    cell_types = sorted(adata.obs["cell_type"].unique())

    contingency = []
    valid_cts = []
    for ct in cell_types:
        row = []
        for sample in samples:
            n = ((adata.obs["cell_type"] == ct) & (adata.obs["sample_id"] == sample)).sum()
            row.append(n)
        if sum(row) > 0:
            contingency.append(row)
            valid_cts.append(ct)

    if len(contingency) < 2:
        log.warning("Not enough cell types for DA analysis")
        return pd.DataFrame(), 0, 1.0

    contingency = np.array(contingency, dtype=float)
    contingency_pseudo = contingency + 0.5  # pseudocount for chi2

    try:
        chi2, p_chi2, _, _ = stats.chi2_contingency(contingency_pseudo)
    except ValueError:
        chi2, p_chi2 = 0, 1.0

    results = []
    for i, ct in enumerate(valid_cts):
        n1 = contingency[i, 0]
        n2 = contingency[i, 1]
        total1 = contingency[:, 0].sum()
        total2 = contingency[:, 1].sum()

        p1 = n1 / total1 if total1 > 0 else 0
        p2 = n2 / total2 if total2 > 0 else 0

        if p1 > 0 and p2 > 0:
            fc = p2 / p1
            log2fc = np.log2(fc)
        elif p1 == 0 and p2 > 0:
            fc = float("inf")
            log2fc = float("inf")
        else:
            fc = 0.0
            log2fc = float("-inf")

        table = np.array([[n1, total1 - n1], [n2, total2 - n2]])
        _, pval = stats.fisher_exact(table)

        results.append({
            "cell_type": ct,
            "n_sample1": int(n1),
            "n_sample2": int(n2),
            "pct_sample1": round(p1 * 100, 2),
            "pct_sample2": round(p2 * 100, 2),
            "fc": round(fc, 4) if np.isfinite(fc) else fc,
            "log2fc": round(log2fc, 4) if np.isfinite(log2fc) else log2fc,
            "pvalue": pval,
            "neg_log10p": round(-np.log10(pval), 4) if pval > 0 else 300,
        })

    return pd.DataFrame(results), chi2, p_chi2


def plot_volcano(df: pd.DataFrame, sample1: str, sample2: str, title: str, output_path: Path) -> None:
    """Plot volcano plot for DA analysis."""
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = []
    for _, row in df.iterrows():
        if row["pvalue"] < 0.001 and abs(row["log2fc"]) > 1:
            colors.append("#E64B35")
        elif row["pvalue"] < 0.05 and abs(row["log2fc"]) > 0.5:
            colors.append("#F39B7F")
        else:
            colors.append("#8491B4")

    ax.scatter(df["log2fc"], df["neg_log10p"], c=colors, s=100, alpha=0.7, edgecolors="black", linewidth=0.5)

    for _, row in df.iterrows():
        if row["pvalue"] < 0.001 and abs(row["log2fc"]) > 0.5:
            ax.annotate(row["cell_type"], (row["log2fc"], row["neg_log10p"]),
                       fontsize=8, ha="center", va="bottom", fontweight="bold")

    ax.axhline(-np.log10(0.001), color="gray", linestyle="--", alpha=0.5, label="p=0.001")
    ax.axvline(-1, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(1, color="gray", linestyle="--", alpha=0.5)

    ax.set_xlabel("log2(Fold Change)", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", output_path)


def plot_da_bar(df: pd.DataFrame, sample1: str, sample2: str, title: str, output_path: Path) -> None:
    """Plot DA bar comparison: proportions + fold change side by side."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    df_sorted = df.sort_values("log2fc", ascending=True)
    x = np.arange(len(df_sorted))
    width = 0.35

    axes[0].barh(x - width/2, df_sorted["pct_sample1"], width, label=sample1, color="#4DBBD5", alpha=0.8)
    axes[0].barh(x + width/2, df_sorted["pct_sample2"], width, label=sample2, color="#E64B35", alpha=0.8)
    axes[0].set_yticks(x)
    axes[0].set_yticklabels(df_sorted["cell_type"], fontsize=9)
    axes[0].set_xlabel("Proportion (%)", fontsize=10)
    axes[0].set_title("Cell Type Proportions", fontsize=12)
    axes[0].legend(fontsize=9)
    axes[0].grid(axis="x", alpha=0.3)

    colors = ["#E64B35" if p < 0.001 else "#F39B7F" if p < 0.05 else "#8491B4"
              for p in df_sorted["pvalue"]]
    axes[1].barh(x, df_sorted["log2fc"], color=colors, alpha=0.8, edgecolor="black", linewidth=0.5)
    axes[1].set_yticks(x)
    axes[1].set_yticklabels(df_sorted["cell_type"], fontsize=9)
    axes[1].set_xlabel("log2(Fold Change)", fontsize=10)
    axes[1].set_title("Fold Change (Sample2 vs Sample1)", fontsize=12)
    axes[1].axvline(0, color="black", linewidth=0.5)
    axes[1].axvline(-1, color="gray", linestyle="--", alpha=0.5)
    axes[1].axvline(1, color="gray", linestyle="--", alpha=0.5)
    axes[1].grid(axis="x", alpha=0.3)

    plt.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", output_path)


def plot_fc_heatmap(da_results: List[Tuple[str, pd.DataFrame]], title: str, output_path: Path) -> None:
    """Plot heatmap of log2FC across comparisons."""
    all_cts: set = set()
    for comp_name, df in da_results:
        all_cts.update(df["cell_type"].tolist())
    all_cts_list = sorted(all_cts)

    matrix = []
    comp_names = []
    for comp_name, df in da_results:
        row = []
        for ct in all_cts_list:
            vals = df[df["cell_type"] == ct]["log2fc"].values
            row.append(float(vals[0]) if len(vals) > 0 else 0)
        matrix.append(row)
        comp_names.append(comp_name)

    matrix_arr = np.array(matrix)

    fig, ax = plt.subplots(figsize=(14, 5))
    im = ax.imshow(matrix_arr, cmap="RdBu_r", aspect="auto", vmin=-5, vmax=5)

    ax.set_xticks(np.arange(len(all_cts_list)))
    ax.set_yticks(np.arange(len(comp_names)))
    ax.set_xticklabels(all_cts_list, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(comp_names, fontsize=9)

    for i in range(len(comp_names)):
        for j in range(len(all_cts_list)):
            val = matrix_arr[i, j]
            if abs(val) > 0.1:
                ax.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=7,
                       color='white' if abs(val) > 2.5 else 'black')

    ax.set_title(title, fontsize=12)
    plt.colorbar(im, ax=ax, label="log2(FC)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", output_path)


# ---------------------------------------------------------------------------
# Process one h5ad group
# ---------------------------------------------------------------------------
def _parse_h5ad_prefix(stem: str) -> str:
    """Extract tissue_method prefix from h5ad filename stem.

    Examples:
        Uterus_cellranger_v3_annotated -> Uterus_cellranger
        ovaries_scTE_v3_annotated      -> ovaries_scTE
    """
    parts = stem.split("_")
    for i, p in enumerate(parts):
        if p in ("cellranger", "scTE"):
            return "_".join(parts[:i + 1])
    return "_".join(parts[:-2]) if len(parts) > 2 else stem


def process_group(h5ad_path: Path, comparisons: List[Tuple[str, str]], output_dir: Path) -> None:
    """Run all comparisons for one h5ad file, output to a named subdirectory."""
    stem = h5ad_path.stem
    group_dir = output_dir / stem
    group_dir.mkdir(parents=True, exist_ok=True)

    prefix = _parse_h5ad_prefix(stem)

    log.info("=" * 60)
    log.info("Group: %s (%d comparisons)", stem, len(comparisons))
    log.info("=" * 60)

    adata = ad.read_h5ad(str(h5ad_path))
    log.info("Loaded %s: %d cells", stem, adata.n_obs)

    # Collect for combined heatmaps
    proportion_comparisons: List[Tuple[str, pd.Series, pd.Series]] = []
    da_all_results: List[Tuple[str, pd.DataFrame]] = []

    for sample1, sample2 in comparisons:
        id1 = sample1.split("-")[1]
        id2 = sample2.split("-")[1]
        comp_name = f"{id1} vs {id2}"
        fname = f"{prefix}_{id1}_vs_{id2}"

        # --- Proportion comparison ---
        log.info("Comparing %s vs %s", sample1, sample2)
        ct1, ct2 = get_cell_type_proportions(adata, sample1, sample2)
        if ct1 is not None:
            proportion_comparisons.append((comp_name, ct1, ct2))
            save_comparison_table(ct1, ct2, sample1, sample2,
                                 group_dir / f"{fname}.tsv")
            plot_comparison(ct1, ct2, sample1, sample2,
                           f"{prefix}: {sample1} vs {sample2}",
                           group_dir / f"{fname}.png")

        # --- DA analysis ---
        log.info("DA analysis: %s vs %s", sample1, sample2)
        da_df, chi2, p_chi2 = perform_da_analysis(adata, sample1, sample2)
        if not da_df.empty:
            da_all_results.append((comp_name, da_df))
            da_df.to_csv(group_dir / f"{prefix}_da_{id1}_vs_{id2}.tsv", sep="\t", index=False)
            log.info("Chi2=%.2f, p=%.2e for %s", chi2, p_chi2, comp_name)

            plot_volcano(da_df, sample1, sample2,
                        f"Volcano: {sample1} vs {sample2}",
                        group_dir / f"{prefix}_volcano_{id1}_vs_{id2}.png")
            plot_da_bar(da_df, sample1, sample2,
                       f"DA: {sample1} vs {sample2}",
                       group_dir / f"{prefix}_bar_{id1}_vs_{id2}.png")

    # --- Combined heatmaps ---
    if proportion_comparisons:
        plot_proportion_heatmap(proportion_comparisons,
                               f"{prefix} Cell Type Proportion Changes",
                               group_dir / f"{prefix}_heatmap.png")

    if da_all_results:
        plot_fc_heatmap(da_all_results,
                        f"{prefix} log2(Fold Change) Heatmap",
                        group_dir / f"{prefix}_heatmap_fc.png")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare cell type proportions between samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Example (4 groups: 2 tissues x 2 quantification methods):

  python compare_cell_types.py \\
    --output-dir /data/results \\
    --group /path/to/Uterus_cellranger.h5ad  zigong-11238-10XSC3 zigong-2200671-10XSC3 \\
                                           zigong-11238-10XSC3 zigong-1411200-10XSC3 \\
    --group /path/to/Uterus_scTE.h5ad       zigong-11238-10XSC3 zigong-2200671-10XSC3 \\
    --group /path/to/ovaries_cellranger.h5ad luanchao-21310-10XSC3 luanchao-11238-10XSC3 \\
    --group /path/to/ovaries_scTE.h5ad       luanchao-21310-10XSC3 luanchao-11238-10XSC3
""")

    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for plots and tables")

    parser.add_argument("--group", nargs="+", action="append", required=True,
                        metavar="ARG",
                        help="H5AD_PATH SAMPLE1 SAMPLE2 [SAMPLE3 SAMPLE4 ...]. "
                             "Repeat for multiple groups. Output subdirectory named after h5ad stem.")

    args = parser.parse_args()

    for i, g in enumerate(args.group):
        if len(g) < 3:
            parser.error(f"--group #{i+1} needs at least H5AD + 2 SAMPLE names, got {len(g)}: {g}")
        if (len(g) - 1) % 2 != 0:
            parser.error(f"--group #{i+1}: sample names must be in pairs (even count), got {len(g)-1} samples after H5AD")

    return args


def main():
    args = parse_args()

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    for group_args in args.group:
        h5ad_path = Path(group_args[0])
        sample_names = group_args[1:]
        comparisons: List[Tuple[str, str]] = [
            (sample_names[j], sample_names[j + 1]) for j in range(0, len(sample_names), 2)
        ]
        process_group(h5ad_path, comparisons, output_dir)

    log.info("=" * 60)
    log.info("All comparisons done!")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
