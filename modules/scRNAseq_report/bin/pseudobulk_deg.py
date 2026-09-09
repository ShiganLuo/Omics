#!/usr/bin/env python3
"""Pseudobulk DEG analysis per cell type with pseudo-replicates (PyDESeq2).

For each cell type, aggregate raw counts by sample, split each sample
into 2 pseudo-replicates, run PyDESeq2 for pairwise comparisons.

Output: per-cell-type DEG TSV + volcano plots (TE highlighted).
"""

import argparse
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")



def pseudobulk_with_replicates(
    adata_ct: ad.AnnData,
    sample_key: str = "sample_id",
    n_splits: int = 2,
) -> tuple:
    """Aggregate raw counts by sample, split into pseudo-replicates.

    Returns:
        counts_df: DataFrame (pseudo-replicates × genes)
        meta_df: DataFrame with 'condition' and 'replicate' columns
    """
    counts = adata_ct.layers["counts"]
    if hasattr(counts, "toarray"):
        counts = counts.toarray()

    gene_names = adata_ct.var_names.tolist()
    sample_ids = adata_ct.obs[sample_key].values

    all_rows = []
    meta_rows = []

    for sample in sorted(set(sample_ids)):
        mask = sample_ids == sample
        cell_indices = np.where(mask)[0]
        np.random.shuffle(cell_indices)
        splits = np.array_split(cell_indices, n_splits)

        for i, split_idx in enumerate(splits):
            sub_counts = counts[split_idx]
            row_sum = np.asarray(sub_counts.sum(axis=0)).flatten()
            rep_name = f"{sample}_rep{i}"
            all_rows.append(row_sum)
            meta_rows.append({
                "sample": rep_name,
                "condition": sample,
                "replicate": f"rep{i}",
            })

    counts_df = pd.DataFrame(all_rows, columns=gene_names)
    counts_df.index = [m["sample"] for m in meta_rows]
    meta_df = pd.DataFrame(meta_rows).set_index("sample")

    return counts_df, meta_df


def run_deg(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample1: str,
    sample2: str,
    cell_type: str,
) -> pd.DataFrame:
    """Run PyDESeq2 for sample1 vs sample2 using pseudo-replicates."""
    # Filter to the two samples' replicates
    mask = meta_df["condition"].isin([sample1, sample2])
    sub_counts = counts_df.loc[mask]
    sub_meta = meta_df.loc[mask]

    if sub_counts.shape[0] < 4:
        return pd.DataFrame()

    # Drop genes with zero counts
    sub_counts = sub_counts.loc[:, sub_counts.sum(axis=0) > 0]
    if sub_counts.shape[1] == 0:
        return pd.DataFrame()

    # PyDESeq2
    dds = DeseqDataSet(
        counts=sub_counts.astype(int),
        metadata=sub_meta,
        design="~condition",
    )
    dds.deseq2()

    # Results
    contrast = ["condition", sample2, sample1]
    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()

    res = stat_res.results_df.copy()
    res["gene"] = res.index
    res["cell_type"] = cell_type
    res["sample1"] = sample1
    res["sample2"] = sample2
    res = res.sort_values("pvalue")

    return res


def plot_volcano(df: pd.DataFrame, title: str, out_png: Path, top_n: int = 5):
    """Volcano plot with top N up/down gene labels."""
    fig, ax = plt.subplots(figsize=(7, 5.5))

    df = df.copy()
    df["neg_log10p"] = -np.log10(df["padj"].clip(lower=1e-300))

    # Classify
    sig_up = (df["padj"] < 0.05) & (df["log2FoldChange"] > 1)
    sig_down = (df["padj"] < 0.05) & (df["log2FoldChange"] < -1)
    ns = ~(sig_up | sig_down)

    # Plot NS
    ax.scatter(
        df.loc[ns, "log2FoldChange"], df.loc[ns, "neg_log10p"],
        c="#CCCCCC", s=6, alpha=0.4, linewidths=0, rasterized=True,
    )
    # Plot sig up
    ax.scatter(
        df.loc[sig_up, "log2FoldChange"], df.loc[sig_up, "neg_log10p"],
        c="#D62728", s=10, alpha=0.6, linewidths=0, rasterized=True,
    )
    # Plot sig down
    ax.scatter(
        df.loc[sig_down, "log2FoldChange"], df.loc[sig_down, "neg_log10p"],
        c="#1F77B4", s=10, alpha=0.6, linewidths=0, rasterized=True,
    )

    # Reference lines
    ax.axhline(-np.log10(0.05), ls="--", c="#888888", lw=0.6, zorder=0)
    ax.axvline(-1, ls="--", c="#888888", lw=0.6, zorder=0)
    ax.axvline(1, ls="--", c="#888888", lw=0.6, zorder=0)

    # Label top N up / down genes with auto-repulsion
    from adjustText import adjust_text

    texts = []
    if sig_up.any():
        top_up = df.loc[sig_up].nlargest(top_n, "neg_log10p")
        for _, row in top_up.iterrows():
            texts.append(ax.text(
                row["log2FoldChange"], row["neg_log10p"], row["gene"],
                fontsize=6, fontstyle="italic", color="#D62728", ha="center", va="center",
            ))
    if sig_down.any():
        top_down = df.loc[sig_down].nlargest(top_n, "neg_log10p")
        for _, row in top_down.iterrows():
            texts.append(ax.text(
                row["log2FoldChange"], row["neg_log10p"], row["gene"],
                fontsize=6, fontstyle="italic", color="#1F77B4", ha="center", va="center",
            ))

    if texts:
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
            force_text=(0.8, 0.8), force_points=(0.5, 0.5),
            expand=(1.2, 1.4), max_move=None,
        )

    # Count
    n_up = sig_up.sum()
    n_down = sig_down.sum()
    n_ns = ns.sum()

    # Legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#D62728",
               markersize=6, label=f"Up ({n_up})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1F77B4",
               markersize=6, label=f"Down ({n_down})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#CCCCCC",
               markersize=6, label=f"NS ({n_ns})"),
    ]
    ax.legend(handles=handles, fontsize=7, loc="upper right",
              frameon=True, framealpha=0.8, edgecolor="#CCCCCC")

    ax.set_xlabel("log$_2$ Fold Change", fontsize=9)
    ax.set_ylabel("-log$_{10}$(padj)", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    fig.savefig(str(out_png), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, nargs="+",
                   help="Annotated h5ad files (with layers['counts'])")
    p.add_argument("--output-dir", required=True,
                   help="Output directory")
    p.add_argument("--min-cells", type=int, default=30,
                   help="Min cells per sample for a cell type (default: 30)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for pseudo-replicate splitting")
    p.add_argument("--top-n", type=int, default=5,
                   help="Number of top up/down genes to label on volcano plot (default: 5)")
    args = p.parse_args()

    np.random.seed(args.seed)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    comparisons = {
        "ovaries": [
            ("luanchao-21310-10XSC3", "luanchao-11238-10XSC3"),
        ],
        "Uterus": [
            ("zigong-21310-10XSC3", "ZIGONG-21224-10XSC3"),
            ("zigong-11238-10XSC3", "zigong-1411200-10XSC3"),
            ("zigong-11238-10XSC3", "zigong-2200671-10XSC3"),
        ],
    }

    for h5ad_path in args.input:
        path = Path(h5ad_path)
        label = path.stem

        if "ovaries" in label.lower():
            tissue = "ovaries"
        elif "uterus" in label.lower():
            tissue = "Uterus"
        else:
            print(f"SKIP (unknown tissue): {path.name}")
            continue

        counter = "Cell Ranger" if "cellranger" in label.lower() else "scTE"
        print(f"\n{'='*60}")
        print(f"Processing: {path.name} ({tissue}, {counter})")
        print(f"{'='*60}")

        adata = ad.read_h5ad(str(path))
        cell_types = sorted(adata.obs["cell_type"].unique().tolist())
        print(f"  Cell types: {len(cell_types)}")

        comp_dir = out_dir / label
        comp_dir.mkdir(exist_ok=True)

        for ct in cell_types:
            adata_ct = adata[adata.obs["cell_type"] == ct].copy()
            sample_counts = adata_ct.obs["sample_id"].value_counts()
            print(f"\n  {ct}: {len(adata_ct)} cells")

            valid_samples = sample_counts[sample_counts >= args.min_cells].index.tolist()
            if len(valid_samples) < 2:
                print(f"    SKIP (<2 samples with >={args.min_cells} cells)")
                continue

            # Pseudobulk with pseudo-replicates
            pb_counts, pb_meta = pseudobulk_with_replicates(adata_ct)
            # Keep only valid samples
            keep_mask = pb_meta["condition"].isin(valid_samples)
            pb_counts = pb_counts.loc[keep_mask]
            pb_meta = pb_meta.loc[keep_mask]

            for s1, s2 in comparisons[tissue]:
                if s1 not in valid_samples or s2 not in valid_samples:
                    print(f"    SKIP {s1.split('-')[1]} vs {s2.split('-')[1]}: sample missing")
                    continue

                comp_name = f"{s1.split('-')[1]}_vs_{s2.split('-')[1]}"
                print(f"    DEG: {comp_name} ...", end=" ", flush=True)

                try:
                    res = run_deg(pb_counts, pb_meta, s1, s2, ct)
                except Exception as e:
                    print(f"ERROR: {e}")
                    continue

                if res.empty:
                    print("EMPTY")
                    continue

                # Save TSV
                tsv_path = comp_dir / f"{ct}_{comp_name}_deg.tsv"
                res.to_csv(tsv_path, sep="\t", index=False)

                sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)]
                print(f"sig={len(sig)}")

                # Volcano
                png_path = comp_dir / f"{ct}_{comp_name}_volcano.png"
                title = f"{ct} — {comp_name} ({counter})"
                plot_volcano(res, title, png_path, args.top_n)

        adata.file.close()

    print(f"\nAll done. Results in: {out_dir}")


if __name__ == "__main__":
    main()
