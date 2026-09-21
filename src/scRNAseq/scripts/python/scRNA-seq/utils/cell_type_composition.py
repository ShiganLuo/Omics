#!/usr/bin/env python3
"""分析各细胞群的样本来源比例。

读取 annotated h5ad，按 cell_type 计算各 sample_id 的占比，
输出 TSV 表格 + stacked bar 图。
"""

import argparse
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

PALETTE = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#00468B", "#42B540", "#0099B4", "#AD002A", "#6A3D9A",
    "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
]


def load_obs(h5ad_path: str) -> pd.DataFrame:
    """只读 obs 的 cell_type 和 sample_id 列。"""
    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs[["cell_type", "sample_id"]].copy()
    adata.file.close()
    return obs


def compute_composition(obs: pd.DataFrame) -> tuple:
    """按 cell_type 计算各 sample 占比（行归一化）。

    Returns:
        counts_df: 绝对数 (cell_type × sample_id)
        pct_df:    百分比 (cell_type × sample_id)
    """
    counts_df = pd.crosstab(obs["cell_type"], obs["sample_id"])
    pct_df = counts_df.div(counts_df.sum(axis=1), axis=0) * 100
    return counts_df, pct_df


def plot_stacked_bar(pct_df: pd.DataFrame, title: str, out_png: Path):
    """每个 cell_type 一个 bar，stacked by sample。"""
    samples = pct_df.columns.tolist()
    cell_types = pct_df.index.tolist()
    colors = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(samples)}

    fig, ax = plt.subplots(figsize=(10, max(4, len(cell_types) * 0.45)))
    bottom = np.zeros(len(cell_types))

    for sample in samples:
        vals = pct_df[sample].values
        ax.barh(cell_types, vals, left=bottom, label=sample,
                color=colors[sample], edgecolor="white", linewidth=0.3)
        # 标注百分比（>5%才显示）
        for i, v in enumerate(vals):
            if v >= 5:
                ax.text(bottom[i] + v / 2, i, f"{v:.0f}%",
                        ha="center", va="center", fontsize=7, color="white",
                        fontweight="bold")
        bottom += vals

    ax.set_xlabel("Sample proportion (%)")
    ax.set_xlim(0, 100)
    ax.legend(title="Sample", bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=8)
    ax.set_title(title, fontsize=11)
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(str(out_png), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}")


def run_one(h5ad_path: str, out_dir: Path):
    path = Path(h5ad_path)
    name = path.stem.replace("_annotated", "").replace("_v2", "").replace("_v3", "")
    label = path.stem  # e.g. Uterus_cellranger_v3_annotated

    print(f"\n{'='*60}")
    print(f"Processing: {path.name}")
    print(f"{'='*60}")

    obs = load_obs(str(path))
    print(f"  Cells: {len(obs)} | Cell types: {obs['cell_type'].nunique()} | Samples: {obs['sample_id'].nunique()}")

    counts_df, pct_df = compute_composition(obs)

    # 输出 TSV
    counts_tsv = out_dir / f"{label}_counts.tsv"
    pct_tsv = out_dir / f"{label}_pct.tsv"
    counts_df.to_csv(counts_tsv, sep="\t")
    pct_df.round(2).to_csv(pct_tsv, sep="\t")
    print(f"  -> {counts_tsv}")
    print(f"  -> {pct_tsv}")

    # 画图
    png = out_dir / f"{label}_sample_composition.png"
    tissue = "Ovary" if "ovaries" in label.lower() else "Uterus"
    counter = "Cell Ranger" if "cellranger" in label.lower() else "scTE"
    title = f"{tissue} ({counter}) — Cell type sample composition"
    plot_stacked_bar(pct_df, title, png)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, nargs="+",
                   help="Annotated h5ad files")
    p.add_argument("--output-dir", required=True,
                   help="Output directory for tables and plots")
    args = p.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for h5ad_path in args.input:
        run_one(h5ad_path, out_dir)

    print(f"\nAll done. Results in: {out_dir}")


if __name__ == "__main__":
    main()
