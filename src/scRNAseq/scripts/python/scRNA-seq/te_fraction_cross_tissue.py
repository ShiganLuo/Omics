#!/usr/bin/env python3
"""TE% cross-tissue summary (Q5): per-cell-type TE fraction changes from the
pseudo-bulk counts tables (diff/pseudobulk).  TE% = TE counts / (TE + gene
counts), aggregated per sample x cell type.

Outputs te_fraction_cross_tissue.png: left = Ovaries TE% change bars
(aged/youth), right = Uterus TE% log2FC heatmap across 3 contrasts.

Usage:
  python te_fraction_cross_tissue.py \
      --root /data/16T/luosg/output/luanchao_mulatta/scRNAseq \
      --out /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/figures_te/te_fraction_cross_tissue.png
"""
import argparse, glob, re, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


def te_pct_table(root, tissue):
    rows = {}
    for f in glob.glob(f"{root}/diff/pseudobulk/{tissue}/scTE/"
                       f"pseudobulk_counts_te/{tissue}_*_counts.csv"):
        ct = re.match(rf"{tissue}_(.+)_counts\.csv", f.split("/")[-1]).group(1)
        te = pd.read_csv(f)
        gn = pd.read_csv(f.replace("counts_te", "counts_gene"))
        meta = ["sample", "cell_type", "n_cells", "n_counts"]
        te["samp"] = te["sample"].str.replace(r"_rep\d+$", "", regex=True)
        gn["samp"] = gn["sample"].str.replace(r"_rep\d+$", "", regex=True)
        t = te.drop(columns=meta).groupby("samp").sum(numeric_only=True).sum(axis=1)
        g = gn.drop(columns=meta).groupby("samp").sum(numeric_only=True).sum(axis=1)
        rows[ct.replace("_", " ")] = (t / (t + g) * 100).round(2)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    ov = te_pct_table(args.root, "ovaries")
    ut = te_pct_table(args.root, "Uterus")
    lfc_ut = pd.DataFrame({
        "sham\u2192OVX":    np.log2(ut.loc["ZIGONG-21224-10XSC3"] / ut.loc["zigong-21310-10XSC3"]),
        "Non-OVX\u2192OVX": np.log2(ut.loc["zigong-1411200-10XSC3"] / ut.loc["zigong-11238-10XSC3"]),
        "OVX\u2192OVX+OT":  np.log2(ut.loc["zigong-2200671-10XSC3"] / ut.loc["zigong-1411200-10XSC3"]),
    }).round(2)
    lfc_ov = (np.log2(ov.loc["luanchao-11238-10XSC3"] /
                      ov.loc["luanchao-21310-10XSC3"])).round(2).dropna().sort_values()

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2),
                             gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axes[0]
    vals, names = lfc_ov.values, lfc_ov.index
    colors = ["#c0392b" if v > 0 else "#2471a3" for v in vals]
    ax.barh(range(len(vals)), vals, color=colors, height=0.62)
    for i, v in enumerate(vals):
        ax.text(v + (0.012 if v >= 0 else -0.012), i, f"{v:+.2f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=9)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(names, fontsize=9.5)
    ax.set_xlim(min(vals) - 0.14, max(vals) + 0.14)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("TE% log2FC (aged / youth)", fontsize=11)
    up, dn = int((vals > 0).sum()), int((vals < 0).sum())
    ax.set_title(f"Ovaries: TE% change \u2014 {up} up / {dn} down\n"
                 "(aged vs youth, per cell type)", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    data = lfc_ut.values
    cmap = LinearSegmentedColormap.from_list("rb", ["#2471a3", "#f5f5f5", "#c0392b"])
    im = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, vmin=-0.5, vmax=0.5,
                   aspect="auto")
    ax.set_xticks(range(3))
    ax.set_xticklabels(lfc_ut.columns, fontsize=10)
    ax.set_yticks(range(len(lfc_ut.index)))
    ax.set_yticklabels(lfc_ut.index, fontsize=9.5)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=9,
                        color="#888")
            else:
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=9.5, color="white" if abs(v) > 0.32 else "black")
    ax.set_title("Uterus: TE% log2FC across 3 contrasts\n"
                 "(red = up in second group)", fontsize=12)
    plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02).set_label("log2FC", fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300, facecolor="white")
    print("saved", args.out)


if __name__ == "__main__":
    main()
