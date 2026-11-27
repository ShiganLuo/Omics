#!/usr/bin/env python3
"""TE fraction cross-tissue summary: per-cell-type TE fraction changes per
comparison, from the pseudo-bulk counts tables (diff/pseudobulk).

Metric and unit: te_fraction = TE counts / (TE + gene counts), aggregated
per group as the ratio of summed counts -- a unitless fraction in [0, 1]
(NOT percent).  The plotted effect is log2FC = log2(exp / ctrl) of
te_fraction (unitless); positive values mean a higher TE fraction in the
experimental group.

Comparisons come from a group-pairs design (DESeq2 group_pairs shape):
each entry names the control group, the experimental group and the sample
list of each side.  Labels are rendered as "<control_group_name> vs
<experimental_group_name>".  A tissue with a single comparison is drawn as
a sorted bar panel; a tissue with several comparisons is drawn as a
heatmap (one column per comparison).

Outputs:
  te_fraction_cross_tissue.png   one panel per tissue (title = tissue name)

Usage:
  python te_fraction_cross_tissue.py \
      --root /data/16T/luosg/output/luanchao_mulatta/scRNAseq \
      --design compare_design.json \
      --out /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/figures_te/te_fraction_cross_tissue.png

  compare_design.json (group_pairs shape):
  {
    "Uterus": {
      "Sham_vs_OVX": {
        "control_group_name": "Sham",
        "experimental_group_name": "OVX",
        "control_samples": ["zigong-21310-10XSC3"],
        "experimental_samples": ["ZIGONG-21224-10XSC3"]
      }
    }
  }
"""
import argparse
import glob
import json
import re
import warnings

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


def te_fraction_table(root, tissue):
    """Per-sample TE fraction and total counts per cell type.

    Parameters:
        root: workflow output root containing diff/pseudobulk/<tissue>/scTE.
        tissue: tissue name (pseudobulk sub-directory name).

    Returns:
        (frac, total): DataFrames indexed by sample, columns = cell types.
        frac = TE / (TE + gene), unitless fraction in [0, 1].
        total = TE + gene counts (aggregation weight for group pooling).
    """
    frac, total = {}, {}
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
        tot = t + g
        frac[ct.replace("_", " ")] = t / tot
        total[ct.replace("_", " ")] = tot
    return pd.DataFrame(frac), pd.DataFrame(total)


def load_design(path):
    """Read the group-pairs design JSON: {tissue: {contrast: entry}}.

    Each entry requires non-empty control_group_name, experimental_group_name,
    control_samples, experimental_samples.  Raises ValueError naming the
    offending key otherwise.
    """
    with open(path) as fh:
        design = json.load(fh)
    for tissue, contrasts in design.items():
        for name, entry in contrasts.items():
            for key in ("control_group_name", "experimental_group_name",
                        "control_samples", "experimental_samples"):
                if not entry.get(key):
                    raise ValueError(
                        f"design[{tissue!r}][{name!r}] missing/empty key {key!r}")
    return design


def group_fraction(frac, total, samples):
    """Group-level te_fraction: ratio of summed counts across samples."""
    missing = [s for s in samples if s not in frac.index]
    if missing:
        raise ValueError(f"samples not in pseudo-bulk tables: {missing}")
    f, t = frac.loc[samples], total.loc[samples]
    return (f * t).sum() / t.sum()


def comparison_lfc(frac, total, entry):
    """log2(exp / ctrl) of te_fraction per cell type (unitless)."""
    ctrl = group_fraction(frac, total, entry["control_samples"])
    exp = group_fraction(frac, total, entry["experimental_samples"])
    lfc = np.log2(exp / ctrl)
    return lfc.replace([np.inf, -np.inf], np.nan)


def plot_bar(ax, lfc, entry):
    """Sorted bar panel for a tissue with a single comparison."""
    label = (f"{entry['control_group_name']} vs "
             f"{entry['experimental_group_name']}")
    lfc = lfc.dropna().sort_values()
    vals, names = lfc.values, lfc.index
    colors = ["#c0392b" if v > 0 else "#2471a3" for v in vals]
    ax.barh(range(len(vals)), vals, color=colors, height=0.62)
    for i, v in enumerate(vals):
        ax.text(v + (0.012 if v >= 0 else -0.012), i, f"{v:+.2f}", va="center",
                ha="left" if v >= 0 else "right", fontsize=9)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(names, fontsize=9.5)
    ax.set_xlim(min(vals) - 0.14, max(vals) + 0.14)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel(f"TE fraction log2FC ({label})", fontsize=11)
    exp = entry["experimental_group_name"]
    ax.text(0.0, 1.015, f"← lower in {exp}", transform=ax.transAxes,
            fontsize=8.5, color="#555", ha="left", va="bottom")
    ax.text(1.0, 1.015, f"higher in {exp} →", transform=ax.transAxes,
            fontsize=8.5, color="#555", ha="right", va="bottom")
    ax.spines[["top", "right"]].set_visible(False)


def plot_heatmap(ax, lfc_df, clim):
    """Heatmap panel for a tissue with several comparisons (one column each)."""
    data = lfc_df.values
    cmap = LinearSegmentedColormap.from_list("rb", ["#2471a3", "#f5f5f5", "#c0392b"])
    im = ax.imshow(np.ma.masked_invalid(data), cmap=cmap,
                   vmin=-clim, vmax=clim, aspect="auto")
    ax.set_xticks(range(lfc_df.shape[1]))
    ax.set_xticklabels([str(i + 1) for i in range(lfc_df.shape[1])], fontsize=11)
    ax.set_yticks(range(lfc_df.shape[0]))
    ax.set_yticklabels(lfc_df.index, fontsize=9.5)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=9,
                        color="#888")
            else:
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=9.5,
                        color="white" if abs(v) > 0.6 * clim else "black")
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("TE fraction log2FC (exp. / ctrl.)", fontsize=10)
    handles = [Line2D([], [], linestyle="none") for _ in lfc_df.columns]
    labels = [f"{i + 1}: {lab}" for i, lab in enumerate(lfc_df.columns)]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.08),
              ncol=1, frameon=False, fontsize=9.5, handlelength=0,
              handletextpad=0, labelspacing=0.35)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    ap = argparse.ArgumentParser(
        description="Per-comparison TE fraction change figure (one panel per tissue).")
    ap.add_argument("-r", "--root", required=True,
                    help="workflow output root (contains diff/pseudobulk)")
    ap.add_argument("-d", "--design", required=True,
                    help="group-pairs design JSON: {tissue: {contrast: "
                         "{control_group_name, experimental_group_name, "
                         "control_samples, experimental_samples}}}")
    ap.add_argument("-o", "--out", required=True, help="output PNG path")
    ap.add_argument("-c", "--clim", type=float, default=None,
                    help="heatmap color limit (symmetric |log2FC| bound; "
                         "default: max |log2FC| in the data)")
    args = ap.parse_args()

    design = load_design(args.design)
    tissues = list(design)
    fig, axes = plt.subplots(1, len(tissues), figsize=(6.6 * len(tissues), 5.2),
                             squeeze=False)
    for ax, tissue in zip(axes[0], tissues):
        frac, total = te_fraction_table(args.root, tissue)
        contrasts = design[tissue]
        labels = [f"{e['control_group_name']} vs {e['experimental_group_name']}"
                  for e in contrasts.values()]
        if len(contrasts) == 1:
            entry = next(iter(contrasts.values()))
            plot_bar(ax, comparison_lfc(frac, total, entry), entry)
        else:
            lfc_df = pd.DataFrame(
                {lab: comparison_lfc(frac, total, entry)
                 for lab, entry in zip(labels, contrasts.values())})
            clim = args.clim
            if clim is None:
                clim = float(np.nanmax(np.abs(lfc_df.values)))
                clim = clim if clim > 0 else 0.5
            plot_heatmap(ax, lfc_df, clim)
        ax.set_title(tissue[0].upper() + tissue[1:], fontsize=12, pad=30)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300, facecolor="white",
                bbox_inches="tight", pad_inches=0.2)
    print("saved", args.out)


if __name__ == "__main__":
    main()
