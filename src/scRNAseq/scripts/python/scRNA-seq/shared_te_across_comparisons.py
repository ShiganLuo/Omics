#!/usr/bin/env python3
"""Shared up/down transposable elements across multiple comparisons.

Reads the scTE pseudo-bulk TE DEG tables produced by modules/scanpy
(``deg_pseudobulk_te``) under the workflow diff tree and finds TEs that are
significantly up- or down-regulated in the SAME direction across multiple
comparisons:

* per tissue: TEs shared by >=2 (and all) contrasts, aggregated over cell
  types (max |log2FC| per contrast) and per cell type (exact intersections);
* cross tissue: TEs shared between the ovaries aging contrast and the
  uterus contrasts.

Outputs (in --out-dir):
  shared_TEs_across_comparisons.csv   per (tissue, cell_type, direction) lists
  shared_TEs_global.csv               TE-level cross-contrast summary
  shared_TEs_summary.png              overview figure (counts + heatmap)

Usage:
  python shared_te_across_comparisons.py \
      --diff-root /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/pseudobulk \
      --out-dir   /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/figures_te
"""

import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CONTRAST_LABELS = {
    "Uterus_zigong-21310-10XSC3_vs_ZIGONG-21224-10XSC3": "sham\u2192OVX",
    "Uterus_zigong-11238-10XSC3_vs_zigong-1411200-10XSC3": "Non-OVX\u2192OVX",
    "Uterus_zigong-1411200-10XSC3_vs_zigong-2200671-10XSC3": "OVX\u2192OVX+OT",
    "ovaries_luanchao-21310-10XSC3_vs_luanchao-11238-10XSC3": "youth\u2192aged",
}


def load_tables(diff_root):
    """tissue -> cell_type -> contrast -> significant-only DataFrame."""
    out = {}
    pat = os.path.join(diff_root, "*", "scTE", "deg_pseudobulk_te",
                       "deg_*", "*", "*", "*_DEG.csv")
    for f in glob.glob(pat):
        rel = os.path.relpath(f, diff_root).split(os.sep)
        tissue, contrast, cell_type = rel[0], rel[4], rel[5].replace("_", " ")
        df = pd.read_csv(f)
        if "significant" in df.columns:
            df = df[df["significant"] == True]  # noqa: E712
        out.setdefault(tissue, {}).setdefault(cell_type, {})[contrast] = df[
            ["gene", "log2FC", "pval_adj"]].copy()
    return out


def te_level_summary(tables):
    """TE -> per-contrast max-|log2FC| direction table (per tissue)."""
    rows = []
    for tissue, ctypes in tables.items():
        for cell_type, contrasts in ctypes.items():
            for contrast, df in contrasts.items():
                for r in df.itertuples():
                    rows.append((tissue, contrast, r.gene, r.log2FC, cell_type))
    long = pd.DataFrame(rows, columns=["tissue", "contrast", "te", "log2FC", "cell_type"])
    if long.empty:
        return long
    idx = long.groupby(["tissue", "te", "contrast"])["log2FC"].apply(
        lambda s: s.loc[s.abs().idxmax()])
    wide = idx.unstack("contrast")
    return wide


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diff-root", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    tables = load_tables(args.diff_root)
    wide = te_level_summary(tables)

    # ---- TE-level sharing across contrasts (per tissue) -------------------
    recs = []
    for tissue in wide.index.get_level_values(0).unique():
        sub = wide.loc[tissue]
        n = sub.shape[1]
        ups = sub.gt(0).sum(axis=1)
        dns = sub.lt(0).sum(axis=1)
        for te in sub.index:
            k_up, k_dn = int(ups[te]), int(dns[te])
            if k_up >= 2 or k_dn >= 2:
                recs.append({
                    "tissue": tissue, "te": te, "n_contrasts": n,
                    "n_up": k_up, "n_down": k_dn,
                    "direction": "up" if k_up >= 2 else "down",
                    "shared_contrasts": "; ".join(
                        CONTRAST_LABELS.get(c, c) for c in sub.columns[
                            (sub.loc[te] > 0 if k_up >= 2 else sub.loc[te] < 0)]),
                    "log2FCs": "; ".join(
                        f"{CONTRAST_LABELS.get(c, c)}:{v:+.2f}"
                        for c, v in sub.loc[te].items()),
                })
    global_df = pd.DataFrame(recs)
    if len(global_df):
        global_df["n_shared"] = global_df[["n_up", "n_down"]].max(axis=1)
        global_df = global_df.sort_values(
            ["tissue", "direction", "n_shared"], ascending=[True, True, False])
    global_df.to_csv(os.path.join(args.out_dir, "shared_TEs_global.csv"), index=False)

    # ---- per-cell-type exact intersections -------------------------------
    ct_recs = []
    for tissue, ctypes in tables.items():
        for cell_type, contrasts in ctypes.items():
            if len(contrasts) < 2:
                continue
            up_sets = {c: set(d[d["log2FC"] > 0]["gene"]) for c, d in contrasts.items()}
            dn_sets = {c: set(d[d["log2FC"] < 0]["gene"]) for c, d in contrasts.items()}
            all_c = list(contrasts)
            for k in (2, len(all_c)):
                if k < 2:
                    continue
                import itertools
                for combo in itertools.combinations(all_c, k):
                    up = set.intersection(*(up_sets[c] for c in combo))
                    dn = set.intersection(*(dn_sets[c] for c in combo))
                    for direction, s in (("up", up), ("down", dn)):
                        for te in sorted(s):
                            ct_recs.append({
                                "tissue": tissue, "cell_type": cell_type,
                                "direction": direction, "n_contrasts": k,
                                "contrasts": "; ".join(CONTRAST_LABELS.get(c, c)
                                                        for c in combo),
                                "te": te,
                            })
    ct_df = pd.DataFrame(ct_recs)
    ct_df.to_csv(os.path.join(args.out_dir, "shared_TEs_across_comparisons.csv"),
                 index=False)

    # ---- figure ---------------------------------------------------------
    fig = plt.figure(figsize=(13.5, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.35])

    ax = fig.add_subplot(gs[0])
    if len(global_df):
        cnt = global_df.groupby(["tissue", "direction", "n_shared"]).size().unstack(
            fill_value=0)
        labels, up2, up3, dn2, dn3 = [], [], [], [], []
        for tissue in sorted(global_df["tissue"].unique()):
            sub = global_df[global_df["tissue"] == tissue]
            labels.append(tissue)
            up2.append(len(sub[(sub.direction == "up") & (sub.n_up == 2)]))
            up3.append(len(sub[(sub.direction == "up") & (sub.n_up >= 3)]))
            dn2.append(len(sub[(sub.direction == "down") & (sub.n_down == 2)]))
            dn3.append(len(sub[(sub.direction == "down") & (sub.n_down >= 3)]))
        x = np.arange(len(labels))
        ax.bar(x - 0.27, up2, 0.18, color="#e57373", label="up in 2")
        ax.bar(x - 0.09, up3, 0.18, color="#c0392b", label="up in 3+")
        ax.bar(x + 0.09, dn2, 0.18, color="#79a9d6", label="down in 2")
        ax.bar(x + 0.27, dn3, 0.18, color="#2471a3", label="down in 3+")
        for xi, vals in zip(x, zip(up2, up3, dn2, dn3)):
            for off, v in zip((-0.27, -0.09, 0.09, 0.27), vals):
                if v:
                    ax.text(xi + off, v + 0.3, str(v), ha="center", fontsize=9)
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10)
        ax.legend(fontsize=9)
    ax.set_ylabel("number of TEs", fontsize=11)
    ax.set_title("TEs shared across comparisons\n(significant, same direction)",
                 fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)

    ax = fig.add_subplot(gs[1])
    top = global_df.sort_values("n_shared", ascending=False).head(24)
    if len(top):
        tissues = list(dict.fromkeys(top["tissue"]))
        all_c = []
        for t in tissues:
            for c in wide.loc[t].columns:
                all_c.append(c)
        mat = []
        for r in top.itertuples():
            row = []
            for c in all_c:
                try:
                    v = wide.loc[(r.tissue, r.te), c]
                except KeyError:
                    v = np.nan
                row.append(v)
            mat.append(row)
        mat = np.array(mat, dtype=float)
        vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1
        im = ax.imshow(np.ma.masked_invalid(mat), cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(all_c)))
        ax.set_xticklabels([CONTRAST_LABELS.get(c, c) for c in all_c],
                           rotation=30, ha="right", fontsize=9)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([f"{r.te} ({r.tissue[:3]})" for r in top.itertuples()],
                           fontsize=8.5)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:+.1f}", ha="center", va="center",
                            fontsize=7.5,
                            color="white" if abs(v) > 0.5 * vmax else "black")
                else:
                    ax.text(j, i, "\u2013", ha="center", va="center",
                            fontsize=8, color="#999")
        plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("log2FC")
    ax.set_title("Top shared TEs \u2014 max |log2FC| per contrast\n"
                 "(\u2013 = not significant in that contrast)", fontsize=12)

    fig.tight_layout()
    png = os.path.join(args.out_dir, "shared_TEs_summary.png")
    fig.savefig(png, dpi=300, facecolor="white")
    print("saved", png)
    print(global_df.to_string(index=False)[:3000])
    print(ct_df.groupby(["tissue", "cell_type", "direction", "n_contrasts"])
          .size().to_string()[:2000])


if __name__ == "__main__":
    main()
