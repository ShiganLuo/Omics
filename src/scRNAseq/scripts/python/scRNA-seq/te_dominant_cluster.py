#!/usr/bin/env python3
"""TE-dominant cluster analysis (Q6, re-run on the current annotation).

For each tissue: compute per-cell TE fraction from the merged scTE counts
(gene + TE features), aggregate per leiden cluster, identify TE-dominant
clusters (mean TE fraction > 35%), then characterise them:

  fig1  UMAP colored by per-cell TE% + TE-dominant clusters outlined
  fig2  TE% per cluster (dorminant highlighted) + per-sample proportion of
        the dominant cluster
  fig3  TE class composition (LTR / LINE / SINE / DNA / other) dominant vs rest
  fig4  marker evidence: senescence / SASP / cell-cycle / lineage markers

Inputs (h5ad): merged (counts layer, genes+TEs) + auto (leiden, cell_type,
X_umap). The TE/gene feature split comes from the pseudo-bulk counts CSV
column lists (pseudobulk_counts_te / pseudobulk_counts_gene).

Usage:
  python te_dominant_cluster.py \
      --root /data/16T/luosg/output/luanchao_mulatta/scRNAseq \
      --out-dir /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/figures_te
"""

import argparse
import os
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

MARKERS = {
    "senescence": ["CDKN1A", "CDKN2A", "GLB1", "LMNB1"],
    "SASP": ["IL6", "CXCL8", "CCL2", "MMP9", "IL1B"],
    "cell_cycle": ["MKI67", "TOP2A", "PCNA"],
}


def te_gene_split(root, tissue):
    import glob
    te_cols, gn_cols = set(), set()
    for f in glob.glob(f"{root}/diff/pseudobulk/{tissue}/scTE/pseudobulk_counts_te/{tissue}_*_counts.csv"):
        te_cols |= set(pd.read_csv(f, nrows=1).columns)
    for f in glob.glob(f"{root}/diff/pseudobulk/{tissue}/scTE/pseudobulk_counts_gene/{tissue}_*_counts.csv"):
        gn_cols |= set(pd.read_csv(f, nrows=1).columns)
    meta = {"sample", "cell_type", "n_cells", "n_counts", "samp"}
    te_cols -= meta
    gn_cols -= meta
    return te_cols, gn_cols


def te_class(name):
    n = name.upper()
    if re.search(r"ALU|B1|B2|ID_|SINE", n):
        return "SINE"
    if re.search(r"^L1|L1[A-Z]|L2|CR1|LINE|LX", n):
        return "LINE"
    if re.search(r"LTR|ERV|MALR|MER[V0-9]|MLT|THE1|THE[0-9]|ERVL|HERV|MACERV|MAcERV", n):
        return "LTR"
    if re.search(r"HAT|TCMAR|TIGGER|MULE|HARBINGER|PIGGYBAC|DNA", n):
        return "DNA"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    summary_rows = []

    for tissue in ("Uterus", "ovaries"):
        merged = sc.read_h5ad(
            f"{args.root}/common/5_combine_h5ad/{tissue}/{tissue}_scTE_merged.h5ad")
        auto = sc.read_h5ad(
            f"{args.root}/common/5_combine_h5ad/{tissue}/{tissue}_scTE_auto.h5ad")
        te_cols, gn_cols = te_gene_split(args.root, tissue)

        common = merged.obs_names.intersection(auto.obs_names)
        merged = merged[common].copy()
        auto = auto[common].copy()

        feats = merged.var_names
        te_idx = np.array([i for i, f in enumerate(feats) if f in te_cols])
        gn_idx = np.array([i for i, f in enumerate(feats) if f in gn_cols])
        X = merged.layers["counts"] if "counts" in merged.layers else merged.X
        X = np.asarray(X.todense() if hasattr(X, "todense") else X)
        te_sum = X[:, te_idx].sum(1)
        gn_sum = X[:, gn_idx].sum(1)
        te_frac = te_sum / np.maximum(te_sum + gn_sum, 1)

        df = pd.DataFrame({
            "te_frac": te_frac,
            "leiden": auto.obs["leiden"].astype(str).values,
            "cell_type": auto.obs["cell_type"].astype(str).values,
            "sample_id": auto.obs["sample_id"].astype(str).values,
        }, index=common)

        cl = df.groupby("leiden").agg(
            te_frac=("te_frac", "mean"), n=("te_frac", "size"),
            cell_type=("cell_type", lambda s: s.value_counts().index[0]))
        cl["te_frac"] *= 100
        dom = cl[cl["te_frac"] > 35].sort_values("te_frac", ascending=False)
        print(f"== {tissue}: {len(cl)} clusters; TE-dominant: "
              f"{[(i, round(r.te_frac,1), r.cell_type, int(r.n)) for i, r in dom.iterrows()]}")
        for i, r in dom.iterrows():
            summary_rows.append((tissue, i, r.cell_type, round(r.te_frac, 1), int(r.n)))

        dom_clusters = set(dom.index)

        # ---- fig1: UMAP by TE% + dominant outline -----------------------
        fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
        emb = auto.obsm["X_umap"]
        ax = axes[0]
        sc_pt = ax.scatter(emb[:, 0], emb[:, 1], c=df["te_frac"] * 100,
                           s=1.2, cmap="RdYlBu_r", vmin=0, vmax=60)
        plt.colorbar(sc_pt, ax=ax, fraction=0.046).set_label("TE% per cell")
        ax.set_title(f"{tissue}: per-cell TE fraction", fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])
        ax = axes[1]
        dom_mask = df["leiden"].isin(dom_clusters).values
        ax.scatter(emb[~dom_mask, 0], emb[~dom_mask, 1], s=1.2, c="#dddddd")
        # red-yellow-blue categorical ramp
        cmap = mcolors.LinearSegmentedColormap.from_list(
            "ryb", ["#c0392b", "#e67e22", "#f4d03f", "#5dade2", "#2471a3"])
        for k, cl_id in enumerate(dom.index):
            m = (df["leiden"] == cl_id).values
            ax.scatter(emb[m, 0], emb[m, 1], s=2.2,
                       color=cmap(k % 10), label=f"cluster {cl_id} "
                       f"({dom.loc[cl_id,'cell_type']}, {dom.loc[cl_id,'te_frac']:.0f}% TE)")
        ax.legend(fontsize=8, markerscale=3, loc="best")
        ax.set_title(f"{tissue}: TE-dominant clusters (>35% TE)", fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/te_dom_{tissue}_umap.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig2a: TE% per cluster -----------------------------------
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        cl_sorted = cl.sort_values("te_frac", ascending=False)
        colors = ["#c0392b" if i in dom_clusters else "#b0b0b0"
                  for i in cl_sorted.index]
        ax.bar(range(len(cl_sorted)), cl_sorted["te_frac"], color=colors)
        gmean = float(cl["te_frac"].mean())
        ax.axhline(gmean, ls="--", lw=1.2, color="#333333",
                   label=f"global mean {gmean:.1f}%")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_xticks(range(len(cl_sorted)))
        ax.set_xticklabels(cl_sorted.index, rotation=90, fontsize=8)
        ax.set_ylabel("mean TE% per cluster")
        ax.set_title(f"{tissue}: TE% per leiden cluster", fontsize=12)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/te_dom_{tissue}_cluster_te.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig2b: per-sample share ---------------------------------
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        n_dom = len(dom.index)
        bw = 0.8 / max(n_dom, 1)
        share_all = []
        for k, cl_id in enumerate(dom.index):
            sub = df[df["leiden"] == cl_id]
            frac = sub.groupby("sample_id").size()
            tot = df.groupby("sample_id").size()
            share = (frac / tot * 100).reindex(tot.index, fill_value=0)
            share_all.extend(share.values.tolist())
            ax.bar(np.arange(len(share)) + (k - (n_dom - 1) / 2) * bw,
                   share.values, bw * 0.92,
                   color=cmap(k % 10), label=f"cluster {cl_id}")
        if share_all:
            ax.axhline(np.mean(share_all), ls="--", lw=1.2, color="#333333",
                       label=f"global mean {np.mean(share_all):.1f}%")
        ax.set_xticks(range(len(tot)))
        ax.set_xticklabels(tot.index, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("% of cells in cluster")
        ax.set_title("TE-dominant cluster share per sample", fontsize=12)
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/te_dom_{tissue}_share.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig3: TE class composition --------------------------------
        fig, ax = plt.subplots(figsize=(8.6, 4.8))
        classes = ["LTR", "LINE", "SINE", "DNA", "other"]
        # red-yellow-blue ramp, consistent with the UMAP panels
        ccolors = ["#c0392b", "#e67e22", "#f4d03f", "#5dade2", "#2471a3"]
        all_te_names = [f for f in feats if f in te_cols]
        cls_of = {f: te_class(f) for f in all_te_names}
        def comp(mask):
            sel = X[mask][:, [i for i, f in enumerate(feats) if f in te_cols]]
            tot = sel.sum(0)
            s = pd.Series(tot, index=[feats[i] for i, f in enumerate(feats) if f in te_cols])
            by = {}
            for f, v in s.items():
                by[cls_of[f]] = by.get(cls_of[f], 0) + v
            t = sum(by.values()) or 1
            return [by.get(c, 0) / t * 100 for c in classes]
        dom_m = df["leiden"].isin(dom_clusters).values
        rest_m = ~dom_m
        xs = np.arange(len(classes))
        ax.bar(xs - 0.18, comp(dom_m), 0.36, color=ccolors,
               label="TE-dominant clusters")
        ax.bar(xs + 0.18, comp(rest_m), 0.36, color=ccolors, alpha=0.45,
               label="other clusters")
        ax.set_xticks(xs); ax.set_xticklabels(classes, fontsize=11)
        ax.set_ylabel("% of TE counts")
        ax.set_title(f"{tissue}: TE class composition", fontsize=12)
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/te_dom_{tissue}_class.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig4: markers ---------------------------------------------
        rows = []
        for group, names in MARKERS.items():
            for g in names:
                hits = [f for f in auto.var_names if f == g or f.startswith(g + ".")]
                if not hits:
                    continue
                gi = auto.var_names.get_loc(hits[0])
                Xa = auto.X[:, gi]
                Xa = np.asarray(Xa.todense() if hasattr(Xa, "todense") else Xa).ravel()
                rows.append((group, g,
                             Xa[dom_m].mean(),
                             Xa[rest_m].mean()))
        if rows:
            mk = pd.DataFrame(rows, columns=["group", "gene", "dom", "rest"])
            fig, ax = plt.subplots(figsize=(10, 4.6))
            xs = np.arange(len(mk))
            ax.bar(xs - 0.2, mk["dom"], 0.4, color="#c0392b",
                   label="TE-dominant clusters")
            ax.bar(xs + 0.2, mk["rest"], 0.4, color="#7f8c8d",
                   label="other clusters")
            gm = float(np.mean(np.concatenate([mk["dom"].to_numpy(),
                                               mk["rest"].to_numpy()])))
            ax.axhline(gm, ls="--", lw=1.2, color="#333333",
                       label=f"global mean {gm:.2f}")
            ax.set_xticks(xs)
            ax.set_xticklabels([f"{r.gene}\n({r.group})" for r in mk.itertuples()],
                               fontsize=8.5)
            ax.set_ylabel("mean log1p expression")
            ax.set_title(f"{tissue}: senescence / SASP / cell-cycle markers",
                         fontsize=12)
            ax.legend(fontsize=9)
            ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            fig.savefig(f"{args.out_dir}/te_dom_{tissue}_markers.png", dpi=300,
                        facecolor="white")
            plt.close(fig)

    pd.DataFrame(summary_rows,
                 columns=["tissue", "cluster", "cell_type", "te_pct",
                          "n_cells"]).to_csv(
        f"{args.out_dir}/te_dominant_clusters.csv", index=False)

    # ---- molecular evidence: top up TEs + ribosomal/histone down --------
    import glob as _glob
    for tissue in ("Uterus", "ovaries"):
        dom_ct = [r[2] for r in summary_rows if r[0] == tissue]
        if not dom_ct:
            continue
        ct_dir = dom_ct[0].replace(" ", "_")
        te_files = _glob.glob(f"{args.root}/diff/pseudobulk/{tissue}/scTE/"
                              f"deg_pseudobulk_te/deg_*/*/{ct_dir}/*_DEG.csv")
        gn_files = _glob.glob(f"{args.root}/diff/pseudobulk/{tissue}/scTE/"
                              f"deg_pseudobulk_gene/deg_*/*/{ct_dir}/*_DEG.csv")
        if not te_files or not gn_files:
            continue
        te = pd.read_csv(te_files[0])
        te = te[te["significant"] == True]  # noqa: E712
        up_te = te[te["log2FC"] > 0].nlargest(6, "log2FC")
        gn = pd.read_csv(gn_files[0])
        gn = gn[gn["significant"] == True]  # noqa: E712
        down = gn[gn["log2FC"] < 0].copy()
        down["kind"] = down["gene"].map(
            lambda g: "ribosomal" if g.upper().startswith(("RPL", "RPS"))
            else ("histone" if g.upper().startswith("HIST") else "other"))
        down_ribo = down[down["kind"] != "other"].nsmallest(6, "log2FC")
        fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.4))
        ax = axes[0]
        ax.barh(range(len(up_te)), up_te["log2FC"], color="#c0392b")
        ax.set_yticks(range(len(up_te)))
        ax.set_yticklabels(up_te["gene"], fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("log2FC")
        ax.set_title(f"{tissue} / {dom_ct[0]}: top up TEs", fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        ax = axes[1]
        if len(down_ribo):
            cols = ["#2471a3" if k == "ribosomal" else "#5dade2"
                    for k in down_ribo["kind"]]
            ax.barh(range(len(down_ribo)), down_ribo["log2FC"], color=cols)
            ax.set_yticks(range(len(down_ribo)))
            ax.set_yticklabels(down_ribo["gene"], fontsize=9)
            ax.invert_yaxis()
            ax.legend(handles=[Patch(color="#2471a3", label="ribosomal"),
                               Patch(color="#5dade2", label="histone")],
                      fontsize=8)
        ax.set_xlabel("log2FC")
        ax.set_title(f"{tissue} / {dom_ct[0]}: down ribosomal/histone",
                     fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/te_dom_{tissue}_evidence.png", dpi=300,
                    facecolor="white")
        plt.close(fig)
    print("done ->", args.out_dir)


if __name__ == "__main__":
    main()
