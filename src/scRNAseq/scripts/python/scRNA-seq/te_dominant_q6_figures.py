#!/usr/bin/env python3
"""Q6 TE-dominant cluster reproduction (old-PPT slides 24-28, new 23-27).

Per tissue (top-TE% leiden cluster = dominant):
  q6_<tissue>_cluster_te.png   per-cluster TE% bars (grey, dominant red,
                               dashed global mean; x = "id (celltype, n=)")
  q6_<tissue>_fam_heatmap.png  top-15 TE families (h5ad var['family_id']) x
                               clusters, value = family-proportion fold-change
                               vs global (family/total-TE share ratio)
  q6_<tissue>_umap.png         2 panels: UMAP by TE% + dominant highlighted
  q6_<tissue>_fam_bar.png      top-8 families: mean TE counts/cell dominant
                               vs others, with fold-change annotations
  q6_<tissue>_share.png        per-sample proportion: % of sample cells and
                               % within lineage (same cell_type)
  q6_<tissue>_deg.png          DESeq2 pseudo-bulk log2FC (dominant vs
                               same-lineage siblings), divergent bars coloured
                               mito / TE-LTR / TE-LINE-DNA / ribosomal / other
  q6_dotplot.png               markers x (TE / siblings / rest) per tissue

TE classification comes from merged-h5ad var['class_id'] / var['family_id']
(never from name matching).  Log2FC layer = PyDESeq2 on 3 pseudo-replicates
per side (cells randomly split), genes expressed in >=25% of cells.

Usage:
  python te_dominant_q6_figures.py \
      --root /data/16T/luosg/output/luanchao_mulatta/scRNAseq \
      --out-dir /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/figures_te
"""

import argparse
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(0)

MARKERS = {
    "SASP / inflammation": ["IL6", "CXCL8", "CXCL1", "CCL2"],
    "senescence / CDKi": ["CDKN1A", "SERPINE1", "GADD45A"],
    "proliferation": ["MKI67", "TOP2A", "CDK1"],
    "stress / apoptosis": ["BAX", "CASP3"],
}
LINEAGE = {
    "ovaries": ["AMH", "CYP19A1", "FSHR", "INHA", "GJA1"],
    "Uterus": ["EPCAM", "KRT19", "MUC1", "PAEP", "FOXA2"],
}


def load(root, tissue):
    merged = sc.read_h5ad(f"{root}/common/5_combine_h5ad/{tissue}/{tissue}_scTE_merged.h5ad")
    auto = sc.read_h5ad(f"{root}/common/5_combine_h5ad/{tissue}/{tissue}_scTE_auto.h5ad")
    common = merged.obs_names.intersection(auto.obs_names)
    merged, auto = merged[common].copy(), auto[common].copy()
    X = merged.layers["counts"] if "counts" in merged.layers else merged.X
    X = np.asarray(X.todense() if hasattr(X, "todense") else X)
    df = pd.DataFrame({
        "leiden": auto.obs["leiden"].astype(str).values,
        "cell_type": auto.obs["cell_type"].astype(str).values,
        "sample_id": auto.obs["sample_id"].astype(str).values,
    }, index=common)
    return merged, auto, X, df


def te_mask(merged):
    return merged.var["gene_type"].astype(str).values == "TE"


def fam_of(merged):
    return merged.var["family_id"].astype(str).values


def cls_of(merged):
    return merged.var["class_id"].astype(str).values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    dom_info = {}

    for tissue in ("ovaries", "Uterus"):
        merged, auto, X, df = load(args.root, tissue)
        is_te = te_mask(merged)
        fam = fam_of(merged)
        cls = cls_of(merged)
        te_sum = X[:, is_te].sum(1)
        tot_sum = X.sum(1)
        df["te_frac"] = te_sum / np.maximum(tot_sum, 1) * 100

        cl = df.groupby("leiden").agg(
            te_frac=("te_frac", "mean"), n=("te_frac", "size"),
            cell_type=("cell_type", lambda s: s.value_counts().index[0]))
        dom = cl["te_frac"].idxmax()
        dom_ct = cl.loc[dom, "cell_type"]
        dom_info[tissue] = (dom, dom_ct, float(cl.loc[dom, "te_frac"]), int(cl.loc[dom, "n"]))
        dom_mask = (df["leiden"] == dom).values

        # ---- fig: per-cluster TE% (grey + red dominant + mean line) ----
        order = cl.sort_values("te_frac", ascending=False)
        fig, ax = plt.subplots(figsize=(11, 5))
        colors = ["#c0392b" if i == dom else "#b0b0b0" for i in order.index]
        ax.bar(range(len(order)), order["te_frac"], color=colors)
        gm = float(cl["te_frac"].mean())
        ax.axhline(gm, ls="--", lw=1.2, color="#333333", label=f"global mean {gm:.1f}%")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(
            [f"{i} ({cl.loc[i,'cell_type']}, n={int(cl.loc[i,'n'])})" for i in order.index],
            rotation=60, ha="right", fontsize=8)
        ax.set_ylabel("TE fraction (%)")
        ax.set_title(f"{tissue} \u2014 Per-cluster TE fraction (raw counts)")
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/q6_{tissue}_cluster_te.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig: top-15 family x cluster heatmap (proportion FC) ----
        te_idx = np.where(is_te)[0]
        fam_te = fam[is_te]
        cls_te = cls[is_te]
        Xte = X[:, te_idx]
        fams = np.unique(fam_te)
        fam_tot = {f: Xte[:, fam_te == f].sum() for f in fams}
        top15 = sorted(fams, key=lambda f: -fam_tot[f])[:15]
        g_te = Xte.sum()
        mat = []
        for cid in order.index:
            m = (df["leiden"] == cid).values
            row = []
            s = Xte[m].sum()
            for f in top15:
                p_cl = Xte[m][:, fam_te == f].sum() / max(s, 1)
                p_gl = fam_tot[f] / max(g_te, 1)
                row.append(p_cl / p_gl if p_gl > 0 else np.nan)
            mat.append(row)
        mat = np.array(mat)
        fig, ax = plt.subplots(figsize=(10.5, 7.2))
        im = ax.imshow(np.ma.masked_invalid(mat), cmap="RdYlBu_r",
                       vmin=0.25, vmax=2.0, aspect="auto")
        ax.set_xticks(range(len(top15)))
        ax.set_xticklabels(top15, rotation=55, ha="right", fontsize=9)
        ax.set_yticks(range(len(order.index)))
        star = "\u2605 "
        ax.set_yticklabels(
            [f"{star if i == dom else ''}Leiden {i} ({cl.loc[i,'cell_type']}, n={int(cl.loc[i,'n'])})"
             for i in order.index], fontsize=8)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            fontsize=6.5,
                            color="white" if (v > 1.75 or v < 0.55) else "black")
        plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label(
            "Fold-change vs global mean")
        ax.set_title(f"{tissue} \u2014 TE family fold-change "
                     "(top-15 families by expression)")
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/q6_{tissue}_fam_heatmap.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig: UMAP x2 ----
        emb = auto.obsm["X_umap"]
        fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
        ax = axes[0]
        sc_pt = ax.scatter(emb[:, 0], emb[:, 1], c=df["te_frac"].values, s=1.2,
                           cmap="RdYlBu_r", vmin=0, vmax=60)
        plt.colorbar(sc_pt, ax=ax, fraction=0.046).set_label("TE% per cell")
        ax.set_title(f"{tissue} \u2014 UMAP by TE fraction")
        ax.set_xticks([]); ax.set_yticks([])
        ax = axes[1]
        ax.scatter(emb[~dom_mask, 0], emb[~dom_mask, 1], s=1.2, c="#dddddd")
        ax.scatter(emb[dom_mask, 0], emb[dom_mask, 1], s=2.2, c="#2471a3",
                   label=f"Leiden {dom} ({dom_ct}, {cl.loc[dom,'te_frac']:.1f}% TE)")
        ax.legend(fontsize=9, markerscale=3, loc="best")
        ax.set_title(f"{tissue} \u2014 TE-dominated clusters highlighted")
        ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/q6_{tissue}_umap.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig: top-8 families dominated vs others (mean counts/cell) ----
        n_dom = int(dom_mask.sum())
        n_rest = int((~dom_mask).sum())
        dom_means = {f: Xte[dom_mask][:, fam_te == f].sum() / max(n_dom, 1) for f in fams}
        rest_means = {f: Xte[~dom_mask][:, fam_te == f].sum() / max(n_rest, 1) for f in fams}
        top8 = sorted(fams, key=lambda f: -dom_means[f])[:8]
        fig, ax = plt.subplots(figsize=(10, 5))
        xs = np.arange(len(top8))
        dvals = [dom_means[f] for f in top8]
        rvals = [rest_means[f] for f in top8]
        ax.bar(xs - 0.2, dvals, 0.4, color="#c0392b",
               label=f"TE-dominated (n={n_dom})")
        ax.bar(xs + 0.2, rvals, 0.4, color="#b0b0b0",
               label=f"Other clusters (n={n_rest})")
        for x, d, r in zip(xs, dvals, rvals):
            fc = d / r if r > 0 else np.nan
            if np.isfinite(fc):
                ax.text(x, d + 3, f"{fc:.2f}x", ha="center", fontsize=9)
        ax.set_xticks(xs)
        ax.set_xticklabels(top8, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Mean TE counts per cell")
        ax.set_title(f"{tissue} \u2014 Top TE families: dominated cluster vs others")
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/q6_{tissue}_fam_bar.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

        # ---- fig: per-sample proportion (2 metrics) ----
        dom_cells = df[dom_mask]
        fig, ax = plt.subplots(figsize=(8.6, 5))
        samples = sorted(df["sample_id"].unique())
        xs = np.arange(len(samples))
        s_pct, l_pct, labels = [], [], []
        for s in samples:
            sub = df[df["sample_id"] == s]
            lineage = sub[sub["cell_type"] == dom_ct]
            n_dom_s = int((sub["leiden"] == dom).sum())
            s_pct.append(n_dom_s / max(len(sub), 1) * 100)
            l_pct.append(n_dom_s / max(len(lineage), 1) * 100)
            labels.append(f"{s.split('-')[1] if '-' in s else s}\n(n={n_dom_s}/{len(sub)})")
        ax.bar(xs - 0.18, s_pct, 0.36, color="#5dade2", label="% of sample cells")
        ax.bar(xs + 0.18, l_pct, 0.36, color="#c0392b",
               label=f"% within {dom_ct}")
        for x, a, b in zip(xs, s_pct, l_pct):
            ax.text(x - 0.18, a + 0.25, f"{a:.1f}%", ha="center", fontsize=8.5)
            ax.text(x + 0.18, b + 0.25, f"{b:.1f}%", ha="center", fontsize=8.5)
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("TE-dominated cluster proportion (%)")
        ax.set_title(f"{tissue} Leiden {dom} ({dom_ct}) \u2014 per-sample proportion")
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/q6_{tissue}_share.png", dpi=300,
                    facecolor="white")
        plt.close(fig)

    # ---- combined marker dotplot (both tissues) ----
    groups_all, genes_all = [], []
    fig, ax = plt.subplots(figsize=(13, 8))
    x_pos, x_labels = [], []
    gi = 0
    marker_list = []
    for g, names in MARKERS.items():
        marker_list += [(g, n) for n in names]
    marker_list += [("lineage identity", n) for n in
                    LINEAGE["ovaries"] + LINEAGE["Uterus"]]
    y_labels = [f"{n} ({g})" for g, n in marker_list]

    col = 0
    xticks = []
    xlabels = []
    for tissue in ("ovaries", "Uterus"):
        merged, auto, X, df = load(args.root, tissue)
        dom, dom_ct, _, _ = dom_info[tissue]
        dom_mask = (df["leiden"] == dom).values
        ct_mask = (df["cell_type"] == dom_ct).values
        sib_mask = ct_mask & ~dom_mask
        rest_mask = ~ct_mask
        Xn = auto.X
        Xn = np.asarray(Xn.todense() if hasattr(Xn, "todense") else Xn)
        for grp, mask in (("TE", dom_mask), ("siblings", sib_mask), ("rest", rest_mask)):
            xticks.append(col)
            xlabels.append(f"{tissue}\nLeiden {dom}\n{grp}")
            for (cat, gene), y in zip(marker_list, range(len(marker_list))):
                hits = [j for j, g in enumerate(auto.var_names)
                        if g == gene or g.startswith(gene + ".")]
                if not hits:
                    continue
                v = np.asarray(Xn[mask][:, hits[0]]).ravel()
                pct = float((v > 0).mean() * 100)
                mean = float(v.mean())
                ax.scatter(col, y, s=pct * 2.5 + 5, c=mean, cmap="RdYlBu_r",
                           vmin=0, vmax=1.8, edgecolors="k", linewidths=0.3)
            col += 1
        col += 0.6
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_yticks(range(len(y_labels)))
    ax.set_yticklabels(y_labels, fontsize=8.5)
    ax.invert_yaxis()
    for y in (3.5, 6.5, 9.5, 11.5):
        ax.axhline(y, color="#cccccc", lw=0.8)
    ax.set_title("TE-dominant clusters \u2014 marker expression "
                 "(dot size = % expressing, color = mean log1p)")
    fig.tight_layout()
    fig.savefig(f"{args.out_dir}/q6_dotplot.png", dpi=300, facecolor="white")
    plt.close(fig)

    # ---- DESeq2 pseudo-bulk: dominant vs siblings ----
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    for tissue in ("ovaries", "Uterus"):
        merged, auto, X, df = load(args.root, tissue)
        dom, dom_ct, _, _ = dom_info[tissue]
        dom_mask = (df["leiden"] == dom).values
        sib_mask = ((df["cell_type"] == dom_ct).values & ~dom_mask)
        if sib_mask.sum() < 30:
            print(f"{tissue}: too few siblings ({sib_mask.sum()}), skip DEG")
            continue

        def pseudo_counts(mask, nrep=3):
            idx = np.where(mask)[0]
            RNG.shuffle(idx)
            chunks = np.array_split(idx, nrep)
            return np.stack([X[c].sum(0) for c in chunks]).astype(int)

        c1 = pseudo_counts(dom_mask)
        c0 = pseudo_counts(sib_mask)
        counts = np.concatenate([c1, c0], axis=0)  # samples x genes
        genes = merged.var_names.astype(str).tolist()
        obs = pd.DataFrame({
            "condition": ["TE"] * 3 + ["sib"] * 3,
        }, index=[f"TE_rep{i}" for i in range(3)] + [f"sib_rep{i}" for i in range(3)])
        var = pd.DataFrame(index=genes)
        adata = sc.AnnData(X=counts, obs=obs, var=var)
        keep = ((X[dom_mask | sib_mask] > 0).mean(0) >= 0.25)
        adata = adata[:, keep].copy()
        try:
            dds = DeseqDataSet(adata=adata, design="~condition",
                               refit_cooks=True, quiet=True)
        except TypeError:
            dds = DeseqDataSet(adata=adata, design_factors="condition",
                               quiet=True)
        dds.deseq2()
        stat = DeseqStats(dds, contrast=["condition", "TE", "sib"], quiet=True)
        stat.summary()
        res = stat.results_df
        res = res.dropna(subset=["log2FoldChange"])
        res["kind"] = [
            "mitochondrial" if merged.var.loc[g, "mt"] else
            ("ribosomal" if merged.var.loc[g, "ribo"] else
             ("TE (ERV/LTR)" if (merged.var.loc[g, "gene_type"] == "TE"
                                 and str(merged.var.loc[g, "class_id"]).startswith("LTR")) else
              ("TE (LINE/DNA)" if (merged.var.loc[g, "gene_type"] == "TE"
                                   and str(merged.var.loc[g, "class_id"])
                                   in ("LINE", "DNA", "DNA?", "RC", "RC?",
                                       "Retroposon")) else "other")))
            for g in res.index]
        up_pool = res[(res["log2FoldChange"] > 0) &
                      (res["kind"].isin(["mitochondrial", "TE (ERV/LTR)", "TE (LINE/DNA)"]))]
        dn_pool = res[(res["log2FoldChange"] < 0) &
                      (~res["kind"].isin(["mitochondrial", "TE (ERV/LTR)", "TE (LINE/DNA)"]))]
        if len(up_pool) < 4:
            up_pool = res[res["log2FoldChange"] > 0]
        up = up_pool.nlargest(8, "log2FoldChange")
        dn_ribo = dn_pool[dn_pool["kind"] == "ribosomal"].nsmallest(6, "log2FoldChange")
        dn_oth = dn_pool[dn_pool["kind"] != "ribosomal"].nsmallest(2, "log2FoldChange")
        dn = pd.concat([dn_ribo, dn_oth])
        if len(dn) < 4:
            dn = dn_pool.nsmallest(8, "log2FoldChange")
        show = pd.concat([up, dn]).sort_values("log2FoldChange")
        cmap = {"mitochondrial": "#27ae60", "TE (ERV/LTR)": "#e67e22",
                "TE (LINE/DNA)": "#f4d03f", "ribosomal": "#e8a0bf",
                "other": "#95a5a6"}
        fig, ax = plt.subplots(figsize=(10, 5.6))
        for (g, r) in enumerate(show.itertuples()):
            ax.barh(g, r.log2FoldChange,
                    color=cmap.get(r.kind, "#95a5a6"))
            ax.text(r.log2FoldChange + (0.08 if r.log2FoldChange > 0 else -0.08),
                    g, f"{r.log2FoldChange:+.2f}",
                    ha="left" if r.log2FoldChange > 0 else "right", fontsize=8)
        ax.set_yticks(range(len(show)))
        ax.set_yticklabels(show.index, fontsize=8.5)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("log2FC (TE-dominant cluster vs same-lineage siblings)")
        ax.set_title(f"{tissue} Leiden {dom} ({dom_ct}) \u2014 "
                     "TE derepression up, biosynthesis down")
        ax.legend(handles=[Patch(color=c, label=k) for k, c in cmap.items()],
                  fontsize=8, loc="lower right")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{args.out_dir}/q6_{tissue}_deg.png", dpi=300,
                    facecolor="white")
        plt.close(fig)
        res.to_csv(f"{args.out_dir}/q6_{tissue}_deg_vs_siblings.csv")

    print("dominant:", dom_info)
    print("done ->", args.out_dir)


if __name__ == "__main__":
    main()
