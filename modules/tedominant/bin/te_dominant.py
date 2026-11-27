#!/usr/bin/env python3
"""TE-dominant cluster identification module engine (v1/v2).

Per tissue, load the scTE-merged (gene+TE features) and scanpy-auto
(annotations) h5ad pair, identify TE-dominant Leiden clusters with the
v2 anti-denominator-collapse gates, and write:

  te_dominant_clusters.csv      per-cluster metric table + v2 verdict
  cluster_burden_<tissue>.png   per-cluster absolute TE burden bars
  umap_te_<tissue>.png          UMAP colored by absolute TE counts + dominant
  markers_dotplot.png           senescence/lineage marker dotplot
  deg_<tissue>.png              pseudo-bulk DEG volcano (dominant vs siblings)

v2 gates (see README §2.2): a cluster is ELIGIBLE as TE-dominant only if
its median absolute TE counts >= 0.5x same-lineage siblings median AND
>= 0.5x tissue median AND n >= min_cells; among eligible clusters the
one with max mean TE% wins. Rejections are logged with the gate numbers.
"""
import argparse
import os
import warnings
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

MARKERS = {
    "SASP / inflammation": ["IL6", "CXCL8", "CXCL1", "CCL2", "IL1A", "IL1B"],
    "senescence / CDKi": ["CDKN1A", "CDKN2B", "CDKN1B", "RB1",
                          "TP53", "GADD45A", "GADD45B", "SERPINE1"],
    "DNA damage / aging": ["LMNB1", "GLB1", "HMGB1", "IGFBP3",
                           "IGFBP7", "FOXO1", "FOXO3", "SIRT1", "SIRT6"],
    "proliferation": ["MKI67", "TOP2A", "CDK1"],
    "stress / apoptosis": ["BAX", "CASP3"],
}
LINEAGE = {
    "granulosa lineage": ["AMH", "CYP19A1", "FSHR", "INHA", "GJA1"],
    "glandular epithelial lineage": ["EPCAM", "KRT19", "MUC1", "PAEP", "FOXA2"],
}


def load(root, tissue):
    """Load the h5ad pair from {root}/{tissue}/{tissue}_scTE_{merged,auto}.h5ad."""
    merged = sc.read_h5ad(f"{root}/{tissue}/{tissue}_scTE_merged.h5ad")
    auto = sc.read_h5ad(f"{root}/{tissue}/{tissue}_scTE_auto.h5ad")
    common = merged.obs_names.intersection(auto.obs_names)
    merged, auto = merged[common].copy(), auto[common].copy()
    return merged, auto


def cell_metrics(merged):
    """Per-cell TE metrics from the integer counts layer."""
    X = merged.layers["counts"] if "counts" in merged.layers else merged.X
    X = np.asarray(X.todense() if hasattr(X, "todense") else X)
    is_te = merged.var["gene_type"].astype(str).values == "TE"
    te = X[:, is_te].sum(1).astype(float)
    tot = X.sum(1).astype(float)
    gene = tot - te
    return {
        "te_abs": te,
        "te_pct": te / np.maximum(tot, 1.0) * 100.0,
        "gene_abs": gene,
        "n_gene": (X[:, ~is_te] > 0).sum(1).astype(float),
    }


def dominant_v2(cm, leiden, cell_type, min_cells, gate_frac, method="v1"):
    """TE-dominant selection; --method v1 (default) = pure max mean TE%
    (legacy), v2 = eligibility gates by absolute TE first, then max TE%.

    Returns (dom, dom_ct, table_df, rejections) where table_df carries every
    per-cluster metric and the eligible/verdict columns, and rejections is a
    list of human-readable gate messages (empty for v1).
    """
    clusters = sorted(set(leiden))
    global_med = float(np.median(cm["te_abs"]))
    rows, rejects = [], []
    for c in clusters:
        m = leiden == c
        ct = Counter(cell_type[m]).most_common(1)[0][0]
        sib = cm["te_abs"][(cell_type == ct) & (leiden != c)]
        sib_med = float(np.median(sib)) if len(sib) >= 30 else float("inf")
        te_abs_med = float(np.median(cm["te_abs"][m]))
        gate_lineage = te_abs_med >= gate_frac * sib_med
        gate_global = te_abs_med >= gate_frac * global_med
        gate_n = int(m.sum()) >= min_cells
        eligible = gate_lineage and gate_global and gate_n
        rows.append({
            "cluster": c, "cell_type": ct, "n_cells": int(m.sum()),
            "te_pct_mean": float(cm["te_pct"][m].mean()),
            "te_abs_median": te_abs_med,
            "gene_abs_median": float(np.median(cm["gene_abs"][m])),
            "n_gene_median": float(np.median(cm["n_gene"][m])),
            "sib_te_abs_median": sib_med,
            "tissue_te_abs_median": global_med,
            "gate_lineage": gate_lineage, "gate_global": gate_global,
            "gate_n": gate_n, "eligible": eligible,
        })
    tab = pd.DataFrame(rows).set_index("cluster")
    old_dom = tab["te_pct_mean"].idxmax()
    if method == "v1":
        return old_dom, tab.loc[old_dom, "cell_type"], tab, rejects
    elig = tab[tab["eligible"]]
    dom = elig["te_pct_mean"].idxmax() if len(elig) else old_dom
    if not len(elig):
        rejects.append(f"WARNING: no eligible cluster; falling back to "
                       f"pure-TE% winner Leiden {old_dom}")
    elif str(old_dom) != str(dom):
        r = tab.loc[old_dom]
        rejects.append(
            f"Leiden {old_dom} rejected: TE% {r.te_pct_mean:.1f} ranked #1 but "
            f"TE abs median {r.te_abs_median:.0f} < {gate_frac}x siblings "
            f"{r.sib_te_abs_median:.0f} (denominator-collapse guard)")
    return dom, tab.loc[dom, "cell_type"], tab, rejects


def fig_burden(tab, dom, tissue, gate_frac, outdir, method="v1"):
    if method == "v1":
        order = tab.sort_values("te_pct_mean", ascending=False)
        fig, ax = plt.subplots(figsize=(11, 5))
        colors = ["#c0392b" if i == dom else "#b0b0b0" for i in order.index]
        ax.bar(range(len(order)), order["te_pct_mean"], color=colors)
        gm = float(tab["te_pct_mean"].mean())
        ax.axhline(gm, ls="--", lw=1.2, color="#333333",
                   label=f"global mean {gm:.1f}%")
        ax.set_ylabel("TE fraction (%)")
        ax.set_title(f"{tissue} \u2014 Per-cluster TE fraction (raw counts)")
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{outdir}/cluster_burden_{tissue}.png", dpi=300,
                    facecolor="white")
        plt.close(fig)
        return
    order = tab.sort_values("te_abs_median", ascending=False)
    fig, ax = plt.subplots(figsize=(11, 5))
    colors = ["#c0392b" if i == dom else "#b0b0b0" for i in order.index]
    ax.bar(range(len(order)), order["te_abs_median"], color=colors)
    gate = gate_frac * float(tab["tissue_te_abs_median"].iloc[0])
    ax.axhline(gate, ls="--", lw=1.2, color="#333333",
               label=f"v2 gate: {gate_frac}x tissue median ({gate:.0f})")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{i} ({tab.loc[i,'cell_type']}, n={int(tab.loc[i,'n_cells'])})"
                        for i in order.index], rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Median TE counts per cell")
    ax.set_title(f"{tissue} — Per-cluster TE burden (absolute; v2-gated dominant in red)")
    ax.legend(fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{outdir}/cluster_burden_{tissue}.png", dpi=300, facecolor="white")
    plt.close(fig)


def fig_umap(auto, cm, dom, dom_ct, te_pct_dom, tissue, outdir, method="v1"):
    emb = auto.obsm["X_umap"]
    dom_mask = (auto.obs["leiden"].astype(str).values == str(dom))
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
    if method == "v1":
        sc_pt = axes[0].scatter(emb[:, 0], emb[:, 1], c=cm["te_pct"], s=1.2,
                                cmap="RdYlBu_r", vmin=0, vmax=60)
        plt.colorbar(sc_pt, ax=axes[0], fraction=0.046).set_label("TE% per cell")
        axes[0].set_title(f"{tissue} \u2014 UMAP by TE fraction")
    else:
        vmax = float(np.percentile(cm["te_abs"], 98))
        sc_pt = axes[0].scatter(emb[:, 0], emb[:, 1], c=cm["te_abs"], s=1.2,
                                cmap="RdYlBu_r", vmin=0, vmax=vmax)
        plt.colorbar(sc_pt, ax=axes[0], fraction=0.046).set_label("TE counts per cell")
        axes[0].set_title(f"{tissue} — UMAP by TE counts (absolute, v2 metric)")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    axes[1].scatter(emb[~dom_mask, 0], emb[~dom_mask, 1], s=1.2, c="#dddddd")
    axes[1].scatter(emb[dom_mask, 0], emb[dom_mask, 1], s=2.2, c="#2471a3",
                    label=f"Leiden {dom} ({dom_ct}, {te_pct_dom:.1f}% TE)")
    axes[1].set_title(f"{tissue} — TE-dominant cluster highlighted")
    axes[1].legend(fontsize=9, markerscale=3, loc="best")
    axes[1].set_xticks([]); axes[1].set_yticks([])
    fig.tight_layout()
    fig.savefig(f"{outdir}/umap_te_{tissue}.png", dpi=300, facecolor="white")
    plt.close(fig)


def fig_dotplot(root, tissues_dom, outdir, dot_args):
    """markers x (TE / siblings / rest) per tissue; color = mean log1p of
    auto.raw.X (full gene set), size = % expressing."""
    fig, ax = plt.subplots(figsize=(13, 11))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    ax.spines[:].set_visible(False); ax.tick_params(length=0)
    marker_list = [(g, n) for g, names in MARKERS.items() for n in names]
    marker_list += [(cat, n) for cat, genes in LINEAGE.items() for n in genes]
    cat_bounds, seen, y_labels = [], set(), []
    for i, (g, _) in enumerate(marker_list):
        if g not in seen:
            cat_bounds.append(i); seen.add(g)
    y_labels = [f"{n} ({g})" for g, n in marker_list]
    xticks, xlabels, col, missing = [], [], 0, set()
    for tissue in tissues_dom:
        merged, auto = load(root, tissue)
        dom, dom_ct = tissues_dom[tissue]
        leiden = auto.obs["leiden"].astype(str).values
        ct = auto.obs["cell_type"].astype(str).values
        dom_mask, ct_mask = leiden == str(dom), ct == dom_ct
        sib_mask, rest_mask = ct_mask & ~dom_mask, ~ct_mask
        ref = auto.raw if auto.raw is not None else auto
        name2idx = {}
        for j, g in enumerate(ref.var_names):
            name2idx.setdefault(g, j)
            if "." in g:
                name2idx.setdefault(g.split(".")[0], j)
        Xn = ref.X
        for grp, mask in (("TE", dom_mask), ("siblings", sib_mask), ("rest", rest_mask)):
            xticks.append(col)
            xlabels.append(f"{tissue}\nLeiden {dom}\n{grp}")
            for (cat, gene), y in zip(marker_list, range(len(marker_list))):
                j = name2idx.get(gene)
                if j is None:
                    missing.add(gene); continue
                v = Xn[mask][:, j]
                v = np.asarray(v.todense() if hasattr(v, "todense") else v).ravel()
                pct, mean = float((v > 0).mean() * 100), float(v.mean())
                ax.scatter(col, y, s=pct * 2.5 + 5, c=mean, cmap="RdYlBu_r",
                           vmin=0, vmax=1.8, edgecolors="k", linewidths=0.3)
            col += 1
        col += 0.6
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_yticks(range(len(y_labels))); ax.set_yticklabels(y_labels, fontsize=8.5)
    ax.invert_yaxis()
    for b in cat_bounds[1:]:
        ax.axhline(b - 0.5, color="#cccccc", lw=0.8)
    fig.tight_layout(rect=(0, 0.02, 0.86, 1))
    import matplotlib.cm as mcm
    from matplotlib.colors import Normalize
    sm = mcm.ScalarMappable(cmap="RdYlBu_r", norm=Normalize(0, 1.8))
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("mean log1p", fontsize=9)
    size_handles = [ax.scatter([], [], s=p * 2.5 + 5, c="#999999", edgecolors="k",
                               linewidths=0.3, label=f"{p}%")
                    for p in (25, 50, 75, 100)]
    ax.legend(handles=size_handles, title="% expressing", loc="upper left",
              bbox_to_anchor=(1.10, 1.0), frameon=False, fontsize=9,
              title_fontsize=9, labelspacing=1.2, borderpad=1.0, handletextpad=1.6)
    fig.savefig(f"{outdir}/markers_dotplot.png", dpi=300, facecolor="white")
    plt.close(fig)
    if missing:
        print("genes not in annotation (skipped):", ", ".join(sorted(missing)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True,
                    help="combine stage dir holding {tissue}_scTE_{merged,auto}.h5ad")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tissues", default="ovaries,Uterus",
                    help="comma list (default ovaries,Uterus)")
    ap.add_argument("--method", choices=["v1", "v2"], default="v1",
                    help="dominant-selection method: v1 = pure max mean TE%% "
                         "(legacy, default); v2 = absolute-TE eligibility "
                         "gates then max TE%% (anti denominator-collapse)")
    ap.add_argument("--min-cells", type=int, default=100,
                    help="cluster size gate (default 100)")
    ap.add_argument("--gate-frac", type=float, default=0.5,
                    help="eligibility = TE abs median >= gate_frac x reference "
                         "median, for both siblings and tissue (default 0.5)")
    ap.add_argument("--skip-dotplot", action="store_true",
                    help="skip the marker dotplot figure")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tissues = [t.strip() for t in args.tissues.split(",") if t.strip()]

    summary, tables, dom_map = [], [], {}
    for tissue in tissues:
        merged, auto = load(args.indir, tissue)
        cm = cell_metrics(merged)
        leiden = auto.obs["leiden"].astype(str).values
        cell_type = auto.obs["cell_type"].astype(str).values
        dom, dom_ct, tab, rejects = dominant_v2(cm, leiden, cell_type,
                                                args.min_cells, args.gate_frac,
                                                method=args.method)
        for msg in rejects:
            print(f"[gate] {tissue}: {msg}")
        dom_map[tissue] = (dom, dom_ct)
        te_pct_dom = float(tab.loc[dom, "te_pct_mean"])
        print(f"{tissue}: dominant = Leiden {dom} ({dom_ct}), "
              f"mean TE% = {te_pct_dom:.2f}, n = {int(tab.loc[dom, 'n_cells'])}")

        tab_out = tab.copy()
        tab_out.insert(0, "tissue", tissue)
        tab_out["dominant"] = tab_out.index.astype(str) == str(dom)
        tables.append(tab_out)

        fig_burden(tab, dom, tissue, args.gate_frac, args.outdir,
                   method=args.method)
        fig_umap(auto, cm, dom, dom_ct, te_pct_dom, tissue, args.outdir,
                 method=args.method)

    pd.concat(tables).to_csv(f"{args.outdir}/te_dominant_clusters.csv")
    print(f"table -> {args.outdir}/te_dominant_clusters.csv")
    if not args.skip_dotplot:
        fig_dotplot(args.indir, dom_map, args.outdir, args)


if __name__ == "__main__":
    main()
