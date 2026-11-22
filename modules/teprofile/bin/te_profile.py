#!/usr/bin/env python3
"""TE profile module engine: three-layer TE-expression characterization.

Layer 1 (baseline): per-CELL-TYPE intrinsic TE features — TE fraction
distribution, depth-normalized absolute TE output, TE family composition
(Alu/L1/L2/MIR/ERV...), gene-program health (detected genes). Answers
"what does each cell type's TE expression normally look like".

Layer 2 (outlier subpopulations): per-cluster deviation from the SAME
cell-type background — fold-enrichment of absolute TE output, TE fraction
shift, family-composition shift (Aitchison distance vs a bootstrap
null of same-size lineage subsets — adaptive, size-aware). Reports EVERY cluster's deviation profile (no winner
election; "no outlier" is a valid outcome). Collapsed clusters (low gene
output, normal TE output) are labelled as such rather than mislabelled
"TE-dominant".

Layer 3 (remodeling): per cell type, young-vs-aged shift of the layer-1
features (baseline up-shift) AND per-cluster size shift (subpopulation
expansion), computed from a sample-level age/design annotation when
available; skipped (with a notice) when no sample grouping is configured.

Outputs (outdir):
  te_profile_baseline.csv        per (tissue, cell_type): feature summaries
  te_profile_clusters.csv        per (tissue, cluster): deviation profile
  te_profile_remodeling.csv      per (tissue, cell_type): layer-3 shifts
  baseline_te_fraction_<t>.png   violin: TE fraction by cell type
  baseline_te_output_<t>.png     violin: depth-normalized TE output by type
  baseline_family_comp_<t>.png   stacked bars: family composition per type
  cluster_deviation_<t>.png      per-cluster deviation panel (3 metrics)
  remodeling_shift_<t>.png       layer-3 shift summary (when grouping set)

Housekeeping normalization: absolute TE output is normalized per cell by a
housekeeping gene set (HK flag; default the intersection of a curated
human/macaque list with the data's gene symbols). te_norm = TE_counts /
max(sum(HK counts), eps). This kills sequencing-depth mixing WITHOUT the
denominator collapse of te_pct (HK genes are not the whole transcriptome,
so gene-program collapse does not deflate the normalizer as strongly).
"""
import argparse
import os
import warnings
from collections import Counter, OrderedDict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

# Curated housekeeping genes (human symbols; macaque orthologs share
# symbols in the merged object's var_names — intersection decides).
HK_CANDIDATES = [
    "ACTB", "GAPDH", "RPL13A", "RPLP0", "HPRT1", "TBP", "PPIA", "B2M",
    "SDHA", "YWHAZ", "UBC", "PGK1", "RPL4", "RPS18", "RPLP1", "FTL",
]


def load(root, tissue):
    merged = sc.read_h5ad(f"{root}/{tissue}/{tissue}_scTE_merged.h5ad")
    auto = sc.read_h5ad(f"{root}/{tissue}/{tissue}_scTE_auto.h5ad")
    common = merged.obs_names.intersection(auto.obs_names)
    merged, auto = merged[common].copy(), auto[common].copy()
    return merged, auto


def cell_metrics(merged):
    X = merged.layers["counts"] if "counts" in merged.layers else merged.X
    X = np.asarray(X.todense() if hasattr(X, "todense") else X)
    is_te = merged.var["gene_type"].astype(str).values == "TE"
    var_names = merged.var_names.astype(str).values
    # housekeeping normalizer: intersection of candidates with gene symbols
    hk_idx = [i for i, g in enumerate(var_names) if (not is_te[i]) and
              any(g == h or g.startswith(h + ".") for h in HK_CANDIDATES)]
    if len(hk_idx) < 5:  # degenerate HK set -> fall back to total gene output
        hk_idx = np.where(~is_te)[0]
        hk_mode = "all_genes"
    else:
        hk_mode = "housekeeping"
    te = X[:, is_te].sum(1).astype(float)
    tot = X.sum(1).astype(float)
    hk = X[:, hk_idx].sum(1).astype(float)
    eps = 1.0
    return {
        "te_abs": te,
        "te_pct": te / np.maximum(tot, 1.0) * 100.0,
        "te_norm": te / np.maximum(hk, eps) * 100.0,   # TE per 100 HK counts
        "gene_abs": tot - te,
        "n_gene": (X[:, ~is_te] > 0).sum(1).astype(float),
        "hk_mode": hk_mode,
        "is_te": is_te,
        "X": X,
    }


def family_shares(cm, merged, mask):
    """CLR-transformed TE family shares for cells in mask (per-cluster use).

    Returns (families, mean_share_df) — mean per-family share of TE output.
    """
    fam = merged.var["family_id"].astype(str).values if "family_id" in merged.var \
        else np.array(["unknown"] * merged.n_vars)
    te_cols = np.where(cm["is_te"])[0]
    fam_of = {j: fam[j] for j in te_cols}
    X = cm["X"][mask][:, te_cols]
    # aggregate to family-level counts per cell
    fams = sorted(set(fam_of.values()))
    fam_idx = {f: i for i, f in enumerate(fams)}
    agg = np.zeros((X.shape[0], len(fams)))
    for j, col in enumerate(te_cols):
        agg[:, fam_idx[fam_of[col]]] += X[:, j]
    share = agg / np.maximum(agg.sum(1, keepdims=True), 1e-9)
    return fams, pd.DataFrame(share, columns=fams)


def clr(x, eps=1e-6):
    x = np.asarray(x, dtype=float) + eps
    return np.log(x / np.exp(np.log(x).mean(axis=-1, keepdims=True)))


def aitchison(a, b):
    """Aitchison distance between two compositional vectors."""
    return float(np.linalg.norm(clr(a) - clr(b)))


def layer1_baseline(merged, auto, cm, tissue):
    ct = auto.obs["cell_type"].astype(str).values
    rows = []
    for t in sorted(set(ct)):
        m = ct == t
        te_fam_shares = family_shares(cm, merged, m)[1].mean().values
        rows.append({
            "tissue": tissue, "cell_type": t, "n_cells": int(m.sum()),
            "te_pct_median": float(np.median(cm["te_pct"][m])),
            "te_pct_q25": float(np.percentile(cm["te_pct"][m], 25)),
            "te_pct_q75": float(np.percentile(cm["te_pct"][m], 75)),
            "te_norm_median": float(np.median(cm["te_norm"][m])),
            "te_norm_q25": float(np.percentile(cm["te_norm"][m], 25)),
            "te_norm_q75": float(np.percentile(cm["te_norm"][m], 75)),
            "n_gene_median": float(np.median(cm["n_gene"][m])),
            "top_family": family_shares(cm, merged, m)[1].mean().idxmax(),
        })
    return pd.DataFrame(rows)


def _share_mean(cm, merged, idx):
    """Mean family-share vector for cells at integer indices idx."""
    _, df = family_shares(cm, merged, idx)
    return df.mean().values


def layer2_clusters(merged, auto, cm, tissue, fe_cut, n_boot, rng_seed):
    """Per-cluster deviation with ADAPTIVE composition threshold.

    The Aitchison distance has no universal absolute scale (it grows with
    family count, depth, and SHRINKS with cluster size via sampling noise),
    so a fixed cutoff is meaningless (real data: every cluster sits >2.0).
    Instead each cluster is compared against a BOOTSTRAP NULL of same-size
    random subsets drawn from its own lineage background: the observed
    distance must exceed the null's quantile (comp_q) to count as a
    composition shift. This auto-calibrates for cluster size — a small
    cluster needs a larger distance to stand out.

    fe display for collapsed clusters: when the gene program collapses the
    HK normalizer itself deflates (HK genes drop too), inflating te_norm;
    the honest display for collapsed rows is fe computed on RAW te counts
    ratio (te_abs_median ratio), labeled fe_raw. Non-collapsed rows keep
    the HK-normalized fe.
    """
    leiden = auto.obs["leiden"].astype(str).values
    ct = auto.obs["cell_type"].astype(str).values
    rng = np.random.default_rng(rng_seed)
    rows = []
    for c in sorted(set(leiden)):
        m = leiden == c
        ct_c = Counter(ct[m]).most_common(1)[0][0]
        bg = (ct == ct_c) & ~m           # same-lineage background
        if bg.sum() < 30:
            bg = ~m                       # fall back to whole tissue
        gene_ratio = float(np.median(cm["n_gene"][m]) / max(np.median(cm["n_gene"][bg]), 1e-9))
        collapsed = gene_ratio < 0.7
        # fe: HK-normalized for healthy clusters; raw-counts ratio when
        # collapsed (HK deflation would otherwise inflate it)
        fe_hk = float(np.median(cm["te_norm"][m]) / max(np.median(cm["te_norm"][bg]), 1e-9))
        fe_raw = float(np.median(cm["te_abs"][m]) / max(np.median(cm["te_abs"][bg]), 1e-9))
        fe = fe_raw if collapsed else fe_hk
        pct_shift = float(np.median(cm["te_pct"][m]) - np.median(cm["te_pct"][bg]))
        # observed composition distance
        bg_idx = np.where(bg)[0]
        c_idx = np.where(m)[0]
        s_bg = _share_mean(cm, merged, bg_idx)
        aD = aitchison(_share_mean(cm, merged, c_idx), s_bg)
        # bootstrap null: same-size random subsets of the background
        n_c = len(c_idx)
        null = np.empty(n_boot)
        for b in range(n_boot):
            sub = rng.choice(bg_idx, size=min(n_c, len(bg_idx)), replace=False)
            null[b] = aitchison(_share_mean(cm, merged, sub), s_bg)
        aD_null_q = float(np.quantile(null, 0.95))
        comp_outlier = bool(aD > aD_null_q)
        te_high = (fe >= fe_cut) and not collapsed
        label = ("collapsed" if collapsed else
                 ("te_outlier" if (te_high and comp_outlier) else
                  ("te_enriched" if te_high else "typical")))
        rows.append({
            "tissue": tissue, "cluster": c, "cell_type": ct_c,
            "n_cells": int(m.sum()),
            "te_norm_median": float(np.median(cm["te_norm"][m])),
            "fe_vs_lineage": fe, "fe_type": "raw" if collapsed else "hk_norm",
            "te_pct_shift": pct_shift,
            "n_gene_ratio": gene_ratio, "family_aitchison": aD,
            "aitchison_null_q95": aD_null_q, "comp_outlier": comp_outlier,
            "collapsed": collapsed, "label": label,
        })
    return pd.DataFrame(rows)


def layer3_remodeling(merged, auto, cm, tissue, group_by):
    """Per cell type: feature shift between sample groups + the dominant
    TE outlier cluster's size shift. group_by maps sample_id -> group."""
    sid = auto.obs["sample_id"].astype(str).values
    ct = auto.obs["cell_type"].astype(str).values
    groups = sorted(set(group_by.values()))
    if len(groups) < 2:
        return None
    g0, g1 = groups[0], groups[-1]        # first vs last (e.g. young vs aged)
    rows = []
    for t in sorted(set(ct)):
        m_t = ct == t
        m0 = m_t & np.isin(sid, [s for s, g in group_by.items() if g == g0])
        m1 = m_t & np.isin(sid, [s for s, g in group_by.items() if g == g1])
        if m0.sum() < 20 or m1.sum() < 20:
            continue
        rows.append({
            "tissue": tissue, "cell_type": t,
            f"te_pct_median[{g0}]": float(np.median(cm["te_pct"][m0])),
            f"te_pct_median[{g1}]": float(np.median(cm["te_pct"][m1])),
            "te_pct_shift": float(np.median(cm["te_pct"][m1]) - np.median(cm["te_pct"][m0])),
            f"te_norm_median[{g0}]": float(np.median(cm["te_norm"][m0])),
            f"te_norm_median[{g1}]": float(np.median(cm["te_norm"][m1])),
            "te_norm_shift": float(np.median(cm["te_norm"][m1]) - np.median(cm["te_norm"][m0])),
            f"n_cells[{g0}]": int(m0.sum()), f"n_cells[{g1}]": int(m1.sum()),
        })
    return pd.DataFrame(rows)


# ------------------------- figures -------------------------

def _violin_by_type(ax, values, ct, title, ylabel):
    order = sorted(set(ct), key=lambda t: -np.median(values[ct == t]))
    data = [values[ct == t] for t in order]
    vp = ax.violinplot(data, showmedians=True, widths=0.8)
    for b in vp["bodies"]:
        b.set_facecolor("#7fa8c9"); b.set_alpha(0.7)
    for key in ("cmedians", "cbars"):
        vp[key].set_color("#2c3e50")
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([f"{t}\n(n={int((ct==t).sum())})" for t in order],
                       rotation=60, ha="right", fontsize=7)
    ax.set_title(title, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)


def fig_baseline(cm, auto, tissue, outdir, hk_mode):
    ct = auto.obs["cell_type"].astype(str).values
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    _violin_by_type(axes[0], cm["te_pct"], ct,
                    f"{tissue} — TE fraction by cell type",
                    "TE fraction (%)")
    safe = np.maximum(cm["te_norm"], 1e-2)  # log-safe floor
    _violin_by_type(axes[1], safe, ct,
                    f"{tissue} — TE output per 100 {'HK' if hk_mode=='housekeeping' else 'gene'} counts (log)",
                    "TE counts / 100 HK counts (log)")
    axes[1].set_yscale("log")
    fig.tight_layout()
    fig.savefig(f"{outdir}/baseline_te_fraction_{tissue}.png", dpi=300,
                facecolor="white")
    plt.close(fig)


def fig_family_comp(cm, merged, auto, tissue, outdir):
    ct = auto.obs["cell_type"].astype(str).values
    types = sorted(set(ct), key=lambda t: -np.median(cm["te_pct"][ct == t]))
    fams, shares_by_type = None, {}
    for t in types:
        fams, df = family_shares(cm, merged, ct == t)
        shares_by_type[t] = df.mean().values
    M = np.vstack([shares_by_type[t] for t in types])
    # top-8 families, rest as 'other'
    top = np.argsort(-M.mean(0))[:8]
    keep = M[:, top]
    other = M.sum(1) - keep.sum(1)
    M2 = np.hstack([keep, other[:, None]])
    labels = [fams[i] for i in top] + ["other"]
    fig, ax = plt.subplots(figsize=(10, max(4, 0.5 * len(types))))
    left = np.zeros(len(types))
    cmap = plt.cm.tab20
    for j, lab in enumerate(labels):
        ax.barh(range(len(types)), M2[:, j], left=left, label=lab,
                color=cmap(j % 20))
        left += M2[:, j]
    ax.set_yticks(range(len(types)))
    ax.set_yticklabels([f"{t} (n={int((ct==t).sum())})" for t in types], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("share of TE output")
    ax.set_title(f"{tissue} — TE family composition by cell type")
    ax.legend(fontsize=7, ncol=1, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), frameon=False,
              title="TE family", title_fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{outdir}/baseline_family_comp_{tissue}.png", dpi=300,
                facecolor="white", bbox_inches="tight")
    plt.close(fig)


def fig_cluster_dev(tab, tissue, outdir, fe_cut):
    order = tab.sort_values("fe_vs_lineage", ascending=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    y = np.arange(len(order))
    labels = [f"{r.cluster} ({r.cell_type[:16]}, n={r.n_cells}"
              f"{', fe:raw' if r.fe_type == 'raw' else ''})"
              for r in order.itertuples()]
    color_map = {"te_outlier": "#c0392b", "te_enriched": "#e67e22",
                 "collapsed": "#7f8c8d", "typical": "#bdc3c7"}
    colors = [color_map.get(r.label, "#bdc3c7") for r in order.itertuples()]
    axes[0].barh(y, order["fe_vs_lineage"], color=colors)
    axes[0].axvline(fe_cut, ls="--", lw=1, color="#333")
    axes[0].set_xlabel("TE output fold-enrichment vs lineage")
    axes[1].barh(y, order["family_aitchison"], color=colors)
    # per-cluster adaptive threshold: bootstrap null q95 (tick markers);
    # a bar is a composition outlier only if it passes its own marker
    axes[1].scatter(order["aitchison_null_q95"], y, marker="|", s=200,
                    color="#333333", linewidths=1.2, zorder=5,
                    label="null q95 (adaptive)")
    axes[1].legend(fontsize=7, loc="lower right")
    axes[1].set_xlabel("family composition shift (Aitchison; | = per-cluster null q95)")
    axes[2].barh(y, order["n_gene_ratio"], color=colors)
    axes[2].axvline(0.7, ls="--", lw=1, color="#333")
    axes[2].set_xlabel("detected-gene ratio vs lineage")
    for ax in axes:
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6)
        ax.invert_yaxis()
        ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in color_map.values()]
    fig.legend(handles, color_map.keys(), fontsize=8, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, -0.02), frameon=False)
    fig.suptitle(f"{tissue} — cluster deviation from lineage background", y=1.02)
    fig.tight_layout()
    fig.savefig(f"{outdir}/cluster_deviation_{tissue}.png", dpi=300,
                facecolor="white", bbox_inches="tight")
    plt.close(fig)


def fig_remodeling(tab, tissue, outdir, groups):
    if tab is None or not len(tab):
        # 占位图:layer 3 在该组织无可比较类型(分组样本不在该组织/细胞数不足)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, f"{tissue}: no cell type comparable between groups\n"
                "(samples missing from groups or <20 cells per group)",
                ha="center", va="center", fontsize=10)
        ax.axis("off")
        fig.savefig(f"{outdir}/remodeling_shift_{tissue}.png", dpi=150,
                    facecolor="white")
        plt.close(fig)
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    y = np.arange(len(tab))
    n0 = tab[f"n_cells[{groups[0]}]"].values
    n1 = tab[f"n_cells[{groups[-1]}]"].values
    labels = [f"{ct} (n={a}→{b})" for ct, a, b in
              zip(tab["cell_type"].values, n0, n1)]
    ax.barh(y - 0.2, tab["te_norm_shift"], height=0.4,
            color="#c0392b", label="Δ TE output (norm)")
    ax.barh(y + 0.2, tab["te_pct_shift"], height=0.4,
            color="#7fa8c9", label="Δ TE fraction")
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False)
    ax.set_title(f"{tissue} — {groups[0]} → {groups[-1]} shift per cell type")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{outdir}/remodeling_shift_{tissue}.png", dpi=300,
                facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tissues", default="ovaries,Uterus")
    ap.add_argument("--fe-cut", type=float, default=1.5,
                    help="layer-2 TE-enrichment cutoff (fold vs lineage; default 1.5)")
    ap.add_argument("--n-boot", type=int, default=200,
                    help="bootstrap resamples for the composition null (default 200)")
    ap.add_argument("--seed", type=int, default=0,
                    help="bootstrap RNG seed (default 0)")
    ap.add_argument("--groups", default="",
                    help="sample grouping for layer 3: 'young:s1,s2;aged:s3,s4'")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tissues = [t.strip() for t in args.tissues.split(",") if t.strip()]
    group_by = {}
    if args.groups:
        for part in args.groups.split(";"):
            g, samples = part.split(":", 1)
            for s in samples.split(","):
                group_by[s.strip()] = g.strip()
        print(f"layer-3 grouping: {group_by}")
    else:
        print("layer-3 grouping: none (remodeling skipped)")

    baselines, clusters, remodels = [], [], []
    for tissue in tissues:
        merged, auto = load(args.indir, tissue)
        cm = cell_metrics(merged)
        print(f"[{tissue}] HK normalizer mode: {cm['hk_mode']}")
        base = layer1_baseline(merged, auto, cm, tissue)
        baselines.append(base)
        print(f"[{tissue}] layer 1: {len(base)} cell types")
        clus = layer2_clusters(merged, auto, cm, tissue, args.fe_cut,
                               args.n_boot, args.seed)
        clusters.append(clus)
        n_out = int((clus["label"] == "te_outlier").sum())
        n_col = int(clus["label"].str.startswith("collapsed").sum())
        print(f"[{tissue}] layer 2: {len(clus)} clusters, "
              f"{n_out} te_outlier, {n_col} collapsed")
        remod = None
        if group_by:
            remod = layer3_remodeling(merged, auto, cm, tissue, group_by)
            if remod is not None:
                remodels.append(remod)
                print(f"[{tissue}] layer 3: {len(remod)} cell types compared")
        fig_baseline(cm, auto, tissue, args.outdir, cm["hk_mode"])
        fig_family_comp(cm, merged, auto, tissue, args.outdir)
        fig_cluster_dev(clus, tissue, args.outdir, args.fe_cut)
        if remod is not None:
            fig_remodeling(remod, tissue, args.outdir, sorted(set(group_by.values())))

    pd.concat(baselines).to_csv(f"{args.outdir}/te_profile_baseline.csv", index=False)
    pd.concat(clusters).to_csv(f"{args.outdir}/te_profile_clusters.csv", index=False)
    if remodels:
        pd.concat(remodels).to_csv(f"{args.outdir}/te_profile_remodeling.csv", index=False)
    print("done ->", args.outdir)


if __name__ == "__main__":
    main()
