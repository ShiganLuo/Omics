#!/usr/bin/env python3
"""Reproduce q6_dotplot.png (Q6 TE-dominant cluster marker dotplot) with an
EXPANDED senescence/aging marker panel.

Standalone reproduction of the dotplot originally produced by
te_dominant_q6_figures.py (which only had 3 senescence markers).  Per tissue
the TE-dominant cluster (highest mean per-cell TE fraction) is compared with
its same-lineage siblings and the rest of the cells.

Groups per tissue (x axis):  TE (dominant cluster) / siblings / rest
Markers (y axis), categories separated by grey lines:
  SASP / inflammation      IL6 CXCL8 CXCL1 CCL2 IL1A IL1B
  senescence / CDKi        CDKN1A CDKN2B CDKN1B RB1 TP53
                           GADD45A GADD45B SERPINE1
                           (CDKN2A absent from the mulatta annotation)
  DNA damage / aging       LMNB1 GLB1 HMGB1 IGFBP3 IGFBP7
                           FOXO1 FOXO3 SIRT1 SIRT6
                           (H2AFX absent from the mulatta annotation)
  proliferation            MKI67 TOP2A CDK1
  stress / apoptosis       BAX CASP3
  lineage identity         granulosa lineage (ovaries dominant cluster):
                           AMH CYP19A1 FSHR INHA GJA1
                           glandular epithelial lineage (Uterus dominant
                           cluster): EPCAM KRT19 MUC1 PAEP FOXA2
Genes absent from the annotation are skipped and listed on stdout.

Dot size = % expressing cells, color = mean log1p expression (auto.raw.X —
the full 19k-gene matrix; auto.X is only the 3000 HVGs and silently lacks
most senescence markers, which is why the original dotplot looked sparse).

Usage (scanpy container):
  apptainer exec /data/16T/Database/env/scanpy/scanpy.sif \
    python te_dominant_q6_figures/reproduce_q6_dotplot.py \
      --root /data/16T/luosg/output/luanchao_mulatta/scRNAseq \
      --out-dir /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/figures_te
"""

import argparse
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
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
    merged = sc.read_h5ad(f"{root}/common/5_combine_h5ad/{tissue}/{tissue}_scTE_merged.h5ad")
    auto = sc.read_h5ad(f"{root}/common/5_combine_h5ad/{tissue}/{tissue}_scTE_auto.h5ad")
    common = merged.obs_names.intersection(auto.obs_names)
    merged, auto = merged[common].copy(), auto[common].copy()
    return merged, auto


def dominant_cluster(merged, auto):
    """Return (dom_leiden, dom_cell_type, leiden, cell_type) with v2 gates.

    v2: eligibility by ABSOLUTE TE counts before ranking by TE% —
    a cluster whose median TE counts is < 0.5× same-lineage siblings
    median (or < 0.5× tissue median) is a gene-program-collapsed cell
    whose high TE% is a denominator artifact, not TE de-repression
    (caught Uterus Leiden 2: TE% 38.5 but TE abs 918 vs siblings 4129).
    """
    X = merged.layers["counts"] if "counts" in merged.layers else merged.X
    is_te = merged.var["gene_type"].astype(str).values == "TE"
    te_sum = np.asarray(X[:, is_te].sum(1)).ravel()
    te_abs = te_sum.astype(float)
    te_frac = te_sum / np.maximum(np.asarray(X.sum(1)).ravel(), 1) * 100
    leiden = auto.obs["leiden"].astype(str).values
    cell_type = auto.obs["cell_type"].astype(str).values
    order = sorted(set(leiden))
    mean_te = {c: te_frac[leiden == c].mean() for c in order}
    med_abs = {c: float(np.median(te_abs[leiden == c])) for c in order}
    n_cells = {c: int((leiden == c).sum()) for c in order}
    dom_ct_of = {c: pd_series_mode(cell_type[leiden == c]) for c in order}
    global_med = float(np.median(te_abs))
    dom = None
    for c in sorted(order, key=lambda c: -mean_te[c]):
        ct = dom_ct_of[c]
        sib = te_abs[(cell_type == ct) & (leiden != c)]
        sib_med = float(np.median(sib)) if len(sib) >= 30 else float("inf")
        if (med_abs[c] >= 0.5 * sib_med and med_abs[c] >= 0.5 * global_med
                and n_cells[c] >= 100):
            dom = c
            break
        if mean_te[c] > 35:
            print(f"[gate] {tissue_label}: Leiden {c} rejected: TE% "
                  f"{mean_te[c]:.1f} but TE abs median {med_abs[c]:.0f} vs "
                  f"siblings {sib_med:.0f} (denominator-collapse guard)")
    if dom is None:  # no eligible cluster: fall back to pure TE% max
        dom = max(order, key=lambda c: mean_te[c])
    ct = dom_ct_of[dom]
    print(f"{tissue_label}: dominant = Leiden {dom} ({ct}), "
          f"mean TE% = {mean_te[dom]:.2f}, n = {n_cells[dom]}")
    return dom, ct, leiden, cell_type


def pd_series_mode(arr):
    vals, counts = np.unique(arr, return_counts=True)
    return str(vals[counts.argmax()])


tissue_label = ""


def main():
    global tissue_label
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="scRNAseq analysis root (contains common/5_combine_h5ad)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-name", default="q6_dotplot.png")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    marker_list = [(g, n) for g, names in MARKERS.items() for n in names]
    marker_list += [(cat, n) for cat, genes in LINEAGE.items() for n in genes]
    cat_bounds = []          # y indices where a new category starts
    seen = set()
    for i, (g, _) in enumerate(marker_list):
        if g not in seen:
            cat_bounds.append(i)
            seen.add(g)
    y_labels = [f"{n} ({g})" for g, n in marker_list]

    fig, ax = plt.subplots(figsize=(13, 11))
    # 无背景、无边框:白底 + 四边 spines 全隐(刻度标签保留)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)
    xticks, xlabels = [], []
    col = 0
    missing = set()
    tissue_info = {}
    for tissue in ("ovaries", "Uterus"):
        tissue_label = tissue
        merged, auto = load(args.root, tissue)
        dom, dom_ct, leiden, cell_type = dominant_cluster(merged, auto)
        tissue_info[tissue] = (dom, dom_ct)
        dom_mask = leiden == dom
        ct_mask = cell_type == dom_ct
        sib_mask = ct_mask & ~dom_mask
        rest_mask = ~ct_mask

        ref = auto.raw if auto.raw is not None else auto
        name2idx = {}
        for j, g in enumerate(ref.var_names):
            name2idx.setdefault(g, j)
            if "." in g:
                name2idx.setdefault(g.split(".")[0], j)

        Xn = ref.X
        for grp, mask in (("TE", dom_mask), ("siblings", sib_mask),
                          ("rest", rest_mask)):
            xticks.append(col)
            xlabels.append(f"{tissue}\nLeiden {dom}\n{grp}")
            for (cat, gene), y in zip(marker_list, range(len(marker_list))):
                j = name2idx.get(gene)
                if j is None:
                    missing.add(gene)
                    continue
                v = Xn[mask][:, j]
                v = np.asarray(v.todense() if hasattr(v, "todense") else v).ravel()
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
    for b in cat_bounds[1:]:
        ax.axhline(b - 0.5, color="#cccccc", lw=0.8)
    ov_dom, ov_ct = tissue_info["ovaries"]
    ut_dom, ut_ct = tissue_info["Uterus"]
    fig.text(0.5, 0.012,
             f"TE = TE-dominant cluster (highest mean TE fraction: "
             f"ovaries Leiden {ov_dom}, {ov_ct} / Uterus Leiden {ut_dom}, {ut_ct})\n"
             f"siblings = other Leiden clusters of the same cell type "
             f"({ov_ct} in ovaries / {ut_ct} in Uterus)\n"
             f"rest = cells of all other cell types "
             f"(not {ov_ct.lower()} in ovaries / not {ut_ct.lower()} in Uterus)",
             ha="center", va="bottom", fontsize=9.5)
    fig.tight_layout(rect=(0, 0.085, 0.80, 1))
    # right-side legends, dot-size column centered on the colorbar axis
    xc = 0.87
    sm = mcm.ScalarMappable(cmap="RdYlBu_r", norm=Normalize(0, 1.8))
    cb_ax = fig.add_axes((xc - 0.009, 0.30, 0.018, 0.28))
    cb = fig.colorbar(sm, cax=cb_ax)
    cb.set_label("mean log1p", fontsize=9)
    s_ax = fig.add_axes((xc - 0.09, 0.62, 0.18, 0.24))
    s_ax.axis("off")
    s_ax.set_xlim(0, 1)
    s_ax.set_ylim(0, 1)
    # title centered over the whole block (dots + % labels), not just the dots
    s_ax.text(0.62, 1.04, "% expressing", ha="center", va="bottom", fontsize=9)
    for y, p in zip((0.82, 0.58, 0.34, 0.10), (100, 75, 50, 25)):
        s_ax.scatter([0.5], [y], s=p * 2.5 + 5, c="#999999",
                     edgecolors="k", linewidths=0.3, clip_on=False)
        s_ax.text(0.66, y, f"{p}%", ha="left", va="center", fontsize=9)
    out = os.path.join(args.out_dir, args.out_name)
    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    if missing:
        print("genes not in annotation (skipped):", ", ".join(sorted(missing)))
    print("done ->", out)


if __name__ == "__main__":
    main()
