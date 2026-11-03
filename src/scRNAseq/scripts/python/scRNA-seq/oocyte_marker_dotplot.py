#!/usr/bin/env python3
"""Marker dotplot identity check for a target cell type (no UMAP).

Quantitative counterpart to the marker UMAP panels: a cell-type x gene
dotplot that decodes a label claim directly —

  rows    = cell types (target pinned on top, rest sorted by target-
            program expression), row label carries n cells;
  columns = marker groups: positive identity signature first, then
            lineage controls the target is expected NOT to express;
  colour  = mean expression scaled per gene across cell types
            (each gene column reads as a relative comparison);
  size    = fraction of expressing cells.

A companion TSV records per-cell-type mean/pct for every gene so every
dot is auditable without reading colour scales.

Usage:
  python oocyte_marker_dotplot.py \
      --auto-h5ad /path/to/ovaries_scTE_auto.h5ad \
      --target-cell-type Oocyte \
      --out-dir /path/to/ovaries_scTE_auto_reports/plot/marker

The default marker set validates the Oocyte (C14d) annotation of the
luanchao_mulatta ovaries dataset.
"""

import argparse
import csv
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

# (group header, role, genes) — role: "positive" = identity signature,
# "negative" = expected absent in the target cell type.
DEFAULT_MARKERS = [
    ("oocyte identity (positive)", "positive",
     ["NLRP5", "NLRP4", "NLRP2", "ZP3", "DDX4", "FIGLA", "NOBOX"]),
    ("proliferation (exp. absent)", "negative", ["MKI67", "TOP2A"]),
    ("granulosa (exp. absent)", "negative", ["AMH", "FST"]),
    ("endothelial (exp. absent)", "negative", ["PECAM1", "VWF"]),
    ("stromal (exp. absent)", "negative", ["COL1A1", "DCN"]),
]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--auto-h5ad", "-i", required=True,
                        help="Annotated scanpy auto h5ad (obs cell_type + raw slot)")
    parser.add_argument("--target-cell-type", "-t", required=True,
                        help="Cell type label to validate (exact obs['cell_type'] value)")
    parser.add_argument("--out-dir", "-o", required=True,
                        help="Output directory (PNG + TSV)")
    parser.add_argument("--type-col", default="cell_type",
                        help="obs column holding cell type labels (default: cell_type)")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Figure DPI (default: 300)")
    return parser.parse_args(argv)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    adata = sc.read_h5ad(args.auto_h5ad)
    src = adata.raw.to_adata() if adata.raw is not None else adata
    labels = adata.obs[args.type_col].astype(str)

    present = set(src.var_names)
    groups = []
    for header, role, genes in DEFAULT_MARKERS:
        keep = [g for g in genes if g in present]
        missing = [g for g in genes if g not in present]
        if missing:
            print(f"WARNING: genes not in h5ad, skipped: {missing}")
        if keep:
            groups.append((header, role, keep))
    genes = [g for _, _, gs in groups for g in gs]
    if not genes:
        raise SystemExit("no marker genes present in the h5ad")

    if args.target_cell_type not in set(labels):
        raise SystemExit(
            f"target cell type '{args.target_cell_type}' not in obs['{args.type_col}']")

    # ── Per-cell-type mean / pct in one matrix pass ──
    ct_counts = labels.value_counts()
    cts = sorted(ct_counts.index.tolist())
    mat = src[:, genes].X
    means = pd.DataFrame(0.0, index=cts, columns=genes)
    pcts = pd.DataFrame(0.0, index=cts, columns=genes)
    for ct in cts:
        mask = (labels == ct).values
        sub = mat[mask]
        sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        means.loc[ct] = sub.mean(axis=0)
        pcts.loc[ct] = (sub > 0).mean(axis=0)

    # Row order: target on top, then descending positive-program mean.
    pos_genes = [g for _, role, gs in groups if role == "positive"
                 for g in gs]
    pos_score = means[pos_genes].mean(axis=1).sort_values(ascending=False)
    row_order = ([args.target_cell_type]
                 + [ct for ct in pos_score.index
                    if ct != args.target_cell_type])
    means, pcts = means.loc[row_order], pcts.loc[row_order]

    # Colour: per-gene 0-1 scaling across cell types.
    span = means.max(axis=0) - means.min(axis=0)
    span[span == 0] = 1.0
    scaled = (means - means.min(axis=0)) / span

    # ── Figure ──
    n_ct, n_g = len(row_order), len(genes)
    fig, ax = plt.subplots(figsize=(0.62 * n_g + 4.2, 0.42 * n_ct + 2.2))

    for yi, ct in enumerate(row_order):
        if ct == args.target_cell_type:
            ax.axhspan(yi - 0.5, yi + 0.5, color="#fdeecf", zorder=0)
        elif yi % 2:
            ax.axhspan(yi - 0.5, yi + 0.5, color="#f5f5f5", zorder=0)
        for xi, g in enumerate(genes):
            ax.scatter(xi, yi, s=16 + 240 * pcts.loc[ct, g],
                       c=[plt.get_cmap("Reds")(scaled.loc[ct, g])],
                       edgecolors="none", zorder=3)

    ax.set_xticks(range(n_g))
    ax.set_xticklabels(genes, rotation=38, ha="right", fontsize=10)
    ax.set_yticks(range(n_ct))
    ax.set_yticklabels([f"{ct} (n={ct_counts[ct]})" for ct in row_order],
                       fontsize=9)
    for tick, ct in zip(ax.get_yticklabels(), row_order):
        if ct == args.target_cell_type:
            tick.set_fontweight("bold")
            tick.set_color("#b3261e")
    ax.set_ylim(n_ct - 0.5, -0.5)
    ax.set_xlim(-0.6, n_g - 0.4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    # Group separators + section headers. The positive block is wide so
    # it takes a centered header above the axis; the negative-control
    # groups are 2 genes wide, where centered/rotated text collides —
    # put their headers BELOW the gene labels instead (top label
    # region stays reserved for the title).
    trans = ax.get_xaxis_transform()
    offset = 0
    bounds = []
    for header, role, gs in groups:
        bounds.append((offset, offset + len(gs) - 1, header))
        offset += len(gs)
        if offset < n_g:
            ax.axvline(offset - 0.5, color="#b0b0b0", lw=1, zorder=1)
    for x0, x1, header in bounds:
        if x1 - x0 >= 4:            # wide group: centered above the axis
            ax.text((x0 + x1) / 2.0, 1.03, header, transform=trans,
                    ha="center", va="bottom", fontsize=8.5, color="#555555")
        else:                         # narrow group: horizontal, below x labels
            # ~22-char header over a 2-gene column: wrap to two lines so
            # neighbours cannot touch (ha=center keeps each centered).
            first = header.split(" (")[0]
            ax.text((x0 + x1) / 2.0, -0.16, f"{first}\n(exp. absent)",
                    transform=trans, ha="center", va="top",
                    fontsize=8, color="#555555", linespacing=1.25)

    ax.set_title(f"{args.target_cell_type} — marker identity check (dotplot)",
                 fontsize=13, fontweight="bold", pad=30)

    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap("Reds"),
                               norm=plt.Normalize(0.0, 1.0))
    cbar = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.02)
    cbar.set_label("mean expression (scaled per gene)", fontsize=8.5)
    cbar.ax.tick_params(labelsize=8)
    fig.text(0.995, 0.012,
             "dot size = % cells expressing; exact values in companion TSV",
             ha="right", fontsize=7.5, color="#777777")

    safe = args.target_cell_type.replace(" ", "_").replace("/", "-")
    out_png = os.path.join(args.out_dir, f"{safe}_marker_dotplot.png")
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"saved: {out_png}")

    # ── Companion TSV (audit table) ──
    out_tsv = os.path.join(args.out_dir, f"{safe}_marker_dotplot.tsv")
    role_of = {g: role for _, role, gs in groups for g in gs}
    group_of = {g: header for header, _, gs in groups for g in gs}
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["cell_type", "n_cells", "gene", "group", "role",
                         "mean_expr", "pct_expr"])
        for ct in row_order:
            for g in genes:
                writer.writerow([ct, int(ct_counts[ct]), g, group_of[g],
                                 role_of[g],
                                 round(float(means.loc[ct, g]), 4),
                                 round(float(pcts.loc[ct, g]), 4)])
    print(f"saved: {out_tsv}")

    # Console verdict summary for the target row.
    tgt = args.target_cell_type
    print(f"\n{tgt} (n={ct_counts[tgt]}):")
    for g in pos_genes:
        print(f"  + {g:<8} mean={means.loc[tgt, g]:.3f}  pct={pcts.loc[tgt, g]:.2%}")
    for g in [g for _, role, gs in groups if role == "negative"
              for g in gs]:
        print(f"  - {g:<8} mean={means.loc[tgt, g]:.3f}  pct={pcts.loc[tgt, g]:.2%}")


if __name__ == "__main__":
    main()
