#!/usr/bin/env python3
"""Prepare pseudo-bulk count matrices for DEG analysis.

Adapted from DEG.py (Omics/scRNAseq/scripts/python/utils/DEG.py).
For each cell type, aggregate single-cell counts into pseudo-bulk samples
(summed UMI per sample), one row per sample, one column per gene.

Output:
  {out_dir}/pseudobulk_counts/{label}_{cell_type}_counts.csv

Usage:
  python pseudobulk_prepare.py \
    --h5ad /path/to/scTE_auto.h5ad \
    --out-dir /path/to/results \
    --label ovaries \
    --min-cells 30
"""
import argparse, os, sys, warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import scanpy as sc


def aggregate_pseudobulk(
    adata,
    cell_type,
    sample_col="sample_id",
    type_col="cell_type",
    min_cells=30,
):
    """Aggregate single cells per sample -> pseudo-bulk count matrix.

    Returns
    -------
    df : pd.DataFrame  rows=samples, cols=genes  (gene columns start after meta cols)
    samples_kept : list[str]
    """
    sub = adata[adata.obs[type_col] == cell_type].copy()
    size_by_sample = sub.obs.groupby(sample_col).size()
    samples_kept = size_by_sample[size_by_sample > min_cells].index.tolist()
    if len(samples_kept) < 2:
        return None, [], []

    raw_sub = sub.raw.to_adata()
    rows = []
    for sample in samples_kept:
        mask = sub.obs[sample_col] == sample
        counts = np.array(raw_sub[mask].X.sum(axis=0)).flatten()
        rows.append({
            "sample": sample,
            "cell_type": cell_type,
            "n_cells": int(mask.sum()),
            **{gene: float(c) for gene, c in zip(raw_sub.var_names, counts)},
        })
    # Ensure 'sample' column is preserved (don't let it be eaten as index)
    df = pd.DataFrame(rows)
    if "sample" in df.columns:
        df = df.set_index("sample")
        df.index.name = "sample"
    return df, samples_kept, list(sub.obs[sample_col].unique())


def run(h5ad_path, out_dir, label, min_cells=30, cell_types=None):
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{label}] Loading {h5ad_path} ...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  shape: {adata.shape}, samples: {adata.obs['sample_id'].nunique()}")

    targets = cell_types or sorted(adata.obs["cell_type"].unique().tolist())
    written = []
    for ct in targets:
        df, kept, all_samples = aggregate_pseudobulk(
            adata, ct, min_cells=min_cells,
        )
        if df is None:
            print(f"  [{ct}] skipped (samples with >{min_cells} cells: {len(all_samples)})")
            continue
        ct_safe = ct.replace(" ", "_").replace("/", "-")
        out_path = os.path.join(out_dir, f"{label}_{ct_safe}_counts.csv")
        df.to_csv(out_path)
        written.append({
            "label": label,
            "cell_type": ct,
            "n_samples": len(kept),
            "n_genes": df.shape[1] - 3,
            "path": out_path,
        })
        print(f"  [{ct}] n_samples={len(kept)} -> {out_path}")
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="Input h5ad file")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--label", required=True, help="Sample-set label (e.g. ovaries, uterus)")
    parser.add_argument("--min-cells", type=int, default=30,
                        help="Minimum cells per sample (default 30)")
    parser.add_argument("--cell-types", nargs="+", default=None,
                        help="Subset of cell types (default all)")
    args = parser.parse_args()

    written = run(args.h5ad, args.out_dir, args.label, args.min_cells, args.cell_types)
    summary_path = os.path.join(args.out_dir, f"{args.label}_pseudobulk_summary.csv")
    pd.DataFrame(written).to_csv(summary_path, index=False)
    print(f"\nSummary saved: {summary_path}")
    print(f"Total files: {len(written)}")


if __name__ == "__main__":
    main()