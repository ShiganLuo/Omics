
"""Scanpy implementation for the standardized scRNA-seq workflow."""
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import scanpy as sc

def read_input(path, sample_paths=None):
    if sample_paths:
        objects = [ad.read_h5ad(item) for item in sample_paths]
        return ad.concat(objects, join="outer", label="batch", keys=[Path(x).stem for x in sample_paths], fill_value=0)
    return ad.read_h5ad(path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True)
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--metrics", default="")
    parser.add_argument("--markers", default="")
    parser.add_argument("--deg", default="")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=6000)
    parser.add_argument("--max-pct-mt", type=float, default=20)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--batch-key", default="")
    parser.add_argument("--trajectory", action="store_true")
    parser.add_argument("--velocity", action="store_true")
    parser.add_argument("--communication", action="store_true")
    parser.add_argument("--cnv", action="store_true")
    args = parser.parse_args()
    adata = read_input(args.input[0], args.input[1:] if len(args.input) > 1 else None)
    if args.mode == "qc":
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        adata = adata[(adata.obs.n_genes_by_counts >= args.min_genes) & (adata.obs.n_genes_by_counts <= args.max_genes)].copy()
        adata = adata[adata.obs.pct_counts_mt <= args.max_pct_mt].copy()
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, batch_key=args.batch_key or None)
        adata.write_h5ad(args.output)
        if args.metrics:
            adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(args.metrics, sep="\t")
    elif args.mode == "cluster":
        if "highly_variable" in adata.var:
            adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=min(args.n_pcs, max(2, adata.n_obs - 1)))
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=min(args.n_pcs, adata.obsm["X_pca"].shape[1]))
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden")
        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
        if args.markers:
            sc.get.rank_genes_groups_df(adata, group=None).to_csv(args.markers, sep="\t", index=False)
        adata.write_h5ad(args.output)
    elif args.mode == "advanced":
        if args.trajectory:
            sc.tl.diffmap(adata)
            sc.tl.dpt(adata)
        if args.velocity and {"spliced", "unspliced"}.issubset(adata.layers):
            import scvelo as scv
            scv.pp.moments(adata, n_pcs=min(args.n_pcs, 30), n_neighbors=args.n_neighbors)
            scv.tl.velocity(adata)
            scv.tl.velocity_graph(adata)
        if args.communication:
            import liana as li
            li.mt.rank_aggregate(adata, groupby="leiden", use_raw=False, verbose=False)
        if args.cnv:
            import infercnvpy as cnv
            cnv.tl.infercnv(adata, reference_key="leiden")
        adata.write_h5ad(args.output)
    elif args.mode == "de":
        group = "condition" if "condition" in adata.obs else "leiden"
        sc.tl.rank_genes_groups(adata, group, method="wilcoxon")
        if args.deg:
            sc.get.rank_genes_groups_df(adata, group=None).to_csv(args.deg, sep="\t", index=False)
        adata.write_h5ad(args.output)
    else:
        raise ValueError(f"Unsupported mode: {args.mode}")

if __name__ == "__main__":
    main()
