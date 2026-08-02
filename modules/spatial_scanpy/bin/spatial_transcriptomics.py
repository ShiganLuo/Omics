
"""Visium/custom h5ad spatial transcriptomics analysis."""
import argparse
from pathlib import Path
import scanpy as sc

def load_data(args):
    if args.input_h5ad:
        adata = sc.read_h5ad(args.input_h5ad)
    elif args.visium_h5:
        adata = sc.read_10x_h5(args.visium_h5)
        adata.var_names_make_unique()
        if args.spatial_dir:
            adata = sc.read_visium(args.spatial_dir, count_file=args.visium_h5, load_images=True)
            adata.var_names_make_unique()
    else:
        raise ValueError("Provide input_h5ad or visium_h5")
    if "spatial" not in adata.obsm:
        raise ValueError("Spatial coordinates are required in adata.obsm['spatial']")
    return adata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True)
    parser.add_argument("--input-h5ad", default="")
    parser.add_argument("--visium-h5", default="")
    parser.add_argument("--spatial-dir", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--metrics", default="")
    parser.add_argument("--genes", default="")
    parser.add_argument("--min-genes", type=int, default=100)
    parser.add_argument("--max-pct-mt", type=float, default=25)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=0.8)
    args = parser.parse_args()
    adata = load_data(args)
    if args.mode == "qc":
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        adata = adata[adata.obs.n_genes_by_counts >= args.min_genes].copy()
        adata = adata[adata.obs.pct_counts_mt <= args.max_pct_mt].copy()
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes)
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
        sc.tl.leiden(adata, resolution=args.resolution, key_added="spatial_leiden")
        sc.tl.rank_genes_groups(adata, "spatial_leiden", method="wilcoxon")
        if args.genes:
            sc.get.rank_genes_groups_df(adata, group=None).to_csv(args.genes, sep="\t", index=False)
        adata.write_h5ad(args.output)
    elif args.mode == "advanced":
        import squidpy as sq
        sq.gr.spatial_neighbors(adata)
        sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names[: min(adata.n_vars, 2000)])
        adata.write_h5ad(args.output)
    else:
        raise ValueError(f"Unsupported mode: {args.mode}")

if __name__ == "__main__":
    main()
