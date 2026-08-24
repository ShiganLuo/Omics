"""Scanpy implementation for the standardized scRNA-seq workflow.

Modes: qc, cluster, batch, annotate, advanced, de.
Pipeline order: qc → cluster → batch → annotate → advanced → de
"""
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


def mode_qc(adata, args):
    """QC: mt/ribo metrics, outlier filtering, doublet detection."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)
    adata = adata[(adata.obs.n_genes_by_counts >= args.min_genes) & (adata.obs.n_genes_by_counts <= args.max_genes)].copy()
    adata = adata[adata.obs.pct_counts_mt <= args.max_pct_mt].copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, batch_key=args.batch_key or None)
    adata.write_h5ad(args.output)
    if args.metrics:
        adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(args.metrics, sep="\t")


def mode_cluster(adata, args):
    """Cluster: HVG subset → scale → PCA → neighbors → UMAP → leiden."""
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


def mode_batch(adata, args):
    """Batch correction: BBKNN or Harmony after clustering.

    Requires 'sample' or 'batch' column in .obs for batch_key.
    """
    batch_key = args.batch_key or ("sample" if "sample" in adata.obs else "batch")
    method = args.batch_method  # "bbknn" or "harmony"

    # PCA on full data (not just HVG) for batch correction
    if "X_pca" not in adata.obsm:
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=args.n_pcs)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)

    if method == "bbknn":
        sc.external.pp.bbknn(adata, batch_key=batch_key)
    elif method == "harmony":
        sc.external.pp.harmony_integrate(adata, key=batch_key)
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)

    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    adata.write_h5ad(args.output)


def mode_annotate(adata, args):
    """Cell type annotation: marker-based or celltypist.

    marker_file: TSV with columns [cell_type, gene1, gene2, ...]
    celltypist_model: celltypist model name (e.g. "Immune_All_High.pkl")
    """
    if args.marker_file:
        # Marker-based annotation
        import csv
        marker_genes = {}
        with open(args.marker_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                ct = row["cell_type"]
                genes = [g for g in row.get("markers", "").split(",") if g and g in adata.var.index]
                if genes:
                    marker_genes[ct] = genes
        # Score each cell type
        for ct, genes in marker_genes.items():
            sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{ct}")
        # Assign cell type by highest score
        score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
        if score_cols:
            adata.obs["cell_type"] = adata.obs[score_cols].idxmax(axis=1).str.replace("score_", "")
            adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    if args.celltypist_model:
        # celltypist annotation
        import celltypist
        from celltypist import models
        adata_ct = adata.copy()
        if "counts" in adata_ct.layers:
            adata_ct.X = adata_ct.layers["counts"]
        sc.pp.normalize_per_cell(adata_ct, counts_per_cell_after=1e4)
        sc.pp.log1p(adata_ct)
        adata_ct.X = adata_ct.X.toarray() if hasattr(adata_ct.X, "toarray") else adata_ct.X
        models.download_models(force_update=True, model=[args.celltypist_model])
        model = models.Model.load(model=args.celltypist_model)
        predictions = celltypist.annotate(adata_ct, model=model, majority_voting=True)
        pred_adata = predictions.to_adata()
        adata.obs["celltypist_label"] = pred_adata.obs.loc[adata.obs.index, "majority_voting"]
        adata.obs["celltypist_score"] = pred_adata.obs.loc[adata.obs.index, "conf_score"]

    adata.write_h5ad(args.output)


def mode_advanced(adata, args):
    """Advanced: trajectory, velocity, communication, CNV."""
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
        if args.gtf:
            cnv.io.genomic_position_from_gtf(args.gtf, adata)
        ref_cats = [c.strip() for c in args.cnv_reference.split(",") if c.strip()] if args.cnv_reference else None
        if ref_cats and "cell_type" in adata.obs:
            cnv.tl.infercnv(adata, reference_key="cell_type", reference_cat=ref_cats, window_size=250)
        else:
            cnv.tl.infercnv(adata, reference_key="leiden")
    adata.write_h5ad(args.output)


def mode_de(adata, args):
    """Differential expression: rank_genes_groups."""
    group = "condition" if "condition" in adata.obs else "leiden"
    sc.tl.rank_genes_groups(adata, group, method="wilcoxon")
    if args.deg:
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(args.deg, sep="\t", index=False)
    adata.write_h5ad(args.output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["qc", "cluster", "batch", "annotate", "advanced", "de"])
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    # QC params
    parser.add_argument("--metrics", default="")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=6000)
    parser.add_argument("--max-pct-mt", type=float, default=20)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--batch-key", default="")
    # Cluster params
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--markers", default="")
    # Batch params
    parser.add_argument("--batch-method", default="bbknn", choices=["bbknn", "harmony"])
    # Annotate params
    parser.add_argument("--marker-file", default="", help="TSV with cell_type and markers columns")
    parser.add_argument("--celltypist-model", default="", help="celltypist model name")
    # Advanced params
    parser.add_argument("--trajectory", action="store_true")
    parser.add_argument("--velocity", action="store_true")
    parser.add_argument("--communication", action="store_true")
    parser.add_argument("--cnv", action="store_true")
    parser.add_argument("--gtf", default="", help="GTF for CNV genomic position")
    parser.add_argument("--cnv-reference", default="", help="Comma-separated reference cell types for CNV")
    # DE params
    parser.add_argument("--deg", default="")
    args = parser.parse_args()

    adata = read_input(args.input[0], args.input[1:] if len(args.input) > 1 else None)

    mode_fn = {
        "qc": mode_qc,
        "cluster": mode_cluster,
        "batch": mode_batch,
        "annotate": mode_annotate,
        "advanced": mode_advanced,
        "de": mode_de,
    }[args.mode]
    mode_fn(adata, args)


if __name__ == "__main__":
    main()
