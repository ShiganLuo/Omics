"""Scanpy implementation for the standardized scRNA-seq workflow.

Modes: qc, merge, cluster, batch, annotate, advanced, de.
Pipeline order: qc(each sample) -> merge(by tissue) -> cluster -> batch -> annotate -> advanced -> de
"""
import argparse
import csv
import json
import logging
import os
import sys
from pathlib import Path
from typing import List, Optional

import anndata as ad
import numpy as np
import scanpy as sc

ad.settings.allow_write_nullable_strings = True

# Lazy import: plotter loaded only when --plot-dir is set
_plotter_cls = None


def _get_plotter():
    """Lazy-load ScanpyPlotter to avoid matplotlib import at startup."""
    global _plotter_cls
    if _plotter_cls is None:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from plot import ScanpyPlotter
        _plotter_cls = ScanpyPlotter
    return _plotter_cls


def setup_logging(level=logging.INFO):
    """Setup logging configuration."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_input(path: str, sample_paths: Optional[List[str]] = None) -> ad.AnnData:
    """Read a single h5ad or concatenate multiple files."""
    if sample_paths:
        objects = [ad.read_h5ad(item) for item in sample_paths]
        return ad.concat(
            objects, join="outer", label="batch",
            keys=[Path(x).stem for x in sample_paths], fill_value=0,
        )
    return ad.read_h5ad(path)


def _make_plotter(plot_dir: str):
    """Instantiate plotter if plot_dir is provided, else return None."""
    if not plot_dir:
        return None
    cls = _get_plotter()
    return cls(plot_dir)


# ---------------------------------------------------------------------------
# QC: per-sample filtering
# ---------------------------------------------------------------------------
def mode_qc(adata: ad.AnnData, args) -> None:
    """QC: mt/ribo metrics, outlier filtering, doublet detection."""
    adata.var["mt"] = np.array(
        adata.var_names.str.upper().str.startswith("MT-")
        | adata.var_names.str.startswith("mt")
    )
    adata.var["ribo"] = np.array(adata.var_names.str.startswith(("RPS", "RPL")))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)

    adata_before = adata.copy()

    n_before = adata.n_obs
    adata = adata[
        (adata.obs.n_genes_by_counts >= args.min_genes)
        & (adata.obs.n_genes_by_counts <= args.max_genes)
    ].copy()
    adata = adata[adata.obs.pct_counts_mt <= args.max_pct_mt].copy()

    # Doublet detection via scrublet
    if args.scrublet:
        try:
            sc.external.pp.scrublet(adata, expected_doublet_rate=args.doublet_rate)
            n_doublet = int(adata.obs["predicted_doublet"].sum())
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            logging.info("Scrublet: removed %d doublets", n_doublet)
        except Exception as exc:
            logging.warning("Scrublet failed (%s), skipping doublet detection", exc)

    n_after = adata.n_obs
    logging.info("QC: %d -> %d cells (%d removed)", n_before, n_after, n_before - n_after)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=args.n_top_genes,
        batch_key=args.batch_key or None,
    )

    adata.uns["qc_stats"] = {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_before - n_after,
        "pct_removed": round((n_before - n_after) / n_before * 100, 2) if n_before > 0 else 0,
    }

    plotter = _make_plotter(getattr(args, "plot_dir", ""))
    if plotter:
        plotter.plot_qc(adata, adata_before)
        del adata_before

    adata.write_h5ad(args.output)
    if args.metrics:
        adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(
            args.metrics, sep="\t"
        )


# ---------------------------------------------------------------------------
# Merge: merge multiple QC'd h5ad by tissue
# ---------------------------------------------------------------------------
def mode_merge(adata: ad.AnnData, args) -> None:
    """Merge multiple QC'd h5ad files into one, preserving sample info."""
    objects = []
    for p in args.input:
        obj = ad.read_h5ad(p)
        sample_id = (
            Path(p).stem
            .replace("_scTE", "")
            .replace("_cellranger", "")
            .replace("_qc", "")
        )
        obj.obs["sample_id"] = sample_id
        objects.append(obj)

    merged = ad.concat(
        objects, join="outer", label="sample_id",
        keys=[o.obs["sample_id"].iloc[0] for o in objects], fill_value=0,
    )

    if all("counts" in o.layers for o in objects):
        merged.layers["counts"] = merged.X.copy()

    for col in ("sample_id", "batch"):
        if col in merged.obs.columns:
            merged.obs[col] = merged.obs[col].astype("category")

    logging.info("Merged %d samples: %d cells x %d genes",
                 len(objects), merged.n_obs, merged.n_vars)

    plotter = _make_plotter(getattr(args, "plot_dir", ""))
    if plotter:
        plotter.plot_merge(merged)

    merged.write_h5ad(args.output)


# ---------------------------------------------------------------------------
# Cluster
# ---------------------------------------------------------------------------
def mode_cluster(adata: ad.AnnData, args) -> None:
    """Cluster: HVG subset -> scale -> PCA -> neighbors -> UMAP -> leiden."""
    if "highly_variable" in adata.var:
        adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(args.n_pcs, max(2, adata.n_obs - 1)))
    sc.pp.neighbors(
        adata, n_neighbors=args.n_neighbors,
        n_pcs=min(args.n_pcs, adata.obsm["X_pca"].shape[1]),
    )
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    if args.markers:
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(
            args.markers, sep="\t", index=False
        )

    plotter = _make_plotter(getattr(args, "plot_dir", ""))
    if plotter:
        plotter.plot_cluster(adata)

    adata.write_h5ad(args.output)


# ---------------------------------------------------------------------------
# Batch correction
# ---------------------------------------------------------------------------
def mode_batch(adata: ad.AnnData, args) -> None:
    """Batch correction: BBKNN or Harmony after clustering."""
    batch_key = args.batch_key or ("sample" if "sample" in adata.obs else "batch")
    method = args.batch_method

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
    sc.tl.leiden(
        adata, resolution=args.resolution, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )

    plotter = _make_plotter(getattr(args, "plot_dir", ""))
    if plotter:
        plotter.plot_batch(adata, method)

    adata.write_h5ad(args.output)


# ---------------------------------------------------------------------------
# Cell type annotation
# ---------------------------------------------------------------------------
def mode_annotate(adata: ad.AnnData, args) -> None:
    """Cell type annotation: marker-based, celltypist, or LLM-assisted."""
    if args.marker_file:
        _annotate_markers(adata, args)
    if args.celltypist_model:
        _annotate_celltypist(adata, args)
    if args.llm_method:
        _annotate_llm(adata, args)

    plotter = _make_plotter(getattr(args, "plot_dir", ""))
    if plotter:
        plotter.plot_annotate(
            adata,
            marker_file=args.marker_file,
            annotate_group=args.annotate_group,
        )

    adata.write_h5ad(args.output)


def _annotate_markers(adata: ad.AnnData, args) -> None:
    """Marker-based annotation using sc.tl.score_genes."""
    marker_genes = {}
    with open(args.marker_file, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ct = row["cell_type"]
            genes = [g for g in row.get("markers", "").split(",") if g and g in adata.var.index]
            if genes:
                marker_genes[ct] = genes

    for ct, genes in marker_genes.items():
        sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{ct}")

    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if score_cols:
        adata.obs["cell_type"] = adata.obs[score_cols].idxmax(axis=1).str.replace("score_", "")
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")


def _annotate_celltypist(adata: ad.AnnData, args) -> None:
    """CellTypist annotation."""
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


def _annotate_llm(adata: ad.AnnData, args) -> None:
    """LLM-assisted cell type annotation."""
    group_key = args.annotate_group or "leiden"
    if group_key not in adata.obs.columns:
        logging.warning("Group key '%s' not found in obs, skipping LLM annotation", group_key)
        return

    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, group_key, method="wilcoxon")

    n_genes = args.llm_top_genes
    cluster_markers = {}
    for cluster in adata.obs[group_key].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        cluster_markers[str(cluster)] = df.head(n_genes)["names"].tolist()

    tissue = adata.uns.get("tissue", args.tissue or "unknown tissue")
    annotation = _llm_annotate_clusters(cluster_markers, tissue, args)

    label_map = {k: v.get("cell_type", "unknown") for k, v in annotation.items()}
    adata.obs["llm_label"] = adata.obs[group_key].map(label_map).astype("category")
    adata.uns["llm_annotation"] = annotation

    logging.info("LLM annotation complete: %d clusters annotated", len(annotation))
    for cluster, info in annotation.items():
        logging.info("  Cluster %s: %s (confidence: %s)",
                     cluster, info.get("cell_type", "?"), info.get("confidence", "?"))


def _llm_annotate_clusters(cluster_markers: dict, tissue: str, args) -> dict:
    """Call LLM to annotate clusters based on marker genes."""
    prompt = _build_annotation_prompt(cluster_markers, tissue)

    if args.llm_method == "openai":
        return _call_openai(prompt, args)
    elif args.llm_method == "ollama":
        return _call_ollama(prompt, args)
    elif args.llm_method == "file":
        output_file = args.output.replace(".h5ad", "_llm_prompt.txt")
        with open(output_file, "w", encoding="utf-8") as fh:
            fh.write(prompt)
        logging.info("LLM prompt written to %s", output_file)
        return {}
    else:
        logging.warning("Unknown LLM method: %s", args.llm_method)
        return {}


def _build_annotation_prompt(cluster_markers: dict, tissue: str) -> str:
    """Build a structured prompt for LLM-based cell type annotation."""
    prompt = f"""You are an expert single-cell RNA-seq analyst. I need you to annotate cell types for clusters from a {tissue} dataset.

For each cluster below, I provide the top differentially expressed genes (DEG).
Please output a JSON object where each key is the cluster ID, and each value has:
- "cell_type": the predicted cell type
- "confidence": "high", "medium", or "low"
- "reasoning": brief explanation of your reasoning
- "canonical_markers": known canonical markers you used for this identification

Important notes:
- Consider species-specific marker conventions
- If a cluster could be multiple cell types, assign the most likely one and note alternatives in reasoning
- For ambiguous clusters, set confidence to "low"
- Output ONLY valid JSON, no markdown formatting

Cluster markers:
"""
    for cluster, genes in cluster_markers.items():
        prompt += f"\nCluster {cluster}: {', '.join(genes)}"
    prompt += "\n\nRespond with a single JSON object:"
    return prompt


def _call_openai(prompt: str, args) -> dict:
    """Call OpenAI API for annotation."""
    try:
        import openai
        client = openai.OpenAI(api_key=args.llm_api_key, base_url=args.llm_base_url)
        response = client.chat.completions.create(
            model=args.llm_model or "gpt-4o",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
            response_format={"type": "json_object"},
        )
        return json.loads(response.choices[0].message.content)
    except Exception as exc:
        logging.error("OpenAI API call failed: %s", exc)
        return {}


def _call_ollama(prompt: str, args) -> dict:
    """Call local Ollama for annotation."""
    try:
        import requests
        url = args.llm_base_url or "http://localhost:11434"
        resp = requests.post(
            f"{url}/api/chat",
            json={
                "model": args.llm_model or "llama3.1",
                "messages": [{"role": "user", "content": prompt}],
                "stream": False,
                "format": "json",
            },
        )
        return json.loads(resp.json()["message"]["content"])
    except Exception as exc:
        logging.error("Ollama API call failed: %s", exc)
        return {}


# ---------------------------------------------------------------------------
# Advanced analysis
# ---------------------------------------------------------------------------
def mode_advanced(adata: ad.AnnData, args) -> None:
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
        ref_cats = (
            [c.strip() for c in args.cnv_reference.split(",") if c.strip()]
            if args.cnv_reference else None
        )
        if ref_cats and "cell_type" in adata.obs:
            cnv.tl.infercnv(adata, reference_key="cell_type",
                            reference_cat=ref_cats, window_size=250)
        else:
            cnv.tl.infercnv(adata, reference_key="leiden")

    plotter = _make_plotter(getattr(args, "plot_dir", ""))
    if plotter:
        plotter.plot_advanced(adata, trajectory=args.trajectory, cnv=args.cnv)

    adata.write_h5ad(args.output)


# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------
def mode_de(adata: ad.AnnData, args) -> None:
    """Differential expression: rank_genes_groups."""
    group = "condition" if "condition" in adata.obs else "leiden"
    sc.tl.rank_genes_groups(adata, group, method="wilcoxon")

    if args.deg:
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(
            args.deg, sep="\t", index=False
        )

    plotter = _make_plotter(getattr(args, "plot_dir", ""))
    if plotter:
        plotter.plot_de(adata, group)

    adata.write_h5ad(args.output)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Scanpy scRNA-seq pipeline")
    parser.add_argument("--mode", required=True,
                        choices=["qc", "merge", "cluster", "batch", "annotate", "advanced", "de"])
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--plot-dir", default="", help="Directory to save plots (optional)")

    # QC params
    parser.add_argument("--metrics", default="")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=6000)
    parser.add_argument("--max-pct-mt", type=float, default=20)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--batch-key", default="")
    parser.add_argument("--scrublet", action="store_true")
    parser.add_argument("--doublet-rate", type=float, default=0.06)

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
    parser.add_argument("--llm-method", default="", choices=["", "openai", "ollama", "file"])
    parser.add_argument("--llm-model", default="")
    parser.add_argument("--llm-api-key", default="")
    parser.add_argument("--llm-base-url", default="")
    parser.add_argument("--llm-top-genes", type=int, default=30)
    parser.add_argument("--annotate-group", default="")
    parser.add_argument("--tissue", default="")

    # Advanced params
    parser.add_argument("--trajectory", action="store_true")
    parser.add_argument("--velocity", action="store_true")
    parser.add_argument("--communication", action="store_true")
    parser.add_argument("--cnv", action="store_true")
    parser.add_argument("--gtf", default="")
    parser.add_argument("--cnv-reference", default="")

    # DE params
    parser.add_argument("--deg", default="")

    args = parser.parse_args()
    adata = read_input(args.input[0], args.input[1:] if len(args.input) > 1 else None)

    mode_fn = {
        "qc": mode_qc,
        "merge": mode_merge,
        "cluster": mode_cluster,
        "batch": mode_batch,
        "annotate": mode_annotate,
        "advanced": mode_advanced,
        "de": mode_de,
    }[args.mode]
    mode_fn(adata, args)


if __name__ == "__main__":
    main()
