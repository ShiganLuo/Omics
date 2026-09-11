"""Scanpy implementation for the standardized scRNA-seq workflow.

Modes: qc, merge, cluster, annotate, auto, advanced, de.
Pipeline order: qc(each sample) -> merge(by tissue) -> cluster -> annotate -> advanced -> de
Auto mode: cluster -> AI annotate -> QC -> filter -> re-cluster (iterative)
"""
import argparse
import csv
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple
from scipy.stats import median_abs_deviation
import anndata as ad
import numpy as np
import pandas as pd
import harmonypy as hm
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


def setup_logging(level: int = logging.INFO) -> None:
    """Configure the root logger with a timestamped format.

    Args:
        level: Logging level (default: INFO).
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_input(path: str, sample_paths: Optional[List[str]] = None) -> ad.AnnData:
    """Read a single h5ad file or concatenate multiple h5ad files.

    Args:
        path: Path to a single h5ad file. Ignored when *sample_paths* is given.
        sample_paths: Optional list of h5ad file paths to concatenate along
            the obs axis (outer join, fill_value=0).

    Returns:
        An AnnData object containing the loaded (or merged) data.
    """
    if sample_paths:
        objects = [ad.read_h5ad(item) for item in sample_paths]
        return ad.concat(
            objects, join="outer", label="batch",
            keys=[Path(x).stem for x in sample_paths], fill_value=0,
        )
    return ad.read_h5ad(path)


def _make_plotter(plot_dir: str):
    """Instantiate a ScanpyPlotter if *plot_dir* is non-empty, else return None.

    Args:
        plot_dir: Directory path for saving plots. Empty string disables plotting.

    Returns:
        A ScanpyPlotter instance or None.
    """
    if not plot_dir:
        return None
    cls = _get_plotter()
    return cls(plot_dir)


def is_outlier(adata: ad.AnnData, metric: str, nmads: int) -> np.ndarray:
    """Detect outlier cells using the median absolute deviation (MAD) method.

    A cell is flagged as an outlier if its *metric* value falls outside
    ``median ± nmads * MAD``.

    Args:
        adata: AnnData object whose ``obs`` contains the *metric* column.
        metric: Column name in ``adata.obs`` to evaluate.
        nmads: Number of MADs from the median to use as the threshold.

    Returns:
        Boolean array of shape ``(n_cells,)`` — True for outlier cells.
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def detect_n_pcs(
    variance_ratio: np.ndarray,
    min_pcs: int = 10,
    max_pcs: int = 100,
    window: int = 5,
    ratio: float = 0.15,
) -> Tuple[int, Dict]:
    """Detect the optimal number of principal components.

    Uses a sliding-window approach to find where the descending
    variance ratio curve starts to plateau.

    Algorithm:
        1. Compute per-PC change: ``delta = abs(diff(variance_ratio))``.
        2. Compute baseline: median of *delta* in the first *min_pcs* PCs.
        3. Slide a window of size *window* from *min_pcs* onward; at
           each position compute the mean *delta* inside the window.
        4. The first position where ``mean_delta < baseline * ratio``
           marks the start of the plateau → that PC index is *n_pcs*.

    Args:
        variance_ratio: 1-D array of per-PC variance ratios.
        min_pcs: Minimum PCs to return (default 10).
        max_pcs: Maximum PCs to return (default 100).
        window: Sliding window size (default 5).
        ratio: Threshold as a fraction of the baseline change rate
            (default 0.15).  Smaller = stricter (more PCs selected).

    Returns:
        Tuple of (recommended_n_pcs, diagnostics) where diagnostics is a
        dict with keys ``delta``, ``window_mean_x``, ``window_mean_y``,
        ``threshold``, ``elbow_pc`` for plotting.
    """
    n = len(variance_ratio)
    if n <= min_pcs:
        diag: Dict = {
            "delta": np.abs(np.diff(variance_ratio)) if n > 1 else np.array([]),
            "window_mean_x": np.array([]),
            "window_mean_y": np.array([]),
            "threshold": 0.0,
            "elbow_pc": n,
        }
        return n, diag

    delta = np.abs(np.diff(variance_ratio))  # length n-1

    # Baseline: median change in the first min_pcs PCs (active decline region)
    baseline_end = min(min_pcs, len(delta))
    baseline = np.median(delta[:baseline_end])
    if baseline == 0:
        diag = {
            "delta": delta,
            "window_mean_x": np.array([]),
            "window_mean_y": np.array([]),
            "threshold": 0.0,
            "elbow_pc": min(max_pcs, n),
        }
        return min(max_pcs, n), diag

    threshold = baseline * ratio

    # Compute sliding window mean for all valid positions
    search_start = max(0, min_pcs - 1)
    wm_list: List = []
    result_pc = min(max_pcs, n)  # fallback

    for i in range(search_start, len(delta) - window + 1):
        wm = float(np.mean(delta[i:i + window]))
        wm_list.append((i, wm))
        if wm < threshold and result_pc == min(max_pcs, n):
            result_pc = max(min_pcs, min(i + 1, max_pcs, n))

    # Build arrays for plotting
    if wm_list:
        wm_x = np.array([w[0] + 1 for w in wm_list])  # 1-indexed PC
        wm_y = np.array([w[1] for w in wm_list])
    else:
        wm_x = np.array([])
        wm_y = np.array([])

    diag = {
        "delta": delta,
        "window_mean_x": wm_x,
        "window_mean_y": wm_y,
        "threshold": threshold,
        "elbow_pc": result_pc,
    }
    return result_pc, diag


# ---------------------------------------------------------------------------
# QC: per-sample filtering
# ---------------------------------------------------------------------------
def mode_qc(
    adata: ad.AnnData,
    output: str,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_pct_mt: float = 20.0,
    use_mad: bool = False,
    scrublet: bool = False,
    doublet_rate: float = 0.06,
    plot_dir: str = "",
    metrics: str = "",
) -> None:
    """Quality-control pipeline for a single sample.

    Steps:
        1. Compute mitochondrial (MT-), ribosomal (RPS/RPL), and
           hemoglobin (HB) gene fraction metrics.
        2. (Optional) Flag outlier cells via MAD-based thresholds on total
           counts, gene counts, and top-20 gene percentage (5 MADs).
        3. Apply hard filters: min_genes, max_genes, max_pct_mt.
        4. Optionally run Scrublet doublet detection and remove predicted
           doublets.
        5. Store a raw counts layer and QC statistics in ``adata.uns``.
        6. Optionally generate QC plots and write a metrics TSV.

    Args:
        adata: Input AnnData object (raw counts expected in ``.X``).
        output: Path to write the filtered h5ad file.
        min_genes: Minimum genes per cell (hard filter, default 200).
        max_genes: Maximum genes per cell (hard filter, default 6000).
        max_pct_mt: Maximum mitochondrial percentage (hard filter, default 20).
        use_mad: If True, apply MAD-based outlier detection before hard
            filters (more permissive strategy). Default False (hard filters only).
        scrublet: If True, run Scrublet doublet detection.
        doublet_rate: Expected doublet rate passed to Scrublet (default 0.06).
        plot_dir: Directory for QC plots. Empty string disables plotting.
        metrics: Path to write a TSV of per-cell QC metrics. Empty string
            skips this step.
    """
    adata.var["mt"] = np.array(
        adata.var_names.str.upper().str.startswith("MT-")
        | adata.var_names.str.startswith("mt")
    )
    adata.var["ribo"] = np.array(adata.var_names.str.startswith(("RPS", "RPL")))
    adata.var["hb"] = np.array(adata.var_names.str.contains(r"^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

    adata_before = adata.copy()

    n_before = adata.n_obs

    # MAD-based outlier detection (optional, more permissive strategy)
    if use_mad:
        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
        )
        n_before_mad = adata.n_obs
        adata = adata[~adata.obs["outlier"]].copy()
        n_after_mad = adata.n_obs
        if n_before_mad != n_after_mad:
            logging.info("MAD outlier detection: %d -> %d cells", n_before_mad, n_after_mad)

    # Hard filters (always applied)
    n_before_hard = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mt)
    ].copy()
    n_after_hard = adata.n_obs
    if n_before_hard != n_after_hard:
        logging.info("Hard filter (genes %d-%d, mt %.1f%%): %d -> %d cells",
                     min_genes, max_genes, max_pct_mt, n_before_hard, n_after_hard)
    # Doublet detection via scrublet
    if scrublet:
        try:
            sc.external.pp.scrublet(adata, expected_doublet_rate=doublet_rate)
            n_doublet = int(adata.obs["predicted_doublet"].sum())
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            logging.info("Scrublet: removed %d doublets", n_doublet)
        except Exception as exc:
            logging.warning("Scrublet failed (%s), skipping doublet detection", exc)

    n_after = adata.n_obs
    logging.info("QC: %d -> %d cells (%d removed)", n_before, n_after, n_before - n_after)

    adata.uns["qc_stats"] = {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_before - n_after,
        "pct_removed": round((n_before - n_after) / n_before * 100, 2) if n_before > 0 else 0,
    }

    plotter = _make_plotter(plot_dir)
    if plotter:
        plotter.plot_qc(adata, adata_before,
                        counts_col="total_counts",
                        genes_col="n_genes_by_counts",
                        mt_col="pct_counts_mt")
        del adata_before

    adata.write_h5ad(output)
    if metrics:
        adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(
            metrics, sep="\t"
        )


# ---------------------------------------------------------------------------
# Merge: merge multiple QC'd h5ad by tissue
# ---------------------------------------------------------------------------
def mode_merge(
    adata: ad.AnnData,
    input_paths: List[str],
    output: str,
    axis: Literal["obs", 0, "var", 1] = "obs",
    plot_dir: str = "",
    te_bed: str = "",
    gene_tsv: str = "",
) -> None:
    """Merge multiple QC'd h5ad files into a single AnnData object.

    Each input file is tagged with a ``sample_id`` derived from its
    filename stem (with common suffixes like ``_scTE`` / ``_cellranger``
    stripped).  Concatenation uses an outer join so that genes present in
    only some samples are filled with zero.

    Args:
        adata: Unused placeholder (kept for uniform ``mode_fn(adata, …)``
            dispatch). The function reads files directly from *input_paths*.
        input_paths: List of h5ad file paths to merge.
        output: Path to write the merged h5ad file.
        axis: Concatenation axis — ``"obs"`` / ``0`` (default) merges cells;
            ``"var"`` / ``1`` merges genes.
        plot_dir: Directory for merge summary plots. Empty string disables
            plotting.
        te_bed: Path to TE BED file for gene_type annotation. Empty string
            skips annotation.
        gene_tsv: Path to gene annotation TSV for gene_type annotation.
            Empty string skips annotation.
    """
    objects = []
    for p in input_paths:
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
        objects, join="outer", label="sample_id", axis=axis,
        keys=[o.obs["sample_id"].iloc[0] for o in objects], merge="same"
    )
    merged.var_names_make_unique()
    merged.obs_names_make_unique()
    merged.layers["counts"] = merged.X.copy()

    for col in ("sample_id", "batch"):
        if col in merged.obs.columns:
            merged.obs[col] = merged.obs[col].astype("category")

    logging.info("Merged %d samples: %d cells x %d genes",
                 len(objects), merged.n_obs, merged.n_vars)

    if te_bed and gene_tsv:
        annotate_gene_type(merged, te_bed, gene_tsv)

    plotter = _make_plotter(plot_dir)
    if plotter:
        plotter.plot_merge(merged, sample_key="sample_id")

    merged.write_h5ad(output)


def _parse_te_bed(bed_path: str) -> set:
    """Parse a TE BED file and return a set of TE names.

    Expected format (no header, tab-separated):
        chrom\\tstart\\tend\\tTE_name

    Args:
        bed_path: Path to the TE BED file.

    Returns:
        Set of unique TE names from column 4.
    """
    te_names: set = set()
    with open(bed_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                te_names.add(parts[3])
    return te_names


def _parse_gene_tsv(tsv_path: str) -> Dict[str, str]:
    """Parse a gene annotation TSV and return gene_name -> gene_type mapping.

    Expected format (with header):
        gene_id\\tgene_name\\tgene_type

    When multiple rows share the same gene_name, the last non-empty
    ``gene_type`` wins.

    Args:
        tsv_path: Path to the gene annotation TSV/CSV.

    Returns:
        Dict mapping gene_name -> gene_type.
    """
    result: Dict[str, str] = {}
    with open(tsv_path, "r", encoding="utf-8") as fh:
        header = fh.readline()  # skip header
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gene_name = parts[1].strip()
            gene_type = parts[2].strip()
            if gene_name:
                result[gene_name] = gene_type
    return result


def annotate_gene_type(
    adata: ad.AnnData,
    te_bed: str,
    gene_tsv: str,
) -> None:
    """Annotate ``adata.var['gene_type']`` using BED/TSV references.

    Genes whose name appears in *te_bed* are labelled ``"TE"``.
    Remaining genes matching *gene_tsv* are labelled with their
    ``gene_type`` value (e.g. ``"protein_coding"``, ``"lncRNA"``).
    Unmatched genes are labelled ``"unknown"``.

    Args:
        adata: AnnData object whose ``var_names`` are gene symbols / TE names.
        te_bed: Path to the TE BED file (chrom, start, end, TE_name).
        gene_tsv: Path to the gene annotation TSV (gene_id, gene_name,
            gene_type).
    """
    logging.info("Parsing TE BED: %s", te_bed)
    te_set = _parse_te_bed(te_bed)
    logging.info("  TE BED: %d unique TE names", len(te_set))

    logging.info("Parsing gene TSV: %s", gene_tsv)
    gene_map = _parse_gene_tsv(gene_tsv)
    logging.info("  Gene TSV: %d entries", len(gene_map))

    gene_names = adata.var_names.tolist()
    types: List[str] = []
    n_te = 0
    n_gene = 0
    n_unknown = 0

    for g in gene_names:
        if g in te_set:
            types.append("TE")
            n_te += 1
        elif g in gene_map:
            types.append(gene_map[g] or "unknown")
            n_gene += 1
        else:
            types.append("unknown")
            n_unknown += 1

    adata.var["gene_type"] = types
    logging.info(
        "Gene type annotation: TE=%d, gene=%d, unknown=%d (total=%d)",
        n_te, n_gene, n_unknown, len(gene_names),
    )


def _filter_te(adata: ad.AnnData) -> ad.AnnData:
    """Remove TE genes from adata if ``gene_type`` column exists.

    Returns a subset with only non-TE genes.  If ``gene_type`` is absent
    the input is returned unchanged.
    """
    if "gene_type" not in adata.var.columns:
        logging.warning("gene_type not in var, skipping TE filter")
        return adata
    n_before = adata.n_vars
    mask = adata.var["gene_type"] != "TE"
    n_after = int(mask.sum())
    logging.info("TE filter: %d -> %d genes (%d TE removed)",
                 n_before, n_after, n_before - n_after)
    return adata[:, mask].copy()


# ---------------------------------------------------------------------------
# Cluster
# ---------------------------------------------------------------------------
def mode_cluster(
    adata: ad.AnnData,
    output: str,
    n_pcs: int = 50,
    n_neighbors: int = 50,
    resolution: float = 0.8,
    n_top_genes: int = 3000,
    batch_method: str = "harmony",
    batch_key: str = "",
    markers: str = "",
    plot_dir: str = "",
    auto_n_pcs: bool = False,
    skip_te: bool = False,
) -> None:
    """Cluster cells: preprocess → batch correct → neighbours → UMAP → Leiden.

    Steps:
        1. Normalise total counts per cell to 10 000 and log-transform.
        2. Store a ``raw`` snapshot (all genes, normalised + log1p) for
           downstream DEG / annotation.
        3. Select the top *n_top_genes* highly variable genes and subset.
        4. Scale HVGs (max_value=10).
        5. Run PCA (up to *n_pcs* components, or 100 if *auto_n_pcs*).
        6. (Optional) Auto-detect optimal n_pcs from the variance ratio
           sliding-window plateau.
        7. Apply batch correction (default: Harmony).
           - **BBKNN**: batch-balanced k-NN graph directly on PCA space.
           - **Harmony**: embed PCA coordinates in a batch-corrected space.
        8. Build a k-NN graph (*n_neighbors*, *n_pcs* PCs).
        9. Compute UMAP embedding (min_dist=0.1, spread=0.8).
        10. Leiden clustering at the given *resolution* (igraph flavour).
        11. Rank differentially expressed genes per cluster (Wilcoxon, use_raw).

    Args:
        adata: Input AnnData (raw counts expected in ``.X``).
        output: Path to write the clustered h5ad file.
        n_pcs: Number of principal components (default 50).
        n_neighbors: Number of neighbours for the k-NN graph (default 50).
        resolution: Leiden clustering resolution (default 0.8).
        n_top_genes: Number of highly variable genes to select (default 3000).
        batch_method: Batch correction method — ``"harmony"``, ``"bbknn"``,
            or ``""`` to skip. Default ``"harmony"``.
        batch_key: Column in ``adata.obs`` identifying batches. Falls back
            to ``"sample_id"`` then ``"batch"`` when empty.
        markers: Path to write a TSV of ranked marker genes. Empty string
            skips this step.
        plot_dir: Directory for cluster plots. Empty string disables plotting.
        auto_n_pcs: Auto-detect optimal n_pcs from PCA variance ratio
            sliding-window plateau (default False).
    """
    # 1. Normalise + log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 2. Save full-gene snapshot for downstream DEG / annotation
    adata.raw = adata.copy()

    # 2.5 Filter TE genes before HVG selection
    if skip_te:
        adata = _filter_te(adata)

    # 3. Detect HVGs (keep all genes for plotting, subset later)
    resolved_batch_key = batch_key or ("sample_id" if "sample_id" in adata.obs else "batch")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat",
        subset=False,
        batch_key=resolved_batch_key if resolved_batch_key in adata.obs else None,
    )

    # 4. Plot HVG (before subsetting — needs all genes as background)
    plotter = _make_plotter(plot_dir)
    if plotter:
        plotter.plot_hvg(adata, n_top_genes=n_top_genes)

    # 5. Subset to HVGs and scale
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)

    # 6. PCA — compute more components when auto-detecting
    pca_comps = max(n_pcs, 100) if auto_n_pcs else n_pcs
    pca_comps = min(pca_comps, max(2, adata.n_obs - 1))
    sc.tl.pca(adata, n_comps=pca_comps)

    # 6. Auto-detect n_pcs if requested
    recommended_n_pcs = n_pcs
    detect_diag: Dict = {}
    if auto_n_pcs:
        variance_ratio = adata.uns["pca"]["variance_ratio"]
        recommended_n_pcs, detect_diag = detect_n_pcs(variance_ratio)
        logging.info("Auto-detected n_pcs: %d", recommended_n_pcs)
        n_pcs = recommended_n_pcs

    # 7. Batch correction
    if batch_method == "bbknn":
        sc.external.pp.bbknn(adata, batch_key=resolved_batch_key)
    elif batch_method == "harmony":
        # Direct harmonypy call to avoid scanpy wrapper .T bug
        # (harmonypy >= 0.1.0 returns Z_corr as (n_cells, n_components),
        #  but old scanpy wrapper still transposes it)
        ho = hm.run_harmony(
            adata.obsm["X_pca"], adata.obs, resolved_batch_key,
        )
        Z = np.asarray(ho.Z_corr)
        if Z.ndim == 1:
            raise ValueError(
                f"harmonypy Z_corr is 1D shape={Z.shape}, expected 2D"
            )
        # Ensure (n_cells, n_components) orientation
        if Z.shape[0] != adata.n_obs:
            Z = Z.T
        adata.obsm["X_pca_harmony"] = Z
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                        use_rep="X_pca_harmony")

    # 8. k-NN graph (skip if already done by BBKNN or Harmony)
    if batch_method not in ("bbknn", "harmony"):
        sc.pp.neighbors(
            adata, n_neighbors=n_neighbors,
            n_pcs=n_pcs,
        )

    # 9. UMAP (tight params for clean clusters)
    sc.tl.umap(adata, min_dist=0.1, spread=0.8)

    # 10. Leiden clustering (igraph flavour)
    sc.tl.leiden(
        adata, resolution=resolution, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )

    # 11. DEG — use raw for full gene coverage
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

    if markers:
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(
            markers, sep="\t", index=False
        )

    if plotter:
        plotter.plot_pca_variance(
            adata,
            n_pcs=n_pcs,
            auto_n_pcs=auto_n_pcs,
            detect_diag=detect_diag if auto_n_pcs else None,
        )
        plotter.plot_cluster(adata, cluster_key="leiden", sample_key=resolved_batch_key)

    adata.write_h5ad(output)

    if markers:
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(
            markers, sep="\t", index=False
        )

    plotter = _make_plotter(plot_dir)
    if plotter:
        plotter.plot_cluster(adata, cluster_key="leiden", sample_key=resolved_batch_key)

    adata.write_h5ad(output)


# ---------------------------------------------------------------------------
# Cell type annotation
# ---------------------------------------------------------------------------
def mode_annotate(
    adata: ad.AnnData,
    output: str,
    marker_file: str = "",
    celltypist_model: str = "",
    llm_method: str = "",
    llm_model: str = "",
    llm_api_key: str = "",
    llm_base_url: str = "",
    llm_top_genes: int = 30,
    annotate_group: str = "",
    tissue: str = "",
    plot_dir: str = "",
) -> None:
    """Cell type annotation dispatcher.

    Runs one or more annotation strategies in sequence:
        1. **Marker-based** — score cells against a TSV of known markers
           (requires *marker_file*).
        2. **CellTypist** — automated annotation with a pre-trained
           CellTypist model (requires *celltypist_model*).
        3. **LLM-assisted** — send top DEGs per cluster to an LLM
           (OpenAI, Ollama, or file) for cell type prediction
           (requires *llm_method*).

    Args:
        adata: Input AnnData with clustering results (e.g. ``leiden``).
        output: Path to write the annotated h5ad file.
        marker_file: Path to a TSV with ``cell_type`` and ``markers`` columns.
            Empty string skips marker-based annotation.
        celltypist_model: CellTypist model name (e.g. ``"Immune_All_Low"``).
            Empty string skips CellTypist annotation.
        llm_method: LLM backend — ``"openai"``, ``"ollama"``, or ``"file"``.
            Empty string skips LLM annotation.
        llm_model: Model identifier for the LLM backend
            (e.g. ``"gpt-4o"`` or ``"llama3.1"``).
        llm_api_key: API key for the OpenAI backend (ignored for others).
        llm_base_url: Base URL for the LLM API endpoint.
        llm_top_genes: Number of top DEGs per cluster to include in the
            LLM prompt (default 30).
        annotate_group: Column in ``adata.obs`` to group clusters by.
            Falls back to ``"leiden"`` when empty.
        tissue: Tissue name included in the LLM prompt for context.
        plot_dir: Directory for annotation plots. Empty string disables
            plotting.
    """
    # Store tissue for downstream access
    if tissue:
        adata.uns["tissue"] = tissue

    if marker_file:
        _annotate_markers(adata, marker_file=marker_file)
    if celltypist_model:
        _annotate_celltypist(adata, celltypist_model=celltypist_model)
    if llm_method:
        _annotate_llm(
            adata,
            llm_method=llm_method,
            llm_model=llm_model,
            llm_api_key=llm_api_key,
            llm_base_url=llm_base_url,
            llm_top_genes=llm_top_genes,
            annotate_group=annotate_group,
            tissue=tissue,
            output=output,
        )

    plotter = _make_plotter(plot_dir)
    if plotter:
        # Build explicit list of annotation columns that were created
        anno_keys: List[str] = []
        if marker_file:
            anno_keys.append("cell_type")
        if celltypist_model:
            anno_keys.append("celltypist_label")
        if llm_method:
            anno_keys.append("llm_label")

        plotter.plot_annotate(
            adata,
            marker_file=marker_file,
            annotate_group=annotate_group or "leiden",
            annotation_keys=anno_keys,
            score_col="celltypist_score",
            has_rank_genes="rank_genes_groups" in adata.uns,
        )

    adata.write_h5ad(output)


def _annotate_markers(adata: ad.AnnData, marker_file: str) -> None:
    """Score cells against known marker gene sets and assign cell types.

    Reads a tab-separated file with columns ``cell_type`` and ``markers``
    (comma-separated gene names).  For each cell type, computes a
    ``score_genes`` enrichment score, then assigns each cell the cell type
    with the highest score.

    Args:
        adata: AnnData object to annotate (modified in place).
        marker_file: Path to the marker gene TSV file.
    """
    marker_genes: Dict[str, List[str]] = {}
    with open(marker_file, encoding="utf-8") as fh:
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


def _annotate_celltypist(adata: ad.AnnData, celltypist_model: str) -> None:
    """Annotate cell types using a pre-trained CellTypist model.

    The function normalises and log-transforms a copy of the data (using
    the ``counts`` layer if available), runs CellTypist with majority
    voting, and writes ``celltypist_label`` and ``celltypist_score``
    columns into ``adata.obs``.

    Args:
        adata: AnnData object to annotate (modified in place).
        celltypist_model: Name of the CellTypist model to load
            (e.g. ``"Immune_All_Low"``).
    """
    import celltypist
    from celltypist import models

    adata_ct = adata.copy()
    if "counts" in adata_ct.layers:
        adata_ct.X = adata_ct.layers["counts"]
    sc.pp.normalize_per_cell(adata_ct, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_ct)
    adata_ct.X = adata_ct.X.toarray() if hasattr(adata_ct.X, "toarray") else adata_ct.X

    models.download_models(force_update=True, model=[celltypist_model])
    model = models.Model.load(model=celltypist_model)
    predictions = celltypist.annotate(adata_ct, model=model, majority_voting=True)
    pred_adata = predictions.to_adata()
    adata.obs["celltypist_label"] = pred_adata.obs.loc[adata.obs.index, "majority_voting"]
    adata.obs["celltypist_score"] = pred_adata.obs.loc[adata.obs.index, "conf_score"]


def _annotate_llm(
    adata: ad.AnnData,
    llm_method: str,
    llm_model: str,
    llm_api_key: str,
    llm_base_url: str,
    llm_top_genes: int,
    annotate_group: str,
    tissue: str,
    output: str,
) -> None:
    """Annotate cell clusters by sending top DEGs to a large language model.

    For each cluster (defined by *annotate_group*), the top
    *llm_top_genes* differentially expressed genes are extracted and
    formatted into a structured prompt.  The prompt is then sent to the
    chosen LLM backend (OpenAI, Ollama, or written to a file for manual
    use).  The returned JSON mapping is stored in ``adata.uns["llm_annotation"]``
    and mapped to ``adata.obs["llm_label"]``.

    Args:
        adata: AnnData object with ``rank_genes_groups`` computed.
        llm_method: LLM backend — ``"openai"``, ``"ollama"``, or ``"file"``.
        llm_model: Model identifier (e.g. ``"gpt-4o"``, ``"llama3.1"``).
        llm_api_key: API key for the OpenAI backend.
        llm_base_url: Base URL for the LLM API endpoint.
        llm_top_genes: Number of top DEGs per cluster to include in the prompt.
        annotate_group: Column in ``adata.obs`` defining cluster membership.
            Falls back to ``"leiden"`` when empty.
        tissue: Tissue name for context in the prompt.
        output: h5ad output path — used to derive the file name for
            ``"file"`` method prompt output.
    """
    group_key = annotate_group or "leiden"
    if group_key not in adata.obs.columns:
        logging.warning("Group key '%s' not found in obs, skipping LLM annotation", group_key)
        return

    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, group_key, method="wilcoxon")

    n_genes = llm_top_genes
    cluster_markers: Dict[str, List[str]] = {}
    for cluster in adata.obs[group_key].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        cluster_markers[str(cluster)] = df.head(n_genes)["names"].tolist()

    resolved_tissue = adata.uns.get("tissue", tissue or "unknown tissue")
    annotation = _llm_annotate_clusters(
        cluster_markers, resolved_tissue,
        llm_method=llm_method,
        llm_model=llm_model,
        llm_api_key=llm_api_key,
        llm_base_url=llm_base_url,
        output=output,
    )

    label_map = {k: v.get("cell_type", "unknown") for k, v in annotation.items()}
    adata.obs["llm_label"] = adata.obs[group_key].map(label_map).astype("category")
    adata.uns["llm_annotation"] = annotation

    logging.info("LLM annotation complete: %d clusters annotated", len(annotation))
    for cluster, info in annotation.items():
        logging.info("  Cluster %s: %s (confidence: %s)",
                     cluster, info.get("cell_type", "?"), info.get("confidence", "?"))


def _llm_annotate_clusters(
    cluster_markers: Dict[str, List[str]],
    tissue: str,
    llm_method: str,
    llm_model: str,
    llm_api_key: str,
    llm_base_url: str,
    output: str,
) -> dict:
    """Dispatch the annotation prompt to the chosen LLM backend.

    Args:
        cluster_markers: Mapping of cluster ID → list of top marker gene names.
        tissue: Tissue name for context in the prompt.
        llm_method: ``"openai"``, ``"ollama"``, or ``"file"``.
        llm_model: Model identifier for the API call.
        llm_api_key: API key for the OpenAI backend.
        llm_base_url: Base URL for the LLM API.
        output: h5ad output path — used to derive the prompt file name
            when *llm_method* is ``"file"``.

    Returns:
        Parsed JSON dict mapping cluster IDs to annotation info, or an
        empty dict on failure / file-only mode.
    """
    prompt = _build_annotation_prompt(cluster_markers, tissue)

    if llm_method == "openai":
        return _call_openai(prompt, llm_model=llm_model, llm_api_key=llm_api_key, llm_base_url=llm_base_url)
    elif llm_method == "ollama":
        return _call_ollama(prompt, llm_model=llm_model, llm_base_url=llm_base_url)
    elif llm_method == "file":
        output_file = output.replace(".h5ad", "_llm_prompt.txt")
        with open(output_file, "w", encoding="utf-8") as fh:
            fh.write(prompt)
        logging.info("LLM prompt written to %s", output_file)
        return {}
    else:
        logging.warning("Unknown LLM method: %s", llm_method)
        return {}


def _build_annotation_prompt(cluster_markers: Dict[str, List[str]], tissue: str) -> str:
    """Build a structured prompt for LLM-based cell type annotation.

    The prompt instructs the LLM to return a JSON object keyed by cluster
    ID, with ``cell_type``, ``confidence``, ``reasoning``, and
    ``canonical_markers`` for each cluster.

    Args:
        cluster_markers: Mapping of cluster ID → list of top marker genes.
        tissue: Tissue name for biological context.

    Returns:
        The fully formatted prompt string.
    """
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


def _call_openai(
    prompt: str,
    llm_model: str,
    llm_api_key: str,
    llm_base_url: str,
) -> dict:
    """Call the OpenAI-compatible chat completions API for cell type annotation.

    Args:
        prompt: The fully formatted annotation prompt.
        llm_model: Model identifier (e.g. ``"gpt-4o"``).
        llm_api_key: API key for authentication.
        llm_base_url: Base URL for the API endpoint.

    Returns:
        Parsed JSON dict of cluster annotations, or an empty dict on failure.
    """
    try:
        import openai
        client = openai.OpenAI(api_key=llm_api_key, base_url=llm_base_url)
        response = client.chat.completions.create(
            model=llm_model or "gpt-4o",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
            response_format={"type": "json_object"},
        )
        return json.loads(response.choices[0].message.content)
    except Exception as exc:
        logging.error("OpenAI API call failed: %s", exc)
        return {}


def _call_ollama(
    prompt: str,
    llm_model: str,
    llm_base_url: str,
) -> dict:
    """Call a local Ollama server for cell type annotation.

    Args:
        prompt: The fully formatted annotation prompt.
        llm_model: Ollama model name (e.g. ``"llama3.1"``).
        llm_base_url: Ollama server URL (default ``http://localhost:11434``).

    Returns:
        Parsed JSON dict of cluster annotations, or an empty dict on failure.
    """
    try:
        import requests
        url = llm_base_url or "http://localhost:11434"
        resp = requests.post(
            f"{url}/api/chat",
            json={
                "model": llm_model or "llama3.1",
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
# PubMed literature search
# ---------------------------------------------------------------------------
def _search_pubmed(
    gene: str,
    cell_type: str,
    tissue: str = "",
    max_results: int = 2,
) -> List[Dict[str, str]]:
    """Search PubMed for references supporting a gene as a cell type marker.

    Uses NCBI E-utilities (no API key, rate-limited to 3 req/s).

    Args:
        gene: Gene symbol to search for.
        cell_type: Cell type the gene is a marker for.
        tissue: Tissue context (optional, narrows search).
        max_results: Maximum number of references to return.

    Returns:
        List of dicts with keys ``pmid``, ``title``, ``year``.
    """
    import urllib.parse
    import urllib.request
    import xml.etree.ElementTree as ET

    query_parts = [f"{gene}", f"{cell_type}"]
    if tissue:
        query_parts.append(tissue)
    query = " AND ".join(query_parts)

    results: List[Dict[str, str]] = []
    try:
        # Step 1: esearch
        params = urllib.parse.urlencode({
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "sort": "relevance",
            "retmode": "json",
        })
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        pmids = data.get("esearchresult", {}).get("idlist", [])
        if not pmids:
            return results

        time.sleep(0.35)  # rate limit

        # Step 2: efetch to get titles
        fetch_params = urllib.parse.urlencode({
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
        })
        fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{fetch_params}"
        with urllib.request.urlopen(fetch_url, timeout=10) as resp:
            xml_text = resp.read().decode()

        root = ET.fromstring(xml_text)
        for article in root.findall(".//PubmedArticle"):
            pmid_el = article.find(".//PMID")
            title_el = article.find(".//ArticleTitle")
            year_el = article.find(".//PubDate/Year")
            if pmid_el is not None:
                pmid = pmid_el.text or ""
                title = title_el.text if title_el is not None else ""
                year = year_el.text if year_el is not None else ""
                # Handle PubDate without Year (use MedlineDate)
                if not year:
                    medline_date = article.find(".//PubDate/MedlineDate")
                    if medline_date is not None and medline_date.text:
                        year = medline_date.text[:4]
                results.append({"pmid": pmid, "title": title[:120], "year": year})

    except Exception as exc:
        logging.warning("PubMed search failed for %s: %s", gene, exc)

    return results


def _build_auto_annotation_prompt(
    cluster_id: str,
    top_genes: List[str],
    tissue: str,
    n_cells: int,
    mean_genes: float,
    mean_counts: float,
    pct_mt: float,
    other_annotations: Optional[Dict[str, Dict]] = None,
    umap_distances: Optional[Dict[str, float]] = None,
) -> str:
    """Build prompt for autonomous cell type annotation of a single cluster.

    The LLM must:
        1. Identify the most likely cell type.
        2. List key marker genes (5-10) that support the identification.
        3. Provide reasoning.
        4. Rate confidence (high/medium/low).
        5. Flag if the cluster looks like low quality / contamination.
        6. Consider if this cluster might be a subcluster of another cell type.

    Args:
        cluster_id: ID of the cluster being annotated.
        top_genes: List of top differentially expressed genes.
        tissue: Tissue name for context.
        n_cells: Number of cells in the cluster.
        mean_genes: Mean genes per cell.
        mean_counts: Mean UMI per cell.
        pct_mt: Mean mitochondrial percentage.
        other_annotations: Dict of other cluster_id -> annotation results.
            Used to help LLM identify if current cluster is a subcluster.
        umap_distances: Dict of other_cluster_id -> UMAP distance to current cluster.
            Used to help LLM judge spatial relationships.
    """
    gene_list = ", ".join(top_genes[:50])

    # Build context about other clusters
    other_context = ""
    if other_annotations:
        other_lines = []
        for cid, ann in other_annotations.items():
            ct = ann.get("cell_type", "Unknown")
            markers = ", ".join(ann.get("key_markers", [])[:5])
            conf = ann.get("confidence", "unknown")
            dist_info = ""
            if umap_distances and cid in umap_distances:
                dist = umap_distances[cid]
                dist_info = f" (UMAP distance: {dist:.1f})"
            other_lines.append(f"  - Cluster {cid}: {ct} (confidence: {conf}, key markers: {markers}){dist_info}")
        other_context = "\n".join(other_lines)

    prompt = f"""You are an expert single-cell RNA-seq analyst. Annotate the following cluster from a {tissue} dataset.

## Cluster {cluster_id}
- Cells: {n_cells}
- Mean genes/cell: {mean_genes:.0f}
- Mean UMI/cell: {mean_counts:.0f}
- MT%: {pct_mt:.1f}%
- Top 50 DEGs (ranked by Wilcoxon): {gene_list}
"""

    if other_context:
        prompt += f"""
## Other clusters already annotated
{other_context}
"""

    prompt += """
## Task
Identify the cell type for this cluster. Consider:
1. **Cell type identification**: What cell type do the marker genes suggest?
2. **Subcluster check**: Is this cluster a subcluster (subtype) of another cluster? Look at marker gene similarity and UMAP distance.
3. **Quality check**: Are there quality issues (low genes, low UMI, high MT%)?
4. **Resolution check**: If markers are very similar to another cluster, this might be over-clustering.

Output a JSON object with:
- "cell_type": predicted cell type name (use standard nomenclature, e.g. "Stromal", "Smooth_muscle", "Endothelial", "Macrophage", "T_cell", "Unknown", etc.)
- "key_markers": list of 5-10 genes from the DEGs that best support this identification
- "reasoning": brief explanation of your reasoning, including:
  * Why this cell type?
  * Is this a subcluster of another cluster? If so, which one?
  * Are there quality concerns?
  * Should this cluster be merged with another cluster?
- "confidence": "high", "medium", or "low"
- "quality_flag": null if cluster is normal, or a string describing the problem
    - "low_quality": mean_genes < 800 or mean_counts < 3000 or MT% > 20
    - "unannotated": top markers are ENSMMUG (unannotated macaque genes)
    - "ribosomal": top markers are RPS/RPL (ribosomal contamination)
    - "te_dominated": top markers are TE elements (ERVK, Alu, etc.)
    - "unknown": cannot determine cell type from markers
- "is_subcluster": boolean - true if this is likely a subcluster of another cell type
- "parent_cluster": cluster_id of the parent cluster if is_subcluster is true, else null
- "should_merge": boolean - true if this cluster should be merged with parent_cluster

## Rules
- ALL cell type annotations MUST have PubMed evidence. Include PMIDs in your reasoning. If no published evidence exists for your predicted cell type, set cell_type to "Unverified" and explain why.
- If most top genes start with "ENSMMUG", set cell_type to "Unannotated" and quality_flag to "unannotated"
- If most top genes are MT- or mt-prefixed, set quality_flag to "low_quality"
- If most top genes are RPS/RPL, set quality_flag to "ribosomal"
- If top genes are TE elements (Alu, ERVK, LINE, SINE, LTR), you MUST search for PubMed evidence. If you find published literature supporting this TE-high population as a real biological group (with PMID), annotate accordingly (e.g. "Alu_high"). If no public evidence exists, set cell_type to "Unverified_TE" and quality_flag to "te_dominated".
- If you truly cannot identify the cell type, set cell_type to "Unknown" and quality_flag to "unknown"
- If marker genes are very similar to another cluster AND UMAP distance is close (< 5.0), consider setting is_subcluster=true and should_merge=true
- If marker genes are similar but UMAP distance is far (> 10.0), this might be over-clustering - still set is_subcluster=true and should_merge=true
- If marker genes are different even though UMAP is close, these are likely distinct cell types
- Output ONLY valid JSON, no markdown
"""
    return prompt


def _ai_annotate_cluster(
    cluster_id: str,
    top_genes: List[str],
    tissue: str,
    n_cells: int,
    mean_genes: float,
    mean_counts: float,
    pct_mt: float,
    llm_method: str,
    llm_model: str,
    llm_api_key: str,
    llm_base_url: str,
    other_annotations: Optional[Dict[str, Dict]] = None,
    umap_distances: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """AI-annotate a single cluster: LLM decision + PubMed references.

    Args:
        cluster_id: ID of the cluster being annotated.
        top_genes: List of top differentially expressed genes.
        tissue: Tissue name for context.
        n_cells: Number of cells in the cluster.
        mean_genes: Mean genes per cell.
        mean_counts: Mean UMI per cell.
        pct_mt: Mean mitochondrial percentage.
        llm_method: LLM backend ("openai" or "ollama").
        llm_model: Model identifier.
        llm_api_key: API key for OpenAI backend.
        llm_base_url: Base URL for LLM API.
        other_annotations: Dict of other cluster_id -> annotation results.
            Used to help LLM identify if current cluster is a subcluster.
        umap_distances: Dict of other_cluster_id -> UMAP distance to current cluster.
            Used to help LLM judge spatial relationships.

    Returns dict with keys: cell_type, key_markers, reasoning, confidence,
    quality_flag, is_subcluster, parent_cluster, should_merge, references.
    """
    prompt = _build_auto_annotation_prompt(
        cluster_id, top_genes, tissue, n_cells, mean_genes, mean_counts, pct_mt,
        other_annotations=other_annotations,
        umap_distances=umap_distances,
    )

    # Call LLM
    if llm_method == "openai":
        raw = _call_openai(prompt, llm_model=llm_model, llm_api_key=llm_api_key, llm_base_url=llm_base_url)
    elif llm_method == "ollama":
        raw = _call_ollama(prompt, llm_model=llm_model, llm_base_url=llm_base_url)
    else:
        logging.warning("Unknown LLM method '%s' for cluster %s", llm_method, cluster_id)
        return {"cell_type": "Unknown", "key_markers": [], "reasoning": "", "confidence": "low",
                "quality_flag": "unknown", "references": {}}

    # Parse result
    annotation: Dict[str, Any] = {
        "cell_type": raw.get("cell_type", "Unknown"),
        "key_markers": raw.get("key_markers", []),
        "reasoning": raw.get("reasoning", ""),
        "confidence": raw.get("confidence", "low"),
        "quality_flag": raw.get("quality_flag"),
        "is_subcluster": raw.get("is_subcluster", False),
        "parent_cluster": raw.get("parent_cluster"),
        "should_merge": raw.get("should_merge", False),
        "references": {},
    }

    # Search PubMed for each key marker
    cell_type = annotation["cell_type"]
    if cell_type not in ("Unknown", "Unannotated", "ERVK_high", "Alu_high"):
        for gene in annotation["key_markers"][:10]:
            refs = _search_pubmed(gene, cell_type, tissue, max_results=2)
            if refs:
                annotation["references"][gene] = refs
            time.sleep(0.35)

    return annotation


# ---------------------------------------------------------------------------
# Cell type separation check
# ---------------------------------------------------------------------------
def _check_cell_type_separation(
    adata: ad.AnnData,
    cell_type_key: str = "cell_type",
    cluster_key: str = "leiden",
    distance_threshold: float = 5.0,
    min_cells_per_cluster: int = 50,
    annotations: Optional[Dict[str, Dict]] = None,
) -> Tuple[bool, Dict[str, Any]]:
    """Check if same cell_type clusters are separated too far in UMAP space.

    When the same cell type appears in multiple clusters that are far apart
    in UMAP space, it may indicate:
    - Over-clustering (resolution too high) → need to merge
    - Annotation error → low-quality cluster mislabeled → need to flag/filter

    Algorithm:
        1. Compute UMAP centroid for each cluster.
        2. For cell types with multiple clusters, compute pairwise distances.
        3. If max distance > threshold, mark as "separated".
        4. Compare key_markers between separated clusters:
           - High overlap (Jaccard > 0.3) → over-clustering → suggest resolution
           - Low overlap (Jaccard ≤ 0.3) → annotation error → flag low-quality cluster
        5. Suggest resolution adjustment based on separation severity.

    Args:
        adata: AnnData with UMAP coordinates, leiden clusters, and cell_type.
        cell_type_key: Column name for cell type annotations.
        cluster_key: Column name for cluster labels.
        distance_threshold: Max allowed UMAP distance between same-type clusters.
        min_cells_per_cluster: Min cells to consider a cluster valid.
        annotations: Dict of cluster_id -> annotation results with key_markers.

    Returns:
        Tuple of (needs_reclustering, diagnostics) where diagnostics is a dict
        with keys: separated_types, cluster_distances, suggested_resolution.
    """
    if "X_umap" not in adata.obsm:
        logging.warning("No UMAP coordinates found, skipping separation check")
        return False, {}

    if cell_type_key not in adata.obs.columns or cluster_key not in adata.obs.columns:
        return False, {}

    umap_coords = adata.obsm["X_umap"]
    obs = adata.obs

    # Step 1: Compute UMAP centroid for each cluster
    cluster_centroids = {}
    cluster_sizes = {}
    for cluster in obs[cluster_key].unique():
        mask = obs[cluster_key] == cluster
        n_cells = int(mask.sum())
        if n_cells >= min_cells_per_cluster:
            centroid = umap_coords[mask.values].mean(axis=0)
            cluster_centroids[cluster] = centroid
            cluster_sizes[cluster] = n_cells

    # Step 2: Group clusters by cell type
    type_to_clusters = {}
    for cluster, centroid in cluster_centroids.items():
        ct = obs.loc[obs[cluster_key] == cluster, cell_type_key].iloc[0]
        if ct not in type_to_clusters:
            type_to_clusters[ct] = []
        type_to_clusters[ct].append({
            "cluster": cluster,
            "centroid": centroid,
            "n_cells": cluster_sizes[cluster],
        })

    # Step 3: Check separation for each cell type with multiple clusters
    separated_types = []  # over-clustering (high marker overlap)
    misannotated = []     # annotation error (low marker overlap)
    cluster_distances = {}

    for ct, clusters_info in type_to_clusters.items():
        if len(clusters_info) < 2:
            continue  # Only one cluster, no separation issue

        # Compute pairwise distances between clusters of same cell type
        max_dist = 0.0
        pair_info = []
        for i in range(len(clusters_info)):
            for j in range(i + 1, len(clusters_info)):
                ci = clusters_info[i]
                cj = clusters_info[j]
                dist = float(np.linalg.norm(ci["centroid"] - cj["centroid"]))
                pair_info.append({
                    "cluster_i": ci["cluster"],
                    "cluster_j": cj["cluster"],
                    "distance": round(dist, 2),
                    "cells_i": ci["n_cells"],
                    "cells_j": cj["n_cells"],
                })
                max_dist = max(max_dist, dist)

        cluster_distances[ct] = {
            "n_clusters": len(clusters_info),
            "max_distance": round(max_dist, 2),
            "pairs": pair_info,
        }

        if max_dist > distance_threshold:
            # Compare marker overlap between the most distant clusters
            # to distinguish over-clustering from annotation error
            if annotations:
                all_markers = []
                for ci in clusters_info:
                    markers = set(annotations.get(ci["cluster"], {}).get("key_markers", []))
                    all_markers.append(markers)

                # Compute average pairwise Jaccard of top markers
                jaccards = []
                for i in range(len(all_markers)):
                    for j in range(i + 1, len(all_markers)):
                        union = all_markers[i] | all_markers[j]
                        inter = all_markers[i] & all_markers[j]
                        j = len(inter) / len(union) if union else 0.0
                        jaccards.append(j)
                avg_jaccard = sum(jaccards) / len(jaccards) if jaccards else 0.0

                if avg_jaccard > 0.3:
                    # High marker overlap → genuine over-clustering
                    separated_types.append({
                        "cell_type": ct,
                        "n_clusters": len(clusters_info),
                        "max_distance": round(max_dist, 2),
                        "marker_jaccard": round(avg_jaccard, 3),
                        "reason": "over_clustering",
                    })
                else:
                    # Low marker overlap → annotation error, flag smallest cluster
                    smallest = min(clusters_info, key=lambda x: x["n_cells"])
                    misannotated.append({
                        "cell_type": ct,
                        "n_clusters": len(clusters_info),
                        "max_distance": round(max_dist, 2),
                        "marker_jaccard": round(avg_jaccard, 3),
                        "flagged_cluster": smallest["cluster"],
                        "flagged_cells": smallest["n_cells"],
                        "reason": "low_marker_overlap",
                    })
            else:
                # No annotations available, treat as over-clustering
                separated_types.append({
                    "cell_type": ct,
                    "n_clusters": len(clusters_info),
                    "max_distance": round(max_dist, 2),
                })

    # Step 4: Determine if re-clustering is needed
    needs_reclustering = len(separated_types) > 0

    # Suggest resolution adjustment
    suggested_resolution = None
    if needs_reclustering:
        # If many cell types are separated, suggest lower resolution
        separation_ratio = len(separated_types) / len(type_to_clusters)
        if separation_ratio > 0.3:  # >30% of cell types are separated
            suggested_resolution = 0.2  # More aggressive merging
        elif separation_ratio > 0.1:  # >10% separated
            suggested_resolution = 0.25
        else:
            suggested_resolution = 0.3  # Keep current or slight adjustment

    diagnostics = {
        "separated_types": separated_types,
        "misannotated": misannotated,
        "cluster_distances": cluster_distances,
        "n_cell_types": len(type_to_clusters),
        "n_separated": len(separated_types),
        "n_misannotated": len(misannotated),
        "suggested_resolution": suggested_resolution,
    }

    if needs_reclustering:
        logging.info("Cell type separation check: %d over-clustering, %d misannotated",
                     len(separated_types), len(misannotated))
        for st in separated_types:
            logging.info("  %s: %d clusters, max dist=%.2f, marker Jaccard=%.3f (over-clustering)",
                         st["cell_type"], st["n_clusters"], st["max_distance"], st.get("marker_jaccard", 0))
    if misannotated:
        for ma in misannotated:
            logging.info("  %s: %d clusters, max dist=%.2f, marker Jaccard=%.3f -> FLAG cluster %s (%d cells, low quality)",
                         ma["cell_type"], ma["n_clusters"], ma["max_distance"],
                         ma["marker_jaccard"], ma["flagged_cluster"], ma["flagged_cells"])

    return needs_reclustering, diagnostics


# ---------------------------------------------------------------------------
# Cluster quality analysis
# ---------------------------------------------------------------------------
def _analyze_cluster_quality(
    adata: ad.AnnData,
    cluster: str,
    annotation: Dict[str, Any],
    min_genes: int = 800,
    min_counts: int = 3000,
    max_pct_mt: float = 20.0,
) -> Dict[str, Any]:
    """Analyze quality of a single cluster and determine if it needs filtering.

    Combines AI-reported quality_flag with programmatic QC checks.

    Returns dict with: cluster, n_cells, mean_genes, mean_counts, pct_mt,
    ai_cell_type, ai_confidence, ai_quality_flag, flags (list), should_filter.
    """
    mask = adata.obs["leiden"] == cluster
    obs = adata.obs[mask]
    n_cells = int(mask.sum())

    mean_genes = float(obs["n_genes_by_counts"].mean()) if "n_genes_by_counts" in obs.columns else 0.0
    mean_counts = float(obs["total_counts"].mean()) if "total_counts" in obs.columns else 0.0
    pct_mt = float(obs["pct_counts_mt"].mean()) if "pct_counts_mt" in obs.columns else 0.0

    flags: List[str] = []

    # Programmatic checks
    if mean_genes < min_genes:
        flags.append("low_genes")
    if mean_counts < min_counts:
        flags.append("low_counts")
    if pct_mt > max_pct_mt:
        flags.append("high_mt")

    # AI-reported quality flag
    ai_flag = annotation.get("quality_flag")
    if ai_flag:
        # te_dominated alone should NOT trigger filtering
        # (scTE data has real ERVK_high/Alu_high populations)
        if ai_flag != "te_dominated":
            flags.append(str(ai_flag))
        else:
            # Only flag te_dominated if combined with other quality issues
            has_quality_issues = any(f in flags for f in ["low_genes", "low_counts", "high_mt"])
            if has_quality_issues:
                flags.append(str(ai_flag))

    # Confidence check
    confidence = annotation.get("confidence", "low")
    cell_type = annotation.get("cell_type", "Unknown")

    # Should filter: any flag present, or Unknown with low confidence
    should_filter = bool(flags) or (cell_type == "Unknown" and confidence == "low")

    report = {
        "cluster": cluster,
        "n_cells": n_cells,
        "mean_genes": round(mean_genes, 1),
        "mean_counts": round(mean_counts, 1),
        "pct_mt": round(pct_mt, 2),
        "ai_cell_type": cell_type,
        "ai_confidence": confidence,
        "ai_quality_flag": ai_flag,
        "flags": flags,
        "should_filter": should_filter,
    }
    return report


def _filter_flagged_cells(
    adata: ad.AnnData,
    quality_reports: List[Dict[str, Any]],
    min_genes: int = 800,
    min_counts: int = 3000,
    max_pct_mt: float = 20.0,
) -> Tuple[ad.AnnData, int]:
    """Filter individual cells within flagged clusters (NOT remove clusters).

    For each flagged cluster, removes cells that fail QC thresholds.
    Cells in non-flagged clusters are kept as-is.

    Args:
        adata: Clustered AnnData.
        quality_reports: List of dicts from _analyze_cluster_quality.
        min_genes: Min genes per cell threshold.
        min_counts: Min UMI per cell threshold.
        max_pct_mt: Max MT% threshold.

    Returns:
        Tuple of (filtered_adata, n_removed).
    """
    flagged_clusters = {r["cluster"] for r in quality_reports if r["should_filter"]}
    if not flagged_clusters:
        return adata, 0

    logging.info("Filtering cells in %d flagged clusters: %s",
                 len(flagged_clusters), sorted(flagged_clusters))

    keep_mask = pd.Series(True, index=adata.obs.index)
    total_removed = 0

    for cluster in flagged_clusters:
        cluster_mask = adata.obs["leiden"] == cluster
        cluster_cells = adata.obs[cluster_mask]

        cell_qc = pd.Series(True, index=cluster_cells.index)
        if "n_genes_by_counts" in adata.obs.columns:
            cell_qc &= cluster_cells["n_genes_by_counts"] >= min_genes
        if "total_counts" in adata.obs.columns:
            cell_qc &= cluster_cells["total_counts"] >= min_counts
        if "pct_counts_mt" in adata.obs.columns:
            cell_qc &= cluster_cells["pct_counts_mt"] <= max_pct_mt

        bad_cells = cluster_cells[~cell_qc].index
        keep_mask[bad_cells] = False
        n_removed = len(bad_cells)
        total_removed += n_removed
        logging.info("  Cluster %s: removed %d / %d cells",
                     cluster, n_removed, len(cluster_cells))

    adata_filtered = adata[keep_mask.values].copy()
    logging.info("Total: %d -> %d cells (%d removed)",
                 adata.n_obs, adata_filtered.n_obs, total_removed)

    return adata_filtered, total_removed


def _write_annotation_report(
    annotations: Dict[str, Dict],
    quality_reports: List[Dict[str, Any]],
    output_path: str,
) -> None:
    """Write per-cluster annotation report as TSV.

    Columns: cluster, cell_type, confidence, key_markers, n_refs,
    flags, should_filter, is_subcluster, parent_cluster, should_merge,
    reasoning.
    """
    rows = []
    for cluster in sorted(annotations.keys(), key=lambda x: int(x) if x.isdigit() else x):
        ann = annotations[cluster]
        qr = next((r for r in quality_reports if r["cluster"] == cluster), {})
        n_refs = sum(len(v) for v in ann.get("references", {}).values())

        # Get merge info if available
        merge_info = ann.get("merge_info", {})
        parent_cluster = merge_info.get("parent_cluster", ann.get("parent_cluster"))
        parent_cell_type = merge_info.get("parent_cell_type")

        rows.append({
            "cluster": cluster,
            "cell_type": ann.get("cell_type", "Unknown"),
            "confidence": ann.get("confidence", "?"),
            "key_markers": ",".join(ann.get("key_markers", [])),
            "n_refs": n_refs,
            "flags": ",".join(qr.get("flags", [])),
            "should_filter": qr.get("should_filter", False),
            "is_subcluster": ann.get("is_subcluster", False),
            "parent_cluster": parent_cluster,
            "parent_cell_type": parent_cell_type,
            "should_merge": ann.get("should_merge", False),
            "reasoning": ann.get("reasoning", "")[:200],
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep="\t", index=False)
    logging.info("Annotation report: %s", output_path)


def _write_references_report(
    annotations: Dict[str, Dict],
    output_path: str,
) -> None:
    """Write per-marker reference report as TSV.

    Columns: cluster, cell_type, gene, pmid, title, year.
    """
    rows = []
    for cluster in sorted(annotations.keys(), key=lambda x: int(x) if x.isdigit() else x):
        ann = annotations[cluster]
        ct = ann.get("cell_type", "Unknown")
        for gene, refs in ann.get("references", {}).items():
            for ref in refs:
                rows.append({
                    "cluster": cluster,
                    "cell_type": ct,
                    "gene": gene,
                    "pmid": ref.get("pmid", ""),
                    "title": ref.get("title", ""),
                    "year": ref.get("year", ""),
                })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep="\t", index=False)
    logging.info("References report: %s (%d entries)", output_path, len(rows))


# ---------------------------------------------------------------------------
# Auto mode: iterative cluster → AI annotate → QC → filter → re-cluster
# ---------------------------------------------------------------------------
def _generate_audit_report(
    ctx: Dict[str, Any],
    report_dir: str,
    llm_method: str,
    llm_model: str,
    llm_api_key: str,
    llm_base_url: str,
) -> None:
    """Call LLM to generate audit report and reproducible script.

    Args:
        ctx: Full iteration context collected during mode_auto.
        report_dir: Directory to write report files.
        llm_method/llm_model/llm_api_key/llm_base_url: LLM config.
    """
    prompt = (
        "你是一位生物信息学审计员。根据以下scRNA-seq auto-mode的执行上下文，"
        "生成两个文件：\n\n"
        "1. audit_report.md — 中文审计报告，包含：\n"
        "   - 流程参数（resolution, batch method, QC阈值等）\n"
        "   - 每轮迭代摘要：聚类结果、细胞类型注释（含置信度和推理依据）、QC标记、过滤决策\n"
        "   - 最终结果：细胞/基因数、迭代次数、停止原因\n"
        "   - 关键决策及其理由\n"
        "   - 注意：所有细胞类型注释必须有PubMed文献PMID支撑。"
        "无法找到文献支撑的注释应标记为'未验证'。\n\n"
        "2. decision_log.sh — LLM决策日志（不是CLI重跑脚本），记录：\n"
        "   - 每轮迭代中LLM做了哪些关键决策\n"
        "   - 标记了哪些cluster、理由是什么\n"
        "   - 过滤了哪些细胞、为什么\n"
        "   - 分辨率是否调整、为什么\n"
        "   - 格式：每条决策一行注释 # [iter N] 决策内容\n\n"
        "返回JSON对象，两个key：\n"
        '  "audit_report": <markdown字符串>,\n'
        '  "decision_log": <shell脚本字符串，用注释记录决策>\n\n'
        f"执行上下文：\n{json.dumps(ctx, indent=2, ensure_ascii=False)}"
    )

    result: Dict[str, str] = {}
    try:
        if llm_method == "openai":
            import openai
            client = openai.OpenAI(api_key=llm_api_key, base_url=llm_base_url)
            response = client.chat.completions.create(
                model=llm_model or "gpt-4o",
                messages=[{"role": "user", "content": prompt}],
                temperature=0.1,
                response_format={"type": "json_object"},
            )
            result = json.loads(response.choices[0].message.content)
        elif llm_method == "ollama":
            import requests
            url = llm_base_url or "http://localhost:11434"
            resp = requests.post(
                f"{url}/api/chat",
                json={"model": llm_model, "messages": [{"role": "user", "content": prompt}],
                       "stream": False, "format": "json"},
                timeout=120,
            )
            resp.raise_for_status()
            result = json.loads(resp.json().get("message", {}).get("content", "{}"))
    except Exception as exc:
        logging.warning("LLM audit report generation failed: %s", exc)

    # Write audit report
    audit_path = os.path.join(report_dir, "audit_report.md")
    with open(audit_path, "w", encoding="utf-8") as f:
        f.write(result.get("audit_report", f"# Audit Report\n\nLLM generation failed. Raw context:\n\n```json\n{json.dumps(ctx, indent=2)}\n```\n"))
    logging.info("Audit report: %s", audit_path)

    # Write reproducible script
    script_path = os.path.join(report_dir, "decision_log.sh")
    with open(script_path, "w", encoding="utf-8") as f:
        f.write(result.get("decision_log", f"#!/bin/bash\n# LLM生成失败。需要手动重建决策日志。\n"))
    os.chmod(script_path, 0o755)
    logging.info("Decision log: %s", script_path)


def mode_auto(
    adata: ad.AnnData,
    output: str,
    input_path: str = "",
    tissue: str = "",
    llm_method: str = "openai",
    llm_model: str = "",
    llm_api_key: str = "",
    llm_base_url: str = "",
    resolution: float = 0.8,
    max_iterations: int = 5,
    min_genes: int = 800,
    min_counts: int = 3000,
    max_pct_mt: float = 20.0,
    n_pcs: int = 50,
    n_neighbors: int = 50,
    n_top_genes: int = 3000,
    batch_method: str = "harmony",
    batch_key: str = "",
    auto_n_pcs: bool = False,
    plot_dir: str = "",
    skip_te: bool = False,
) -> None:
    """Fully autonomous: cluster → AI annotate → QC → filter → re-cluster.

    Iterates until no problematic clusters remain or *max_iterations* is
    reached.  Each iteration:
        1. Cluster (same logic as mode_cluster).
        2. AI-annotate each cluster (LLM + PubMed).
        3. Analyze cluster quality.
        4. If flagged clusters exist: filter cells (not clusters), repeat.
        5. If clean: save final h5ad + reports.

    Args:
        adata: Merged AnnData (from mode_merge).
        output: Path to write the final annotated h5ad.
        tissue: Tissue context for annotation prompts.
        llm_method: LLM backend (``"openai"`` or ``"ollama"``).
        llm_model: Model identifier.
        llm_api_key: API key for OpenAI.
        llm_base_url: Base URL for LLM API.
        resolution: Leiden clustering resolution (default 0.8).
        max_iterations: Max QC refinement iterations (default 3).
        min_genes: Min genes per cell for QC filtering.
        min_counts: Min UMI per cell for QC filtering.
        max_pct_mt: Max MT% for QC filtering.
        n_pcs: Number of PCA components.
        n_neighbors: Number of neighbours.
        n_top_genes: Number of HVGs.
        batch_method: Batch correction method.
        batch_key: Batch key column.
        auto_n_pcs: Auto-detect n_pcs.
        plot_dir: Directory for plots.
    """
    output_dir = os.path.dirname(output) or "."
    output_stem = Path(output).stem
    report_dir = os.path.join(output_dir, f"{output_stem}_reports")
    os.makedirs(report_dir, exist_ok=True)

    # ── Iteration context for LLM-generated audit report ──
    ctx: Dict[str, Any] = {
        "input_path": input_path or "(adata_loaded)",
        "output": output,
        "tissue": tissue,
        "params": {
            "llm_method": llm_method, "llm_model": llm_model,
            "resolution": resolution, "max_iterations": max_iterations,
            "min_genes": min_genes, "min_counts": min_counts,
            "max_pct_mt": max_pct_mt, "n_pcs": n_pcs,
            "n_neighbors": n_neighbors, "n_top_genes": n_top_genes,
            "batch_method": batch_method, "batch_key": batch_key,
            "auto_n_pcs": auto_n_pcs, "skip_te": skip_te,
        },
        "initial_cells": adata.n_obs,
        "initial_genes": adata.n_vars,
        "iterations": [],
    }
    iter_ctx: Dict[str, Any] = {}  # current iteration context

    # Save raw counts for re-clustering
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    raw_adata = adata.copy()

    iteration = 0
    while iteration < max_iterations:
        iteration += 1
        logging.info("=" * 60)
        logging.info("AUTO ITERATION %d / %d (%d cells)",
                     iteration, max_iterations, adata.n_obs)
        logging.info("=" * 60)

        # ── Step 1: Cluster ──
        logging.info("[Step 1] Clustering (resolution=%.2f)...", resolution)
        resolved_batch_key = batch_key or ("sample_id" if "sample_id" in adata.obs else "batch")

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata.copy()

        if skip_te:
            adata = _filter_te(adata)

        try:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top_genes, flavor="seurat", subset=False,
                batch_key=resolved_batch_key if resolved_batch_key in adata.obs else None,
            )
        except ValueError:
            logging.warning("HVG with batch_key failed (likely too few cells per batch after filtering). Retrying without batch_key.")
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top_genes, flavor="seurat", subset=False,
            )
        adata = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata, max_value=10)

        pca_comps = max(n_pcs, 100) if auto_n_pcs else n_pcs
        pca_comps = min(pca_comps, max(2, adata.n_obs - 1))
        sc.tl.pca(adata, n_comps=pca_comps)

        if auto_n_pcs:
            variance_ratio = adata.uns["pca"]["variance_ratio"]
            n_pcs_detected, _ = detect_n_pcs(variance_ratio)
            n_pcs = n_pcs_detected
            logging.info("Auto-detected n_pcs: %d", n_pcs)

        if batch_method == "harmony":
            ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, resolved_batch_key)
            Z = np.asarray(ho.Z_corr)
            if Z.ndim == 1:
                raise ValueError(f"harmonypy Z_corr is 1D shape={Z.shape}")
            if Z.shape[0] != adata.n_obs:
                Z = Z.T
            adata.obsm["X_pca_harmony"] = Z
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                            use_rep="X_pca_harmony")
        elif batch_method == "bbknn":
            sc.external.pp.bbknn(adata, batch_key=resolved_batch_key)
        else:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

        sc.tl.umap(adata, min_dist=0.1, spread=0.8)
        sc.tl.leiden(adata, resolution=resolution, key_added="leiden",
                     flavor="igraph", n_iterations=2, directed=False)
        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

        n_clusters = len(adata.obs["leiden"].unique())
        logging.info("Clustering done: %d clusters", n_clusters)
        iter_ctx = {
            "iteration": iteration,
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "resolution": resolution,
            "n_clusters": n_clusters,
            "batch_method": batch_method,
            "n_hvg": int(adata.var["highly_variable"].sum()),
            "n_pcs_used": n_pcs,
            "cluster_sizes": {str(k): int(v) for k, v in
                              adata.obs["leiden"].value_counts().items()},
        }

        # ── Step 2: AI annotate each cluster ──
        logging.info("[Step 2] AI annotation...")
        annotations: Dict[str, Dict] = {}

        # Compute UMAP centroids for distance calculation
        umap_centroids = {}
        if "X_umap" in adata.obsm:
            for cluster in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
                mask = adata.obs["leiden"] == cluster
                if int(mask.sum()) >= 50:  # Min cells for valid centroid
                    umap_centroids[cluster] = adata.obsm["X_umap"][mask.values].mean(axis=0)

        for cluster in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
            mask = adata.obs["leiden"] == cluster
            n_cells = int(mask.sum())
            obs = adata.obs[mask]
            mean_g = float(obs["n_genes_by_counts"].mean()) if "n_genes_by_counts" in obs.columns else 0.0
            mean_c = float(obs["total_counts"].mean()) if "total_counts" in obs.columns else 0.0
            pmt = float(obs["pct_counts_mt"].mean()) if "pct_counts_mt" in obs.columns else 0.0

            marker_df = sc.get.rank_genes_groups_df(adata, group=cluster)
            top_genes = marker_df.head(50)["names"].tolist()

            # Compute UMAP distances to other clusters
            umap_distances = {}
            if cluster in umap_centroids:
                for other_cluster, other_centroid in umap_centroids.items():
                    if other_cluster != cluster:
                        dist = float(np.linalg.norm(umap_centroids[cluster] - other_centroid))
                        umap_distances[other_cluster] = dist

            logging.info("  Cluster %s (%d cells): annotating...", cluster, n_cells)
            ann = _ai_annotate_cluster(
                cluster, top_genes, tissue, n_cells, mean_g, mean_c, pmt,
                llm_method=llm_method, llm_model=llm_model,
                llm_api_key=llm_api_key, llm_base_url=llm_base_url,
                other_annotations=annotations,
                umap_distances=umap_distances,
            )
            annotations[cluster] = ann
            logging.info("    -> %s (confidence: %s, %d refs, is_subcluster: %s, should_merge: %s)",
                         ann["cell_type"], ann["confidence"],
                         sum(len(v) for v in ann.get("references", {}).values()),
                         ann.get("is_subcluster", False),
                         ann.get("should_merge", False))

            time.sleep(0.5)  # rate limit between clusters

        # Apply annotations to adata
        cluster_to_ct = {cl: ann["cell_type"] for cl, ann in annotations.items()}
        adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_ct).astype("category")
        iter_ctx["annotations"] = {cl: {"cell_type": ann["cell_type"],
                                         "confidence": ann["confidence"],
                                         "key_markers": ann.get("key_markers", []),
                                         "reasoning": ann.get("reasoning", "")}
                                    for cl, ann in annotations.items()}
        iter_ctx["cell_types"] = sorted(set(cluster_to_ct.values()))

        # ── Step 3: Quality analysis ──
        logging.info("[Step 3] Quality analysis...")
        quality_reports = []
        for cluster in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
            ann = annotations.get(cluster, {})
            qr = _analyze_cluster_quality(
                adata, cluster, ann,
                min_genes=min_genes, min_counts=min_counts, max_pct_mt=max_pct_mt,
            )
            quality_reports.append(qr)
            if qr["should_filter"]:
                logging.info("  Cluster %s FLAGGED: %s (%s)",
                             cluster, qr["ai_cell_type"], ",".join(qr["flags"]))

        flagged = [r for r in quality_reports if r["should_filter"]]
        logging.info("Flagged clusters: %d / %d", len(flagged), n_clusters)
        iter_ctx["n_flagged"] = len(flagged)
        iter_ctx["flagged_clusters"] = [{"cluster": r["cluster"],
                                          "cell_type": r["ai_cell_type"],
                                          "flags": r["flags"]}
                                         for r in flagged]
        iter_ctx["quality_reports"] = [{"cluster": r["cluster"],
                                         "cell_type": r["ai_cell_type"],
                                         "should_filter": r["should_filter"],
                                         "flags": r["flags"]}
                                        for r in quality_reports]

                                        # ── Step 3.25: Check should_merge flags from LLM ──
        merge_candidates = []
        for cluster, ann in annotations.items():
            if ann.get("should_merge") and ann.get("parent_cluster"):
                parent = ann["parent_cluster"]
                if parent in annotations:
                    merge_candidates.append({
                        "cluster": cluster,
                        "parent": parent,
                        "cell_type": ann["cell_type"],
                        "parent_cell_type": annotations[parent]["cell_type"],
                        "is_subcluster": ann.get("is_subcluster", False),
                    })
                    logging.info("  Cluster %s (%s) -> should merge with Cluster %s (%s)",
                                 cluster, ann["cell_type"], parent, annotations[parent]["cell_type"])

        if merge_candidates:
            logging.info("LLM identified %d clusters that should be merged", len(merge_candidates))
            # Store merge info in annotations for report
            for mc in merge_candidates:
                annotations[mc["cluster"]]["merge_info"] = {
                    "parent_cluster": mc["parent"],
                    "parent_cell_type": mc["parent_cell_type"],
                }
        else:
            logging.info("No clusters identified for merging by LLM")

        # ── Step 3.5: Check cell type separation ──
        logging.info("[Step 3.5] Checking cell type separation...")
        misannotated: List[Dict] = []
        needs_reclustering, separation_diag = _check_cell_type_separation(
            adata,
            cell_type_key="cell_type",
            cluster_key="leiden",
            distance_threshold=5.0,
            annotations=annotations,
        )

        if needs_reclustering and separation_diag.get("suggested_resolution"):
            suggested_res = separation_diag["suggested_resolution"]
            if suggested_res < resolution:
                logging.info("Cell type separation detected. Adjusting resolution: %.2f -> %.2f",
                             resolution, suggested_res)
                resolution = suggested_res
                # If separation is severe and no quality issues, re-cluster immediately
                if not flagged:
                    logging.info("No quality issues. Re-clustering with lower resolution...")
                    continue  # Re-cluster with adjusted resolution
                else:
                    logging.info("Quality issues present. Will filter first, then re-cluster.")

        # Handle misannotated clusters (low marker overlap → annotation error)
        misannotated = separation_diag.get("misannotated", [])
        for ma in misannotated:
            flag_cluster = ma["flagged_cluster"]
            logging.info("Misannotation detected: cluster %s labeled '%s' but markers diverge (Jaccard=%.3f). Flagging for filter.",
                         flag_cluster, ma["cell_type"], ma["marker_jaccard"])
            # Add to quality_reports if not already flagged
            if not any(r["cluster"] == flag_cluster for r in quality_reports if r["should_filter"]):
                quality_reports.append({
                    "cluster": flag_cluster,
                    "ai_cell_type": ma["cell_type"],
                    "should_filter": True,
                    "flags": ["misannotated_low_marker_overlap"],
                })
                flagged = [r for r in quality_reports if r["should_filter"]]
        iter_ctx["misannotated"] = misannotated

        # ── Step 4: Filter or finish ──
        if not flagged:
            logging.info("All clusters clean. Saving final results.")
            iter_ctx["outcome"] = "clean"
            ctx["iterations"].append(iter_ctx)
            break

        # Filter cells within flagged clusters (even on last iteration)
        logging.info("[Step 4] Filtering cells in flagged clusters...")
        adata_filtered, n_removed = _filter_flagged_cells(
            adata, quality_reports,
            min_genes=min_genes, min_counts=min_counts, max_pct_mt=max_pct_mt,
        )

        if n_removed == 0:
            logging.info("No cells removed after QC. Stopping iterations.")
            iter_ctx["outcome"] = "no_cells_removed"
            ctx["iterations"].append(iter_ctx)
            break

        # Stop after filtering on last iteration (don't re-cluster)
        if iteration >= max_iterations:
            logging.warning("Max iterations (%d) reached after filtering. Saving current state.", max_iterations)
            iter_ctx["outcome"] = "max_iterations"
            iter_ctx["filter_n_removed"] = n_removed
            adata = adata_filtered  # Use filtered data for final output
            ctx["iterations"].append(iter_ctx)
            break

        # Prepare raw counts for re-clustering
        # Use the raw_adata (pre-processing) and subset to filtered cells
        filtered_indices = adata_filtered.obs.index
        raw_filtered = raw_adata[raw_adata.obs.index.isin(filtered_indices)].copy()

        # Preserve sample_id
        if "sample_id" in adata_filtered.obs.columns:
            raw_filtered.obs["sample_id"] = adata_filtered.obs.loc[
                raw_filtered.obs.index, "sample_id"
            ]

        adata = raw_filtered
        logging.info("Iteration %d complete. %d cells remaining.", iteration, adata.n_obs)
        iter_ctx["filter_n_removed"] = n_removed
        iter_ctx["filter_cells_remaining"] = adata.n_obs
        ctx["iterations"].append(iter_ctx)

    # ── Final: save reports and h5ad ──
    logging.info("Writing final reports...")
    _write_annotation_report(annotations, quality_reports,
                            os.path.join(report_dir, "annotation_report.tsv"))
    _write_references_report(annotations,
                            os.path.join(report_dir, "references.tsv"))

    # Plot if requested
    plotter = _make_plotter(plot_dir)
    if plotter:
        plotter.plot_cluster(adata, cluster_key="leiden", sample_key=resolved_batch_key)
        
        # Generate annotation plots (UMAP cell_type, dotplot, etc.)
        annotation_keys = []
        if "cell_type" in adata.obs.columns:
            annotation_keys.append("cell_type")
        if "llm_label" in adata.obs.columns:
            annotation_keys.append("llm_label")
        
        has_rank_genes = "rank_genes_groups" in adata.uns
        plotter.plot_annotate(
            adata,
            marker_file="",  # No marker file in auto mode
            annotate_group="leiden",
            annotation_keys=annotation_keys,
            has_rank_genes=has_rank_genes,
        )

    adata.write_h5ad(output)
    ctx["final_cells"] = adata.n_obs
    ctx["final_genes"] = adata.n_vars
    logging.info("Final annotated h5ad: %s", output)

    # ── LLM-generated audit report & reproducible script ──
    _generate_audit_report(ctx, report_dir, llm_method, llm_model,
                           llm_api_key, llm_base_url)

    logging.info("Reports: %s", report_dir)
    logging.info("AUTO MODE COMPLETE")


# ---------------------------------------------------------------------------
# Advanced analysis
# ---------------------------------------------------------------------------
def mode_advanced(
    adata: ad.AnnData,
    output: str,
    n_pcs: int = 50,
    n_neighbors: int = 50,
    trajectory: bool = False,
    velocity: bool = False,
    communication: bool = False,
    cnv: bool = False,
    gtf: str = "",
    cnv_reference: str = "",
    plot_dir: str = "",
) -> None:
    """Run advanced downstream analyses on a clustered AnnData.

    Available analyses (all optional, toggled by flags):
        - **Trajectory**: diffusion map + diffusion pseudotime (DPT).
        - **RNA velocity**: scVelo stochastic model (requires ``spliced``
          and ``unspliced`` layers).
        - **Cell–cell communication**: LIANA rank-aggregate ligand–receptor
          analysis.
        - **Copy-number variation**: inferCNVpy with optional GTF annotation
          and reference cell types.

    Args:
        adata: Clustered AnnData object.
        output: Path to write the annotated h5ad file.
        n_pcs: Number of PCA components for velocity computation (default 50).
        n_neighbors: Number of neighbours for velocity moments (default 50).
        trajectory: If True, compute diffusion map and pseudotime.
        velocity: If True, run scVelo RNA velocity analysis.
        communication: If True, run LIANA cell–cell communication.
        cnv: If True, run inferCNVpy copy-number inference.
        gtf: Path to a GTF file for genomic coordinate annotation
            (used with CNV). Empty string skips GTF loading.
        cnv_reference: Comma-separated cell type names to use as normal
            reference for CNV. Empty string uses ``leiden`` clusters.
        plot_dir: Directory for advanced analysis plots. Empty string
            disables plotting.
    """
    if trajectory:
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)

    if velocity and {"spliced", "unspliced"}.issubset(adata.layers):
        import scvelo as scv
        scv.pp.moments(adata, n_pcs=min(n_pcs, 30), n_neighbors=n_neighbors)
        scv.tl.velocity(adata)
        scv.tl.velocity_graph(adata)

    if communication:
        import liana as li
        li.mt.rank_aggregate(adata, groupby="leiden", use_raw=False, verbose=False)

    if cnv:
        import infercnvpy as cnv_mod
        if gtf:
            cnv_mod.io.genomic_position_from_gtf(gtf, adata)
        ref_cats = (
            [c.strip() for c in cnv_reference.split(",") if c.strip()]
            if cnv_reference else None
        )
        if ref_cats and "cell_type" in adata.obs:
            cnv_mod.tl.infercnv(adata, reference_key="cell_type",
                            reference_cat=ref_cats, window_size=250)
        else:
            cnv_mod.tl.infercnv(adata, reference_key="leiden")

    plotter = _make_plotter(plot_dir)
    if plotter:
        # Detect annotation key from adata for trajectory/CNV coloring
        _anno_key = None
        for k in ("cell_type", "celltypist_label", "llm_label"):
            if k in adata.obs.columns:
                _anno_key = k
                break
        plotter.plot_advanced(adata, trajectory=trajectory, cnv=cnv,
                              annotation_key=_anno_key, cluster_key="leiden")

    adata.write_h5ad(output)


# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------
def mode_de(
    adata: ad.AnnData,
    output: str,
    deg: str = "",
    plot_dir: str = "",
) -> None:
    """Differential expression analysis between conditions or clusters.

    Uses ``rank_genes_groups`` (Wilcoxon) to identify DEGs.  Grouping is
    by ``condition`` if present in ``adata.obs``, otherwise by ``leiden``.

    Args:
        adata: Clustered AnnData object.
        output: Path to write the h5ad file with DE results.
        deg: Path to write a TSV of ranked DEGs. Empty string skips export.
        plot_dir: Directory for DE plots. Empty string disables plotting.
    """
    group = "condition" if "condition" in adata.obs else "leiden"
    sc.tl.rank_genes_groups(adata, group, method="wilcoxon")

    if deg:
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(
            deg, sep="\t", index=False
        )

    plotter = _make_plotter(plot_dir)
    if plotter:
        plotter.plot_de(adata, group)

    adata.write_h5ad(output)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    """Parse CLI arguments and dispatch to the requested mode function."""
    setup_logging()
    parser = argparse.ArgumentParser(description="Scanpy scRNA-seq pipeline")
    parser.add_argument("--mode", required=True,
                        choices=["qc", "merge", "cluster", "annotate", "auto", "advanced", "de"])
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--plot-dir", default="", help="Directory to save plots (optional)")

    # QC params
    parser.add_argument("--metrics", default="", help="Path to write per-cell QC metrics TSV")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=6000)
    parser.add_argument("--max-pct-mt", type=float, default=20)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--scrublet", action="store_true", help="Run Scrublet doublet detection")
    parser.add_argument("--doublet-rate", type=float, default=0.06, help="Expected doublet rate for Scrublet")
    parser.add_argument("--use-mad", action="store_true", help="Use MAD-based outlier detection (more permissive)")

    # Cluster params
    parser.add_argument("--n-pcs", type=int, default=50, help="Number of principal components")
    parser.add_argument("--n-neighbors", type=int, default=50, help="Number of k-NN neighbours")
    parser.add_argument("--resolution", type=float, default=0.8, help="Leiden clustering resolution")
    parser.add_argument("--markers", default="", help="Path to write ranked marker gene TSV")
    parser.add_argument("--auto-n-pcs", action="store_true", help="Auto-detect optimal n_pcs from PCA variance ratio")
    parser.add_argument("--skip-te", action="store_true", help="Exclude TE genes before HVG selection and clustering")

    # Batch params (integrated into cluster mode)
    parser.add_argument("--batch-method", default="harmony", choices=["harmony", "bbknn", ""],
                        help="Batch correction method (default: harmony, empty to skip)")
    parser.add_argument("--batch-key", default="", help="Column in obs identifying batches")

    # Annotate params
    parser.add_argument("--marker-file", default="", help="TSV with cell_type and markers columns")
    parser.add_argument("--celltypist-model", default="", help="CellTypist model name")
    parser.add_argument("--llm-method", default="", choices=["", "openai", "ollama", "file"],
                        help="LLM backend for annotation")
    parser.add_argument("--llm-model", default="", help="LLM model identifier")
    parser.add_argument("--llm-api-key", default="", help="API key for OpenAI backend")
    parser.add_argument("--llm-base-url", default="", help="Base URL for LLM API")
    parser.add_argument("--llm-top-genes", type=int, default=30, help="Top DEGs per cluster for LLM prompt")
    parser.add_argument("--annotate-group", default="", help="Obs column for cluster grouping")
    parser.add_argument("--tissue", default="", help="Tissue name for LLM prompt context")

    # Auto mode params
    parser.add_argument("--max-iterations", type=int, default=3, help="Max QC refinement iterations for auto mode")
    parser.add_argument("--min-counts", type=int, default=3000, help="Min UMI per cell for auto mode QC")

    # Advanced params
    parser.add_argument("--trajectory", action="store_true", help="Compute diffusion map + DPT")
    parser.add_argument("--velocity", action="store_true", help="Run scVelo RNA velocity")
    parser.add_argument("--communication", action="store_true", help="Run LIANA cell-cell communication")
    parser.add_argument("--cnv", action="store_true", help="Run inferCNVpy CNV inference")
    parser.add_argument("--gtf", default="", help="GTF file for CNV genomic annotation")
    parser.add_argument("--cnv-reference", default="", help="Comma-separated reference cell types for CNV")

    # DE params
    parser.add_argument("--deg", default="", help="Path to write DEG results TSV")

    # Gene type annotation params (for merge mode)
    parser.add_argument("--te-bed", default="", help="Path to TE BED file for gene_type annotation")
    parser.add_argument("--gene-tsv", default="", help="Path to gene annotation TSV for gene_type annotation")

    args = parser.parse_args()
    adata = read_input(args.input[0], args.input[1:] if len(args.input) > 1 else None)

    if args.mode == "qc":
        mode_qc(
            adata,
            output=args.output,
            min_genes=args.min_genes,
            max_genes=args.max_genes,
            max_pct_mt=args.max_pct_mt,
            use_mad=args.use_mad,
            scrublet=args.scrublet,
            doublet_rate=args.doublet_rate,
            plot_dir=args.plot_dir,
            metrics=args.metrics,
        )
    elif args.mode == "merge":
        mode_merge(
            adata,
            input_paths=args.input,
            output=args.output,
            plot_dir=args.plot_dir,
            te_bed=args.te_bed,
            gene_tsv=args.gene_tsv,
        )
    elif args.mode == "cluster":
        mode_cluster(
            adata,
            output=args.output,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            resolution=args.resolution,
            n_top_genes=args.n_top_genes,
            batch_method=args.batch_method,
            batch_key=args.batch_key,
            markers=args.markers,
            plot_dir=args.plot_dir,
            auto_n_pcs=args.auto_n_pcs,
            skip_te=args.skip_te,
        )
    elif args.mode == "annotate":
        mode_annotate(
            adata,
            output=args.output,
            marker_file=args.marker_file,
            celltypist_model=args.celltypist_model,
            llm_method=args.llm_method,
            llm_model=args.llm_model,
            llm_api_key=args.llm_api_key,
            llm_base_url=args.llm_base_url,
            llm_top_genes=args.llm_top_genes,
            annotate_group=args.annotate_group,
            tissue=args.tissue,
            plot_dir=args.plot_dir,
        )
    elif args.mode == "auto":
        mode_auto(
            adata,
            output=args.output,
            input_path=args.input[0],
            tissue=args.tissue,
            llm_method=args.llm_method,
            llm_model=args.llm_model,
            llm_api_key=args.llm_api_key,
            llm_base_url=args.llm_base_url,
            resolution=args.resolution,
            max_iterations=args.max_iterations,
            min_genes=args.min_genes,
            min_counts=args.min_counts,
            max_pct_mt=args.max_pct_mt,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            n_top_genes=args.n_top_genes,
            batch_method=args.batch_method,
            batch_key=args.batch_key,
            auto_n_pcs=args.auto_n_pcs,
            plot_dir=args.plot_dir,
            skip_te=args.skip_te,
        )
    elif args.mode == "advanced":
        mode_advanced(
            adata,
            output=args.output,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            trajectory=args.trajectory,
            velocity=args.velocity,
            communication=args.communication,
            cnv=args.cnv,
            gtf=args.gtf,
            cnv_reference=args.cnv_reference,
            plot_dir=args.plot_dir,
        )
    elif args.mode == "de":
        mode_de(
            adata,
            output=args.output,
            deg=args.deg,
            plot_dir=args.plot_dir,
        )


if __name__ == "__main__":
    main()
