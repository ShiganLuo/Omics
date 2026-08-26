"""Visualization module for Scanpy scRNA-seq pipeline.

Provides a ScanpyPlotter class that generates publication-quality figures
for each pipeline stage: QC, merge, cluster, batch, annotate, advanced, de.

Usage:
    from plot import ScanpyPlotter
    plotter = ScanpyPlotter(plot_dir="/path/to/plots")
    plotter.plot_qc(adata, adata_before)
"""

import csv
import logging
import os
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default style
# ---------------------------------------------------------------------------
sc.set_figure_params(dpi=80, facecolor="white", frameon=False)


class ScanpyPlotter:
    """Generates pipeline stage plots for scanpy scRNA-seq analysis.

    Args:
        plot_dir: Output directory for PNG figures.
        dpi: Resolution for saved figures.
    """

    def __init__(self, plot_dir: str, dpi: int = 300) -> None:
        self.plot_dir = plot_dir
        self.dpi = dpi
        os.makedirs(plot_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _save(self, filename: str) -> None:
        """Save current matplotlib figure and release memory."""
        path = os.path.join(self.plot_dir, filename)
        plt.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close("all")

    @staticmethod
    def _detect_anno_key(adata: AnnData) -> Optional[str]:
        """Return the first available cell-type annotation column."""
        for key in ("cell_type", "celltypist_label", "llm_label", "major_celltype"):
            if key in adata.obs.columns:
                return key
        return None

    @staticmethod
    def _detect_sample_key(adata: AnnData) -> Optional[str]:
        """Return the first available sample identifier column."""
        for key in ("sample_id", "sample", "batch"):
            if key in adata.obs.columns:
                return key
        return None

    @staticmethod
    def _umap(adata: AnnData, color: str, *,
              ax: plt.Axes, title: str = "", legend_loc: str = "right margin",
              legend_fontsize: int = 8, legend_fontoutline: int = 2) -> None:
        """Thin wrapper around sc.pl.umap with standard styling."""
        sc.pl.umap(
            adata, color=color, frameon=False, show=False, ax=ax,
            title=title or color,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            legend_fontoutline=legend_fontoutline,
        )

    # ------------------------------------------------------------------
    # QC
    # ------------------------------------------------------------------
    def plot_qc(self, adata: AnnData, adata_before: AnnData) -> None:
        """QC diagnostics: counts distribution, mt%, scatter before/after, HVG.

        Args:
            adata: Filtered & normalised AnnData.
            adata_pre-filter AnnData (with QC metrics computed).
        """
        # 1. Total counts distribution
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.histplot(adata_before.obs["total_counts"], bins=100, kde=False, ax=ax)
        ax.set_title("Total Counts Distribution")
        self._save("qc_total_counts_distribution.png")

        # 2. Violin: pct_counts_mt
        fig, ax = plt.subplots(figsize=(4, 6))
        sc.pl.violin(adata_before, "pct_counts_mt", show=False, ax=ax)
        ax.set_title("% Mitochondrial Counts")
        self._save("qc_pct_counts_mt_violin.png")

        # 3. Scatter: before filtering
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.scatter(
            adata_before, "total_counts", "n_genes_by_counts",
            color="pct_counts_mt", show=False, ax=ax,
        )
        ax.set_title("Before Filtering")
        self._save("qc_scatter_before_filter.png")

        # 4. Scatter: after filtering
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.scatter(
            adata, "total_counts", "n_genes_by_counts",
            color="pct_counts_mt", show=False, ax=ax,
        )
        ax.set_title("After Filtering")
        self._save("qc_scatter_after_filter.png")

        # 5. Highly variable genes (does not accept ax parameter)
        sc.pl.highly_variable_genes(adata, show=False)
        self._save("qc_highly_variable_genes.png")

        logger.info("QC plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------
    def plot_merge(self, merged: AnnData) -> None:
        """Post-merge overview: per-sample cell counts bar chart.

        Args:
            merged: Concatenated AnnData with ``sample_id`` in .obs.
        """
        sample_key = "sample_id" if "sample_id" in merged.obs.columns else "sample"
        if sample_key not in merged.obs.columns:
            logger.warning("No sample column found, skipping merge plots")
            return

        counts = merged.obs[sample_key].value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(max(6, len(counts) * 0.8), 5))
        counts.plot(kind="bar", ax=ax, color="steelblue")
        ax.set_title("Cell Counts per Sample")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Number of Cells")
        plt.xticks(rotation=45, ha="right")
        self._save("merge_cell_counts_per_sample.png")

        logger.info("Merge plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Cluster
    # ------------------------------------------------------------------
    def plot_cluster(self, adata: AnnData) -> None:
        """Clustering results: PCA variance, UMAP by leiden & sample.

        Args:
            adata: Clustered AnnData with ``X_pca``, ``X_umap``, ``leiden``.
        """
        # 1. PCA variance ratio (does not accept ax parameter)
        sc.pl.pca_variance_ratio(adata, show=False)
        self._save("cluster_pca_variance_ratio.png")

        # 2. UMAP — leiden
        fig, ax = plt.subplots(figsize=(8, 6))
        self._umap(adata, "leiden", ax=ax, title="Leiden Clusters",
                   legend_loc="on data", legend_fontsize=10)
        self._save("cluster_umap_leiden.png")

        # 3. UMAP — sample
        sample_key = self._detect_sample_key(adata)
        if sample_key:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._umap(adata, sample_key, ax=ax, title="Samples")
            self._save("cluster_umap_sample.png")

        logger.info("Cluster plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Batch correction
    # ------------------------------------------------------------------
    def plot_batch(self, adata: AnnData, method: str) -> None:
        """Batch correction results: UMAP by sample & leiden, side-by-side.

        Args:
            adata: Batch-corrected AnnData.
            method: Correction method name (``bbknn`` / ``harmony``).
        """
        sample_key = self._detect_sample_key(adata)
        title_suffix = method.upper()

        # 1. UMAP — sample (batch effect)
        if sample_key:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._umap(adata, sample_key, ax=ax, title=f"After {title_suffix} — Samples")
            self._save("batch_umap_sample.png")

        # 2. UMAP — leiden
        fig, ax = plt.subplots(figsize=(8, 6))
        self._umap(adata, "leiden", ax=ax, title=f"After {title_suffix} — Leiden",
                   legend_loc="on data", legend_fontsize=10)
        self._save("batch_umap_leiden.png")

        # 3. Side-by-side comparison
        if sample_key:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            plt.subplots_adjust(wspace=0.4)
            self._umap(adata, sample_key, ax=ax1, title="Sample")
            self._umap(adata, "leiden", ax=ax2, title="Leiden",
                       legend_loc="on data", legend_fontsize=10)
            self._save("batch_umap_side_by_side.png")

        logger.info("Batch correction plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Annotation
    # ------------------------------------------------------------------
    def plot_annotate(self, adata: AnnData,
                      marker_file: str = "",
                      annotate_group: str = "") -> None:
        """Cell type annotation: marker dotplot, DEG overview, UMAP labels.

        Args:
            adata: Annotated AnnData.
            marker_file: Path to marker TSV (if used).
            annotate_group: Obs column used for grouping (default ``leiden``).
        """
        group = annotate_group or ("leiden" if "leiden" in adata.obs.columns else None)

        # 1. Marker dotplot (from marker file)
        if marker_file and group:
            marker_dict = self._load_marker_dict(marker_file, adata)
            if marker_dict:
                sc.pl.dotplot(
                    adata, var_names=marker_dict, groupby=group,
                    standard_scale="var", show=False,
                )
                self._save("annotate_marker_dotplot.png")

        # 2. Rank genes groups overview
        if group and "rank_genes_groups" in adata.uns:
            sc.pl.rank_genes_groups(adata, n_genes=5, show=False)
            self._save("annotate_rank_genes_groups.png")

            sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, show=False)
            self._save("annotate_deg_dotplot.png")

        # 3. UMAP — each annotation column
        for col in ("cell_type", "celltypist_label", "llm_label"):
            if col in adata.obs.columns:
                fig, ax = plt.subplots(figsize=(8, 6))
                self._umap(adata, col, ax=ax, legend_loc="on data",
                           legend_fontsize=8)
                self._save(f"annotate_umap_{col}.png")

        # 4. CellTypist confidence
        if "celltypist_score" in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._umap(adata, "celltypist_score", ax=ax,
                       title="CellTypist Confidence Score")
            self._save("annotate_celltypist_confidence.png")

        logger.info("Annotation plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Advanced (trajectory / CNV)
    # ------------------------------------------------------------------
    def plot_advanced(self, adata: AnnData, *,
                      trajectory: bool = False,
                      cnv: bool = False) -> None:
        """Advanced analysis: diffusion map / pseudotime, CNV heatmaps.

        Args:
            adata: AnnData with advanced results.
            trajectory: Whether trajectory (diffmap/dpt) was run.
            cnv: Whether CNV (infercnvpy) was run.
        """
        anno_key = self._detect_anno_key(adata)
        leiden_key = "leiden" if "leiden" in adata.obs.columns else None

        # --- Trajectory ---
        if trajectory and "X_diffmap" in adata.obsm:
            color = anno_key or leiden_key or "dpt_pseudotime"
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.diffmap(adata, color=color, components=["2, 3"],
                          show=False, ax=ax, title="Diffusion Map")
            self._save("advanced_diffmap.png")

            if "dpt_pseudotime" in adata.obs.columns:
                fig, ax = plt.subplots(figsize=(8, 6))
                sc.pl.diffmap(adata, color="dpt_pseudotime",
                              components=["2, 3"], show=False, ax=ax,
                              title="Pseudotime")
                self._save("advanced_pseudotime.png")

        # --- CNV ---
        if cnv:
            self._plot_cnv(adata, anno_key, leiden_key)

        logger.info("Advanced analysis plots saved to %s", self.plot_dir)

    def _plot_cnv(self, adata: AnnData,
                  anno_key: Optional[str],
                  leiden_key: Optional[str]) -> None:
        """Generate CNV-specific plots (requires infercnvpy)."""
        try:
            import infercnvpy as cnv
        except ImportError:
            logger.warning("infercnvpy not installed, skipping CNV plots")
            return

        groupby = anno_key or leiden_key

        # 1. Chromosome heatmap by annotation
        if groupby:
            cnv.pl.chromosome_heatmap(adata, groupby=groupby, show=False)
            self._save("advanced_cnv_heatmap.png")

        # 2. CNV clustering + UMAP panels
        cnv.tl.pca(adata)
        cnv.pp.neighbors(adata)
        cnv.tl.leiden(adata)
        cnv.tl.umap(adata)
        cnv.tl.cnv_score(adata)

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes[1, 1].axis("off")
        cnv.pl.umap(adata, color="cnv_leiden", legend_loc="on data",
                    legend_fontoutline=2, ax=axes[0, 0], show=False)
        axes[0, 0].set_title("CNV Leiden")
        cnv.pl.umap(adata, color="cnv_score", ax=axes[0, 1], show=False)
        axes[0, 1].set_title("CNV Score")
        if anno_key:
            cnv.pl.umap(adata, color=anno_key, ax=axes[1, 0], show=False)
            axes[1, 0].set_title(anno_key)
        fig.suptitle("CNV Analysis", fontsize=14, y=0.98)
        self._save("advanced_cnv_umap.png")

        # 3. CNV heatmap by cnv_leiden
        cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden",
                                  dendrogram=True, show=False)
        self._save("advanced_cnv_heatmap_leiden.png")

    # ------------------------------------------------------------------
    # Differential expression
    # ------------------------------------------------------------------
    def plot_de(self, adata: AnnData, group: str) -> None:
        """Differential expression: rank_genes overview, volcano plot.

        Args:
            adata: AnnData with ``rank_genes_groups`` in .uns.
            group: Grouping key (e.g. ``leiden`` or ``condition``).
        """
        # 1. Rank genes groups overview
        sc.pl.rank_genes_groups(adata, n_genes=10, show=False)
        self._save("de_rank_genes_groups.png")

        # 2. Rank genes groups dotplot
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
        self._save("de_rank_genes_dotplot.png")

        # 3. Volcano plot (first group)
        try:
            first_group = adata.obs[group].cat.categories[0]
            self._volcano(adata, group=first_group)
        except Exception as exc:
            logger.warning("Volcano plot failed: %s", exc)

        logger.info("DE plots saved to %s", self.plot_dir)

    def _volcano(self, adata: AnnData, group: str) -> None:
        """Draw a volcano plot for one group vs rest."""
        df = sc.get.rank_genes_groups_df(adata, group=group)
        pval_col = next(
            (c for c in ("pvals_adj", "pval_adj") if c in df.columns), None
        )
        if pval_col is None or "logfoldchanges" not in df.columns:
            return

        df = df.copy()
        df["-log10pval"] = -np.log10(df[pval_col].clip(lower=1e-300))
        sig = (df[pval_col] < 0.05) & (df["logfoldchanges"].abs() > 1)
        df["significant"] = sig

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(df.loc[~sig, "logfoldchanges"], df.loc[~sig, "-log10pval"],
                   c="grey", s=5, alpha=0.5, label="NS")
        ax.scatter(df.loc[sig, "logfoldchanges"], df.loc[sig, "-log10pval"],
                   c="red", s=10, alpha=0.7, label="Significant")

        # Label top 10 most significant genes
        top = df.loc[sig].sort_values(pval_col).head(10)
        for _, row in top.iterrows():
            ax.annotate(
                str(row["names"]),
                (float(row["logfoldchanges"]), float(row["-log10pval"])),
                fontsize=7, alpha=0.8,
            )

        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10(adjusted p-value)")
        ax.set_title(f"Volcano Plot — {group}")
        ax.legend()
        ax.axhline(-np.log10(0.05), ls="--", color="grey", alpha=0.5)
        ax.axvline(-1, ls="--", color="grey", alpha=0.5)
        ax.axvline(1, ls="--", color="grey", alpha=0.5)
        self._save(f"de_volcano_{group}.png")

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------
    @staticmethod
    def _load_marker_dict(marker_file: str, adata: AnnData) -> Dict[str, List[str]]:
        """Parse marker TSV into {cell_type: [gene, ...]} (genes present in adata)."""
        markers: Dict[str, List[str]] = {}
        with open(marker_file, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ct = row["cell_type"]
                genes = [
                    g.strip() for g in row.get("markers", "").split(",")
                    if g.strip() and g.strip() in adata.var.index
                ]
                if genes:
                    markers[ct] = genes
        return markers
