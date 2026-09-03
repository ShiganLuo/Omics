"""Visualization module for Scanpy scRNA-seq pipeline.

Provides a ScanpyPlotter class that generates publication-quality figures
for each pipeline stage: QC, merge, cluster, batch, annotate, advanced, de.

Usage:
    from plot import ScanpyPlotter
    plotter = ScanpyPlotter(plot_dir="/path/to/plots")
    plotter.plot_qc(adata, adata_before, ...)
"""

import csv
import logging
import os
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
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
    def plot_qc(
        self,
        adata: AnnData,
        adata_before: AnnData,
        counts_col: str = "total_counts",
        genes_col: str = "n_genes_by_counts",
        mt_col: str = "pct_counts_mt",
    ) -> None:
        """QC diagnostics: counts distribution, mt%, scatter before/after.

        Generates:
            - qc_total_counts_distribution.png: histogram of total counts.
            - qc_pct_counts_mt_violin.png: violin plot of mitochondrial %.
            - qc_scatter_filter_comparison.png: counts vs genes before/after QC.

        Args:
            adata: Filtered AnnData (after QC filtering).
            adata_before: Pre-filter AnnData (with QC metrics computed).
            counts_col: Column in ``adata.obs`` for total UMI counts
                (default ``"total_counts"``).
            genes_col: Column in ``adata.obs`` for gene counts per cell
                (default ``"n_genes_by_counts"``).
            mt_col: Column in ``adata.obs`` for mitochondrial fraction
                (default ``"pct_counts_mt"``).
        """
        # 1. Total counts distribution
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.histplot(adata_before.obs[counts_col], bins=100, kde=False, ax=ax)
        ax.set_title("Total Counts Distribution")
        self._save("qc_total_counts_distribution.png")

        # 2. Violin: mt%
        fig, ax = plt.subplots(figsize=(4, 6))
        sc.pl.violin(adata_before, mt_col, show=False, ax=ax)
        ax.set_title(f"% Mitochondrial Counts ({mt_col})")
        self._save("qc_pct_counts_mt_violin.png")

        # Compute unified axis limits for before/after comparison
        all_counts = np.concatenate([adata_before.obs[counts_col].values,
                                     adata.obs[counts_col].values])  # type: ignore[arg-type]
        all_genes = np.concatenate([adata_before.obs[genes_col].values,
                                    adata.obs[genes_col].values])  # type: ignore[arg-type]
        all_mt = np.concatenate([adata_before.obs[mt_col].values,
                                 adata.obs[mt_col].values])  # type: ignore[arg-type]
        x_min, x_max = float(all_counts.min()), float(all_counts.max())
        y_min, y_max = float(all_genes.min()), float(all_genes.max())
        c_min, c_max = float(all_mt.min()), float(all_mt.max())
        x_margin = (x_max - x_min) * 0.05
        y_margin = (y_max - y_min) * 0.05
        x_lim = (x_min - x_margin, x_max + x_margin)
        y_lim = (y_min - y_margin, y_max + y_margin)

        # 3+4. Scatter: before vs after side by side
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        plt.subplots_adjust(wspace=0.35)

        for ax, src, title in [(ax1, adata_before, "Before Filtering"),
                                (ax2, adata, "After Filtering")]:
            sc.pl.scatter(src, counts_col, genes_col,
                          color=mt_col, show=False, ax=ax)
            ax.set_title(title)
            ax.set_xlim(x_lim)
            ax.set_ylim(y_lim)
            if ax.collections:
                ax.collections[0].set_clim(c_min, c_max)
                old_cbar = ax.collections[0].colorbar
                if old_cbar is not None:
                    old_cbar.remove()

        # Single colorbar on the right edge
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(ax2.collections[0], cax=cax)
        cbar.set_label(mt_col)

        self._save("qc_scatter_filter_comparison.png")

        logger.info("QC plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------
    def plot_merge(self, adata: AnnData, sample_key: str = "sample_id") -> None:
        """Post-merge overview: per-sample cell counts bar chart.

        Args:
            adata: Concatenated AnnData.
            sample_key: Column in ``adata.obs`` identifying samples
                (default ``"sample_id"``).
        """
        if sample_key not in adata.obs.columns:
            logger.warning("Column '%s' not found in obs, skipping merge plots", sample_key)
            return

        counts = adata.obs[sample_key].value_counts().sort_index()
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
    def plot_cluster(
        self,
        adata: AnnData,
        cluster_key: str = "leiden",
        sample_key: str = "sample_id",
    ) -> None:
        """Clustering results: PCA variance, HVG, UMAP by cluster & sample.

        Generates:
            - cluster_pca_variance_ratio.png: variance explained per PC.
            - cluster_highly_variable_genes.png: top HVGs highlighted.
            - cluster_umap_<cluster_key>.png: UMAP coloured by clusters.
            - cluster_umap_<sample_key>.png: UMAP coloured by sample.

        Args:
            adata: Clustered AnnData with ``X_pca``, ``X_umap``, and
                cluster labels in ``adata.obs[cluster_key]``.
            cluster_key: Column in ``adata.obs`` for cluster labels
                (default ``"leiden"``).
            sample_key: Column in ``adata.obs`` for sample identifiers
                (default ``"sample"``).
        """
        sc.pl.pca_variance_ratio(adata, show=False)
        self._save("cluster_pca_variance_ratio.png")

        if "highly_variable" in adata.var.columns:
            try:
                sc.pl.highly_variable_genes(adata, show=False)
                self._save("cluster_highly_variable_genes.png")
            except KeyError:
                logger.warning("HVG metadata incomplete, skipping HVG plot")
        
        fig, ax = plt.subplots(figsize=(8, 6))
        self._umap(adata, cluster_key, ax=ax,
                   title=f"{cluster_key} Clusters",
                   legend_loc="on data", legend_fontsize=10)
        self._save(f"cluster_umap_{cluster_key}.png")

        if sample_key in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._umap(adata, sample_key, ax=ax, title="Samples")
            self._save(f"cluster_umap_{sample_key}.png")

        logger.info("Cluster plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Batch correction
    # ------------------------------------------------------------------
    def plot_batch(
        self,
        adata: AnnData,
        method: str,
        cluster_key: str = "leiden",
        sample_key: str = "sample",
    ) -> None:
        """Batch correction results: UMAP by sample & cluster, side-by-side.

        Generates:
            - batch_umap_<sample_key>.png: UMAP coloured by sample.
            - batch_umap_<cluster_key>.png: UMAP coloured by clusters.
            - batch_umap_side_by_side.png: sample + cluster side by side.

        Args:
            adata: Batch-corrected AnnData with ``X_umap``.
            method: Correction method name (``"bbknn"`` / ``"harmony"``).
            cluster_key: Column in ``adata.obs`` for cluster labels
                (default ``"leiden"``).
            sample_key: Column in ``adata.obs`` for sample identifiers
                (default ``"sample"``).
        """
        title_suffix = method.upper()
        has_sample = sample_key in adata.obs.columns

        # 1. UMAP — sample (batch effect)
        if has_sample:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._umap(adata, sample_key, ax=ax,
                       title=f"After {title_suffix} — Samples")
            self._save(f"batch_umap_{sample_key}.png")

        # 2. UMAP — cluster
        fig, ax = plt.subplots(figsize=(8, 6))
        self._umap(adata, cluster_key, ax=ax,
                   title=f"After {title_suffix} — {cluster_key}",
                   legend_loc="on data", legend_fontsize=10)
        self._save(f"batch_umap_{cluster_key}.png")

        # 3. Side-by-side comparison
        if has_sample:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            plt.subplots_adjust(wspace=0.4)
            self._umap(adata, sample_key, ax=ax1, title="Sample")
            self._umap(adata, cluster_key, ax=ax2, title=cluster_key,
                       legend_loc="on data", legend_fontsize=10)
            self._save("batch_umap_side_by_side.png")

        logger.info("Batch correction plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Annotation
    # ------------------------------------------------------------------
    def plot_annotate(
        self,
        adata: AnnData,
        marker_file: str = "",
        annotate_group: str = "leiden",
        annotation_keys: Optional[List[str]] = None,
        score_col: str = "celltypist_score",
        has_rank_genes: bool = False,
    ) -> None:
        """Cell type annotation: marker dotplot, DEG overview, UMAP labels.

        Generates (when data available):
            - annotate_marker_dotplot.png: dotplot of marker gene expression.
            - annotate_rank_genes_groups.png: top DEGs per cluster.
            - annotate_deg_dotplot.png: DEG dotplot per cluster.
            - annotate_umap_<col>.png: UMAP for each annotation column.
            - annotate_<score_col>.png: confidence score UMAP.

        Args:
            adata: Annotated AnnData with clustering and optional annotation.
            marker_file: Path to marker TSV (empty string skips marker dotplot).
            annotate_group: Column in ``adata.obs`` for grouping
                (default ``"leiden"``).
            annotation_keys: List of ``adata.obs`` columns to plot as UMAP
                colourings (e.g. ``["cell_type", "celltypist_label"]``).
                Empty/None skips annotation UMAPs.
            score_col: Column in ``adata.obs`` for annotation confidence
                scores (default ``"celltypist_score"``).
            has_rank_genes: Whether ``rank_genes_groups`` has been computed
                and stored in ``adata.uns`` (default False).
        """
        if annotation_keys is None:
            annotation_keys = []

        # 1. Marker dotplot (from marker file)
        if marker_file and annotate_group:
            marker_dict = self._load_marker_dict(marker_file, adata)
            if marker_dict:
                sc.pl.dotplot(
                    adata, var_names=marker_dict, groupby=annotate_group,
                    standard_scale="var", show=False,
                )
                self._save("annotate_marker_dotplot.png")

        # 2. Rank genes groups overview
        if has_rank_genes:
            sc.pl.rank_genes_groups(adata, n_genes=5, show=False)
            self._save("annotate_rank_genes_groups.png")

            sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, show=False)
            self._save("annotate_deg_dotplot.png")

        # 3. UMAP — each annotation column
        for col in annotation_keys:
            if col in adata.obs.columns:
                fig, ax = plt.subplots(figsize=(8, 6))
                self._umap(adata, col, ax=ax, legend_loc="on data",
                           legend_fontsize=8)
                self._save(f"annotate_umap_{col}.png")

        # 4. Confidence score
        if score_col in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._umap(adata, score_col, ax=ax,
                       title=f"Confidence Score ({score_col})")
            self._save(f"annotate_{score_col}.png")

        logger.info("Annotation plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Advanced (trajectory / CNV)
    # ------------------------------------------------------------------
    def plot_advanced(
        self,
        adata: AnnData,
        trajectory: bool = False,
        cnv: bool = False,
        annotation_key: Optional[str] = None,
        cluster_key: str = "leiden",
        pseudotime_col: str = "dpt_pseudotime",
    ) -> None:
        """Advanced analysis: diffusion map / pseudotime, CNV heatmaps.

        Generates (when enabled):
            - advanced_diffmap.png: diffusion map coloured by annotation.
            - advanced_pseudotime.png: diffusion map coloured by pseudotime.
            - advanced_cnv_heatmap.png: CNV chromosome heatmap by annotation.
            - advanced_cnv_umap.png: CNV leiden/score/annotation UMAP panels.
            - advanced_cnv_heatmap_leiden.png: CNV heatmap by CNV-leiden.

        Args:
            adata: AnnData with advanced analysis results.
            trajectory: Whether trajectory (diffmap/dpt) was computed.
            cnv: Whether CNV (infercnvpy) was computed.
            annotation_key: Column in ``adata.obs`` for cell type annotation
                used as colour in trajectory/CNV plots. None falls back to
                *cluster_key*.
            cluster_key: Column in ``adata.obs`` for cluster labels
                (default ``"leiden"``).
            pseudotime_col: Column in ``adata.obs`` for pseudotime values
                (default ``"dpt_pseudotime"``).
        """
        color = annotation_key or cluster_key

        # --- Trajectory ---
        if trajectory and "X_diffmap" in adata.obsm:
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.diffmap(adata, color=color, components=["2, 3"],
                          show=False, ax=ax, title="Diffusion Map")
            self._save("advanced_diffmap.png")

            if pseudotime_col in adata.obs.columns:
                fig, ax = plt.subplots(figsize=(8, 6))
                sc.pl.diffmap(adata, color=pseudotime_col,
                              components=["2, 3"], show=False, ax=ax,
                              title="Pseudotime")
                self._save("advanced_pseudotime.png")

        # --- CNV ---
        if cnv:
            self._plot_cnv(adata, annotation_key=annotation_key,
                           cluster_key=cluster_key)

        logger.info("Advanced analysis plots saved to %s", self.plot_dir)

    def _plot_cnv(
        self,
        adata: AnnData,
        annotation_key: Optional[str] = None,
        cluster_key: str = "leiden",
    ) -> None:
        """Generate CNV-specific plots (requires infercnvpy).

        Args:
            adata: AnnData with CNV results from infercnvpy.
            annotation_key: Column in ``adata.obs`` for cell type annotation.
                None falls back to *cluster_key*.
            cluster_key: Column in ``adata.obs`` for cluster labels.
        """
        try:
            import infercnvpy as cnv
        except ImportError:
            logger.warning("infercnvpy not installed, skipping CNV plots")
            return

        groupby = annotation_key or cluster_key

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
        if annotation_key:
            cnv.pl.umap(adata, color=annotation_key, ax=axes[1, 0], show=False)
            axes[1, 0].set_title(annotation_key)
        fig.suptitle("CNV Analysis", fontsize=14, y=0.98)
        self._save("advanced_cnv_umap.png")

        # 3. CNV heatmap by cnv_leiden
        cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden",
                                  dendrogram=True, show=False)
        self._save("advanced_cnv_heatmap_leiden.png")

    # ------------------------------------------------------------------
    # Differential expression
    # ------------------------------------------------------------------
    def plot_de(
        self,
        adata: AnnData,
        group: str,
        volcano_group: Optional[str] = None,
    ) -> None:
        """Differential expression: rank_genes overview, volcano plot.

        Generates:
            - de_rank_genes_groups.png: top DEGs per group.
            - de_rank_genes_dotplot.png: DEG dotplot per group.
            - de_volcano_<group>.png: volcano plot for one group vs rest.

        Args:
            adata: AnnData with ``rank_genes_groups`` in ``adata.uns``.
            group: Grouping key used for DE (e.g. ``"leiden"`` or
                ``"condition"``).
            volcano_group: Specific group category to plot in the volcano.
                None uses the first category in ``adata.obs[group]``.
        """
        # 1. Rank genes groups overview
        sc.pl.rank_genes_groups(adata, n_genes=10, show=False)
        self._save("de_rank_genes_groups.png")

        # 2. Rank genes groups dotplot
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
        self._save("de_rank_genes_dotplot.png")

        # 3. Volcano plot
        vg = volcano_group
        if vg is None:
            try:
                vg = adata.obs[group].cat.categories[0]
            except (KeyError, IndexError) as exc:
                logger.warning("Cannot determine volcano group: %s", exc)
                vg = None

        if vg is not None:
            try:
                self._volcano(adata, group=vg)
            except Exception as exc:
                logger.warning("Volcano plot failed: %s", exc)

        logger.info("DE plots saved to %s", self.plot_dir)

    def _volcano(self, adata: AnnData, group: str) -> None:
        """Draw a volcano plot for one group vs rest.

        Args:
            adata: AnnData with ``rank_genes_groups`` computed.
            group: The specific group category to plot.
        """
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
        """Parse marker TSV into {cell_type: [gene, ...]} (genes present in adata).

        Args:
            marker_file: Path to a tab-separated file with ``cell_type`` and
                ``markers`` (comma-separated gene names) columns.
            adata: AnnData whose ``var.index`` is used to filter valid genes.

        Returns:
            Dict mapping cell type names to lists of valid marker gene names.
        """
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
