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
    # ------------------------------------------------------------------
    # HVG
    # ------------------------------------------------------------------
    def plot_hvg(
        self,
        adata: AnnData,
        n_top_genes: int = 3000,
    ) -> None:
        """Plot highly variable genes with all genes visible.

        Two subplots:
            Left:  mean expression vs normalized dispersion scatter
                   (grey = non-HVG, blue = HVG) with selected count.
            Right: sorted normalized dispersion elbow curve.

        Must be called BEFORE subsetting to HVGs so all genes are visible.

        Args:
            adata: AnnData with HVG metadata in ``.var`` (all genes).
            n_top_genes: Number of HVGs selected for downstream.
        """
        if "highly_variable" not in adata.var.columns:
            logger.warning("No HVG metadata found, skipping HVG plot")
            return

        hvg_mask = adata.var["highly_variable"].values
        has_means = "means" in adata.var.columns
        has_disp = "dispersions_norm" in adata.var.columns

        if not has_means or not has_disp:
            logger.warning("HVG metadata incomplete (means/dispersions_norm missing)")
            return

        means = adata.var["means"].values
        disp_norm = adata.var["dispersions_norm"].values

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # --- Left: mean vs normalized dispersion scatter ---
        ax = axes[0]
        ax.scatter(means[~hvg_mask], disp_norm[~hvg_mask],
                   c="lightgrey", s=2, alpha=0.3, label="Non-HVG",
                   rasterized=True)
        ax.scatter(means[hvg_mask], disp_norm[hvg_mask],
                   c="steelblue", s=4, alpha=0.7, label="HVG",
                   rasterized=True)

        # Threshold line at the dispersion cutoff
        hvg_dispersions = np.sort(disp_norm[hvg_mask])[::-1]
        actual_n = min(n_top_genes, len(hvg_dispersions))
        if actual_n > 0:
            threshold = hvg_dispersions[actual_n - 1]
            ax.axhline(y=threshold, color="red", linestyle="--",
                       alpha=0.7,
                       label=f"Cutoff ({actual_n} genes)")

        ax.set_xlabel("Mean expression")
        ax.set_ylabel("Normalized dispersion")
        ax.set_title(
            f"Highly Variable Genes ({int(hvg_mask.sum())} / {len(hvg_mask)})"
        )
        ax.legend(fontsize=8, loc="upper right")

        # --- Right: sorted dispersion elbow plot ---
        ax = axes[1]
        valid_mask = np.isfinite(disp_norm)
        sorted_disp = np.sort(disp_norm[valid_mask])[::-1]
        ax.plot(sorted_disp, linewidth=1, color="steelblue")

        if actual_n > 0 and actual_n < len(sorted_disp):
            ax.axvline(x=actual_n, color="red",
                        linestyle="--", alpha=0.7,
                        label=f"Selected: {actual_n}")
            ax.plot(actual_n, sorted_disp[actual_n], "ro", markersize=6)
            ax.legend(fontsize=8)

        ax.set_xlabel("Gene rank (by normalized dispersion)")
        ax.set_ylabel("Normalized dispersion")
        ax.set_title("Sorted Dispersion (Elbow)")
        ax.set_xlim(
            0,
            min(len(sorted_disp), max(actual_n * 2, 10000)),
        )

        fig.tight_layout()
        self._save("cluster_highly_variable_genes.png")
        logger.info("HVG plot saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Cluster
    # ------------------------------------------------------------------
    def plot_cluster(
        self,
        adata: AnnData,
        cluster_key: str = "leiden",
        sample_key: str = "sample_id",
    ) -> None:
        """Clustering results: UMAP by cluster & sample.

        Note: PCA variance ratio and HVG plots are generated separately
        by ``plot_pca_variance()`` and ``plot_hvg()`` respectively,
        which are called earlier in the pipeline.

        Generates:
            - cluster_umap_<cluster_key>.png: UMAP coloured by clusters.
            - cluster_umap_<sample_key>.png: UMAP coloured by sample.

        Args:
            adata: Clustered AnnData with ``X_umap`` and cluster labels
                in ``adata.obs[cluster_key]``.
            cluster_key: Column in ``adata.obs`` for cluster labels
                (default ``"leiden"``).
            sample_key: Column in ``adata.obs`` for sample identifiers
                (default ``"sample_id"``).
        """
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
    # PCA variance
    # ------------------------------------------------------------------
    def plot_pca_variance(
        self,
        adata: AnnData,
        n_pcs: int = 50,
        auto_n_pcs: bool = False,
        detect_diag: Optional[Dict] = None,
    ) -> None:
        """Plot PCA variance ratio.

        When *auto_n_pcs* is True and *detect_diag* is provided, generates
        two subplots:
            Left:  scree plot with each PC labeled.
            Right: sliding-window mean Δ curve + threshold + elbow marker.

        When *auto_n_pcs* is False, generates one scree plot with the
        selected *n_pcs* PCs annotated.

        Args:
            adata: AnnData with PCA results in ``.uns["pca"]``.
            n_pcs: Number of PCs selected for downstream analysis.
            auto_n_pcs: Whether auto-detection was used.
            detect_diag: Diagnostics dict from ``detect_n_pcs()`` with
                keys ``window_mean_x``, ``window_mean_y``, ``threshold``,
                ``elbow_pc``.  Only used when *auto_n_pcs* is True.
        """
        if "pca" not in adata.uns:
            logger.warning("No PCA results found, skipping PCA variance plot")
            return

        variance_ratio = adata.uns["pca"]["variance_ratio"]
        n_pcs_total = len(variance_ratio)

        has_diag = (
            auto_n_pcs
            and detect_diag is not None
            and len(detect_diag.get("delta", [])) > 0
        )

        if has_diag:
            fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 5))
        else:
            fig, ax_left = plt.subplots(figsize=(8, 5))

        x = np.arange(1, n_pcs_total + 1)

        # --- Scree plot with PC labels ---
        ax_left.scatter(x, variance_ratio, s=20, c="steelblue", zorder=3)
        ax_left.plot(x, variance_ratio, linewidth=0.8, alpha=0.4, c="steelblue")

        for xi, yi in zip(x, variance_ratio):
            ax_left.annotate(f"PC{xi}", (xi, yi), fontsize=6,
                             ha="center", va="bottom", alpha=0.7)

        if n_pcs <= n_pcs_total:
            ax_left.axvline(x=n_pcs, color="red", linestyle="--", alpha=0.7,
                            label=f"Selected: {n_pcs} PCs")
            ax_left.legend(fontsize=8)

        ax_left.set_xlabel("ranking")
        ax_left.set_ylabel("variance ratio")
        ax_left.set_title("variance ratio")

        # --- Right: relative change rate + stability ---
        if has_diag:
            delta = detect_diag["delta"]       # now relative change
            threshold = detect_diag["threshold"]
            std_threshold = detect_diag.get("std_threshold", 0.03)
            elbow_pc = detect_diag["elbow_pc"]
            wm_x = detect_diag.get("window_mean_x", np.array([]))
            wm_y = detect_diag.get("window_mean_y", np.array([]))
            ws_y = detect_diag.get("window_std_y", np.array([]))

            delta_x = np.arange(2, len(delta) + 2)

            ax_right.plot(delta_x, delta * 100, linewidth=0.8, c="steelblue",
                          alpha=0.5, label="relative change (%)")
            if len(wm_x) > 0:
                ax_right.plot(wm_x, wm_y * 100, linewidth=1.5, c="steelblue",
                              marker="o", markersize=3,
                              label=f"window median")
            ax_right.axhline(y=threshold * 100, color="red", linestyle="--",
                             alpha=0.7,
                             label=f"rel threshold ({threshold*100:.0f}%)")

            # Mark elbow point
            if 2 <= elbow_pc <= len(delta) + 1:
                ax_right.plot(elbow_pc, delta[elbow_pc - 2] * 100, "ro",
                              markersize=8, zorder=5,
                              label=f"elbow: PC{elbow_pc}")

            ax_right.set_xlabel("PC index")
            ax_right.set_ylabel("relative change (%)")
            ax_right.set_title("Per-PC relative change (plateau = low + stable)")
            ax_right.legend(fontsize=7)

        fig.tight_layout()
        self._save("cluster_pca_variance_ratio.png")
        logger.info("PCA variance plot saved to %s", self.plot_dir)

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
                fig, ax = plt.subplots(figsize=(10, 6))
                self._umap(adata, col, ax=ax, legend_loc="right margin",
                           legend_fontsize=9, legend_fontoutline=1)
                self._save(f"annotate_umap_{col}.png")

        # 4. Confidence score
        if score_col in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(8, 6))
            self._umap(adata, score_col, ax=ax,
                       title=f"Confidence Score ({score_col})")
            self._save(f"annotate_{score_col}.png")

        logger.info("Annotation plots saved to %s", self.plot_dir)

    # ------------------------------------------------------------------
    # Marker expression on UMAP (post-annotation)
    # ------------------------------------------------------------------
    def plot_markers(
        self,
        adata: AnnData,
        annotations: Dict[str, Dict],
        cluster_key: str = "leiden",
    ) -> int:
        """Plot marker gene expression on UMAP, grouped by cell type.

        For each unique cell type, plots ALL key_markers from all its
        clusters on the UMAP with cluster centroid labels.  Multiple
        clusters sharing the same cell type are highlighted together.

        Also generates:
        - A combined panel (one top marker per cell type).
        - A dotplot grouped by cell type.

        Output directory: ``<plot_dir>/marker/``.

        Args:
            adata: Annotated AnnData (must have ``X_umap`` in ``obsm``).
            annotations: Dict mapping cluster_id -> annotation dict
                (must contain ``cell_type`` and ``key_markers``).
            cluster_key: Column in ``adata.obs`` for cluster labels.

        Returns:
            Total number of unique markers actually plotted.
        """
        marker_dir = os.path.join(self.plot_dir, "marker")
        os.makedirs(marker_dir, exist_ok=True)

        # Use raw for expression values
        if adata.raw is not None:
            adata_raw = adata.raw.to_adata()
            adata_raw.obsm = adata.obsm
            adata_raw.uns = adata.uns
        else:
            adata_raw = adata

        all_genes = set(adata_raw.var_names)
        umap = adata.obsm["X_umap"]

        # Compute cluster centroids
        centroids: Dict[str, np.ndarray] = {}
        for cid in sorted(annotations.keys(), key=lambda x: int(x)):
            mask = adata.obs[cluster_key] == cid
            if mask.sum() > 0:
                centroids[cid] = umap[mask.values].mean(axis=0)

        # Group clusters by cell type
        ct_to_clusters: Dict[str, List[str]] = {}
        for cid in sorted(annotations.keys(), key=lambda x: int(x)):
            ct = annotations[cid]["cell_type"]
            ct_to_clusters.setdefault(ct, []).append(cid)

        def _get_expr(gene: str) -> np.ndarray:
            x = adata_raw[:, gene].X
            if hasattr(x, "toarray"):
                return x.toarray().flatten()
            return np.asarray(x).flatten()

        # ── Cluster palette helpers ──────────────────────────────────────
        # `adata.uns[f"{cluster_key}_colors"]` is the canonical palette written
        # by scanpy when ``sc.pl.umap(..., color=cluster_key)`` is first
        # invoked.  We reuse it for the "Cluster ref" subplot so the colours
        # match the main UMAP.  Falls back to matplotlib's tab20 cycle when
        # the palette is missing or shorter than the number of clusters.
        _tab_palette = tuple(
            plt.get_cmap("tab20")(i) for i in range(20)
        )

        def _to_hex(col) -> str:
            """Coerce a matplotlib colour (str or RGBA tuple) to a hex string."""
            if isinstance(col, str):
                return col
            try:
                return "#" + "".join(
                    f"{int(round(c * 255)):02x}" for c in col[:3]
                )
            except (TypeError, ValueError):
                return "#1f77b4"

        def _resolve_cluster_colors(
            _adata, _key: str,
        ) -> Dict[str, str]:
            obs_col = _adata.obs[_key]
            cats = list(obs_col.cat.categories) \
                if hasattr(obs_col, "cat") \
                else sorted(obs_col.astype(str).unique())
            stored = _adata.uns.get(f"{_key}_colors")
            if stored is None or len(stored) < len(cats):
                stored = [_tab_palette[i % len(_tab_palette)] for i in range(len(cats))]
            return {cat: _to_hex(stored[i]) for i, cat in enumerate(cats)}

        def _default_cluster_color(_cid: str) -> str:
            """Tab-palette fallback when `_cid` is not in the resolved map."""
            try:
                idx = int(_cid)
            except (TypeError, ValueError):
                idx = 0
            return _to_hex(_tab_palette[idx % len(_tab_palette)])

        def _annotate_centroids(ax: plt.Axes, highlight_cids: List[str]) -> None:
            for _cid, (cx, cy) in centroids.items():
                is_target = _cid in highlight_cids
                ax.annotate(
                    _cid, (cx, cy),
                    fontsize=9 if is_target else 7,
                    fontweight="bold" if is_target else "normal",
                    color="red" if is_target else "black",
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15",
                              facecolor="white", alpha=0.0, edgecolor="none"),
                )

        # ── Per-cell-type marker panels ──
        total_unique_markers = 0
        seen_genes: set = set()
        all_top_markers: List[str] = []   # one top marker per cell type
        all_top_cts: List[str] = []

        for ct, cluster_ids in ct_to_clusters.items():
            # Collect ALL unique markers from all clusters of this cell type
            markers: List[str] = []
            seen_ct_genes: set = set()
            for cid in cluster_ids:
                for g in annotations[cid].get("key_markers", []):
                    if g in all_genes and g not in seen_ct_genes:
                        markers.append(g)
                        seen_ct_genes.add(g)

            if not markers:
                continue

            # Subplot grid: markers + 1 UMAP reference, up to 5 cols
            n = len(markers)
            total_panels = n + 1  # +1 for UMAP cluster reference
            ncols = min(5, total_panels)
            nrows = (total_panels + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(4 * ncols, 4 * nrows))
            axes_flat = np.asarray(axes).flatten() if total_panels > 1 else [axes]

            for ax, gene in zip(axes_flat, markers):
                expr = _get_expr(gene)
                ax.scatter(umap[:, 0], umap[:, 1], c="#e0e0e0",
                           s=2, alpha=0.3, edgecolors="none")
                sort_idx = np.argsort(expr)
                sm = ax.scatter(
                    umap[sort_idx, 0], umap[sort_idx, 1],
                    c=expr[sort_idx], cmap="Reds", s=3,
                    alpha=0.7, edgecolors="none",
                )
                _annotate_centroids(ax, cluster_ids)
                ax.set_title(gene, fontsize=11, fontweight="bold")
                ax.set_xticks([])
                ax.set_yticks([])
                plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04,
                             label="expression")

            # Last subplot: UMAP cluster reference — colour each cluster with
            # the same palette as the main UMAP (adata.uns[leiden_colors]).
            # Falling back to tab10 cycle keeps things readable if the colours
            # weren't stored.
            ax_ref = axes_flat[n]
            ax_ref.scatter(umap[:, 0], umap[:, 1], c="#e0e0e0",
                           s=2, alpha=0.3, edgecolors="none")
            cluster_color_map = _resolve_cluster_colors(adata, cluster_key)
            for _cid in cluster_ids:
                mask = adata.obs[cluster_key] == _cid
                if mask.sum() == 0:
                    continue
                ax_ref.scatter(umap[mask.values, 0], umap[mask.values, 1],
                               c=cluster_color_map.get(
                                   _cid, _default_cluster_color(_cid)),
                               s=12, alpha=0.85, edgecolors="none",
                               linewidths=0,
                               label=f"C{_cid}")
            for _cid, (cx, cy) in centroids.items():
                is_target = _cid in cluster_ids
                ax_ref.annotate(
                    _cid, (cx, cy),
                    fontsize=9 if is_target else 7,
                    fontweight="bold" if is_target else "normal",
                    color="red" if is_target else "black",
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15",
                              facecolor="white", alpha=0.0, edgecolor="none"),
                )
            ax_ref.set_title("Cluster ref", fontsize=11, fontweight="bold")
            ax_ref.set_xticks([])
            ax_ref.set_yticks([])
            ax_ref.legend(loc="best", fontsize=7, framealpha=0.8,
                          labelcolor="black")

            # Hide unused axes
            for j in range(total_panels, len(axes_flat)):
                axes_flat[j].set_visible(False)

            # Title: cell type + cluster ids
            cid_str = ", ".join(f"C{c}" for c in cluster_ids)
            fig.suptitle(f"{ct}  [{cid_str}]",
                         fontsize=14, fontweight="bold", y=1.02)
            plt.tight_layout()
            safe_ct = ct.replace("/", "_").replace(" ", "_")
            fig.savefig(os.path.join(marker_dir, f"{safe_ct}.png"),
                        dpi=self.dpi, bbox_inches="tight")
            plt.close(fig)

            for g in markers:
                if g not in seen_genes:
                    seen_genes.add(g)
                    total_unique_markers += 1

            # Collect top marker for combined panel
            for cid in cluster_ids:
                for g in annotations[cid].get("key_markers", []):
                    if g in all_genes and g not in all_top_markers:
                        all_top_markers.append(g)
                        all_top_cts.append(ct)
                        break
                if len(all_top_cts) > 0 and all_top_cts[-1] == ct:
                    break

        # ── Combined panel: one top marker per cell type ──
        if all_top_markers:
            n = len(all_top_markers)
            ncols = min(5, n)
            nrows = (n + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(4 * ncols, 4 * nrows))
            axes_flat = np.asarray(axes).flatten() if n > 1 else [axes]

            for i, (gene, ct, ax) in enumerate(
                    zip(all_top_markers, all_top_cts, axes_flat)):
                expr = _get_expr(gene)
                ax.scatter(umap[:, 0], umap[:, 1], c="#e0e0e0",
                           s=1, alpha=0.2, edgecolors="none")
                sort_idx = np.argsort(expr)
                ax.scatter(
                    umap[sort_idx, 0], umap[sort_idx, 1],
                    c=expr[sort_idx], cmap="Reds", s=2,
                    alpha=0.7, edgecolors="none",
                )
                # Highlight all clusters for this cell type
                cids = ct_to_clusters[ct]
                _annotate_centroids(ax, cids)
                cid_str = ",".join(cids)
                ax.set_title(f"C{cid_str} {gene}", fontsize=8)
                ax.set_xticks([])
                ax.set_yticks([])

            for j in range(i + 1, len(axes_flat)):
                axes_flat[j].set_visible(False)

            fig.suptitle("Top marker expression on UMAP (per cell type)",
                         fontsize=14, fontweight="bold")
            plt.tight_layout()
            fig.savefig(os.path.join(marker_dir, "all_cell_types_top_marker.png"),
                        dpi=self.dpi, bbox_inches="tight")
            plt.close(fig)

        # ── Dotplot grouped by cell_type ──
        # Use ALL markers, grouped by cell type label
        dotplot_markers: List[str] = []
        seen_dotplot: set = set()
        # Build ordered cell type labels for groupby
        type_count: Dict[str, int] = {}
        for ann in annotations.values():
            type_count[ann["cell_type"]] = type_count.get(ann["cell_type"], 0) + 1

        label_map: Dict[str, str] = {}
        for cid in sorted(annotations.keys(), key=lambda x: int(x)):
            ct = annotations[cid]["cell_type"]
            if type_count[ct] > 1:
                label_map[cid] = f"{ct} (C{cid})"
            else:
                label_map[cid] = ct

        adata.obs["_marker_ct_label"] = adata.obs[cluster_key].map(label_map).astype("category")
        adata_raw.obs["_marker_ct_label"] = adata.obs["_marker_ct_label"]

        for ct, cluster_ids in ct_to_clusters.items():
            for cid in cluster_ids:
                for g in annotations[cid].get("key_markers", []):
                    if g in all_genes and g not in seen_dotplot:
                        dotplot_markers.append(g)
                        seen_dotplot.add(g)

        if dotplot_markers:
            sc.pl.dotplot(
                adata_raw, var_names=dotplot_markers,
                groupby="_marker_ct_label", standard_scale="var",
                show=False,
            )
            plt.savefig(os.path.join(marker_dir, "marker_dotplot.png"),
                        dpi=self.dpi, bbox_inches="tight")
            plt.close("all")

        # Clean up temporary columns from adata.obs
        for col in ("_marker_ct_label", "_marker_cell_type_label"):
            if col in adata.obs.columns:
                adata.obs.drop(columns=[col], inplace=True)

        logger.info("Marker plots saved to %s (%d unique markers, %d cell types)",
                     marker_dir, total_unique_markers, len(ct_to_clusters))
        return total_unique_markers

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
