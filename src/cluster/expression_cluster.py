"""
Gene expression matrix clustering.

Three-class architecture:
  ExpressionPreprocessor — load, filter, transform
  DistanceCalculator     — distance metrics + clustering algorithms
  ClusterPlotter         — heatmap, UMAP, dendrogram, elbow

Typical usage:
    prep = ExpressionPreprocessor(min_mean=1.0).load("matrix.tsv")
    calc = DistanceCalculator(metric="cosine")
    result = calc.hierarchical(prep.scaled, n_clusters=5)
    plotter = ClusterPlotter()
    plotter.heatmap(prep.raw, result.labels, "heatmap.png")
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ClusterResult:
    """Container for clustering output.

    Attributes
    ----------
    labels : np.ndarray
        Cluster label for each gene (0-indexed for kmeans, 1-indexed for others).
    Z : np.ndarray | None
        Linkage matrix (only for hierarchical), used by dendrogram.
    method : str
        Clustering method that produced this result.
    metric : str
        Distance metric used.
    n_clusters : int
        Number of clusters found.
    """
    labels: np.ndarray
    Z: np.ndarray | None = None
    method: str = ""
    metric: str = ""
    n_clusters: int = 0

    def to_dataframe(self, gene_names: pd.Index) -> pd.DataFrame:
        """Return DataFrame with columns [Gene, Cluster], sorted by Cluster."""
        df = pd.DataFrame({"Gene": gene_names, "Cluster": self.labels})
        return df.sort_values("Cluster").reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Cluster size summary table."""
        _, counts = np.unique(self.labels, return_counts=True)
        clusters = sorted(set(self.labels))
        return pd.DataFrame({"Cluster": clusters, "n_genes": counts})


# ---------------------------------------------------------------------------
# Class 1: ExpressionPreprocessor
# ---------------------------------------------------------------------------

class ExpressionPreprocessor:
    """Load, filter, and transform gene expression matrices.

    Parameters
    ----------
    min_mean : float
        Minimum mean expression to retain a gene.
    min_samples : int
        Minimum number of samples exceeding min_mean.
    log_transform : bool
        Apply log2(x + 1) before scaling.
    quantile_normalize : bool
        Quantile-normalize samples after log transform (removes
        distributional batch effects between datasets).
    combat : bool
        Apply ComBat batch correction (empirical Bayes). Requires
        batch labels set via from_dataframe(df, batch=...) or
        set_batch(). Overrides quantile_normalize when both are True.
    zscore : bool
        Z-score scale each gene (row) after log transform.

    Attributes (after processing)
    -----------------------------
    scaled : pd.DataFrame
        Z-scored matrix ready for clustering (genes x samples).
    raw : pd.DataFrame
        Filtered but unscaled matrix (for heatmap display).
    """

    def __init__(
        self,
        min_mean: float = 1.0,
        min_samples: int = 1,
        log_transform: bool = True,
        quantile_normalize: bool = False,
        combat: bool = False,
        zscore: bool = True,
    ):
        self.min_mean = min_mean
        self.min_samples = min_samples
        self.log_transform = log_transform
        self.quantile_normalize = quantile_normalize
        self.combat = combat
        self.zscore = zscore
        self._df: pd.DataFrame | None = None
        self._scaled: pd.DataFrame | None = None
        self._raw: pd.DataFrame | None = None
        self._batch: pd.Series | None = None

    # ---- loading ----

    def load(
        self,
        path: str | Path,
        sep: str = "\t",
        index_col: int = 0,
        transpose: bool = False,
    ) -> ExpressionPreprocessor:
        """Load expression matrix from file.

        Parameters
        ----------
        path : str or Path
            Path to TSV/CSV matrix.
        sep : str
            Column separator.
        index_col : int
            Column to use as row index.
        transpose : bool
            Transpose after loading (when rows are samples).

        Returns
        -------
        self
        """
        self._df = pd.read_csv(path, sep=sep, index_col=index_col)
        if transpose:
            self._df = self._df.T
        logger.info(f"Loaded matrix: {self._df.shape[0]} genes x {self._df.shape[1]} samples")
        self._process()
        return self

    def from_dataframe(
        self, df: pd.DataFrame, batch: pd.Series | None = None,
    ) -> ExpressionPreprocessor:
        """Load from an existing DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Expression matrix (genes x samples).
        batch : pd.Series, optional
            Batch labels for each sample (index must match df.columns).
            Required for ComBat correction.

        Returns
        -------
        self
        """
        self._df = df.copy()
        if batch is not None:
            self._batch = batch
        logger.info(f"Loaded matrix: {self._df.shape[0]} genes x {self._df.shape[1]} samples")
        self._process()
        return self

    def set_batch(self, batch: pd.Series) -> ExpressionPreprocessor:
        """Set batch labels for ComBat correction.

        Parameters
        ----------
        batch : pd.Series
            Batch label per sample. Index must match columns of the matrix.

        Returns
        -------
        self
        """
        self._batch = batch
        return self

    # ---- internal pipeline ----

    def _process(self) -> None:
        """Run filter → log → quantile_normalize → z-score pipeline."""
        df = self._filter_low_expression(self._df, self.min_mean, self.min_samples)

        if self.log_transform:
            df = np.log2(df + 1)
            logger.info("Applied log2(x + 1) transform")

        if self.combat:
            if self._batch is None:
                raise ValueError("ComBat requires batch labels. Use from_dataframe(df, batch=...) or set_batch()")
            df = self._combat_correct(df, self._batch)
            logger.info(f"Applied ComBat batch correction ({self._batch.nunique()} batches)")
        elif self.quantile_normalize:
            df = self._quantile_normalize(df)
            logger.info("Applied quantile normalization")

        self._raw = df.copy()

        if self.zscore:
            df = df.subtract(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)
            df = df.fillna(0)
            logger.info("Applied z-score scaling per gene")

        self._scaled = df

    @staticmethod
    def _quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
        """Quantile-normalize samples so all share the same value distribution.

        For each sample, rank its values, then replace each value with the
        mean of all samples' values at that rank. This eliminates
        distributional differences between datasets (batch effect).

        Parameters
        ----------
        df : pd.DataFrame
            Expression matrix (genes x samples).

        Returns
        -------
        pd.DataFrame
            Quantile-normalized matrix with same shape and index/columns.
        """
        # sort each column, compute rank means across all columns
        sorted_vals = np.sort(df.values, axis=0)
        rank_means = sorted_vals.mean(axis=1)

        # rank each column (1-based), then map ranks to means
        ranks = df.rank(method="min").astype(int).values
        # clip to valid range (rank is 1-based, index is 0-based)
        ranks = np.clip(ranks, 1, len(rank_means))
        df_qn = pd.DataFrame(
            rank_means[ranks - 1],
            index=df.index, columns=df.columns,
        )
        return df_qn

    @staticmethod
    def _combat_correct(
        df: pd.DataFrame, batch: pd.Series,
    ) -> pd.DataFrame:
        """ComBat batch correction (empirical Bayes).

        Implements the ComBat algorithm (Johnson et al. 2007) in pure
        numpy/pandas. Handles both additive (location) and multiplicative
        (scale) batch effects with empirical Bayes shrinkage.

        Parameters
        ----------
        df : pd.DataFrame
            Expression matrix (genes x samples), typically after log2.
        batch : pd.Series
            Batch label per sample. Index must match df.columns.

        Returns
        -------
        pd.DataFrame
            Batch-corrected matrix with same shape.
        """
        batch = batch.loc[df.columns]
        batches = batch.unique()
        n_batch = len(batches)
        n_genes = df.shape[0]
        n_samples = df.shape[1]

        if n_batch < 2:
            logger.warning("Only 1 batch; skipping correction")
            return df

        # design matrix (intercept only, no biological covariates)
        batch_design = pd.get_dummies(batch).loc[df.columns].values.astype(float)  # (n_samples, n_batch)

        # --- 1. OLS estimation of batch effects ---
        B_hat = np.linalg.lstsq(batch_design, df.values.T, rcond=None)[0]  # (n_batch, n_genes)
        grand_mean = B_hat.mean(axis=0)  # (n_genes,)

        # residual variance per gene
        fitted = batch_design @ B_hat  # (n_samples, n_genes)
        var_pooled = ((df.values.T - fitted) ** 2).mean(axis=0)  # (n_genes,)

        # --- 2. Standardize ---
        # batch-specific mean and std
        gamma_hat = B_hat - grand_mean  # (n_batch, n_genes) additive effect
        delta_hat = np.zeros((n_batch, n_genes))
        for i, b in enumerate(batches):
            idx = batch == b
            batch_data = df.values[:, idx.values].T  # (n_samples_in_batch, n_genes)
            delta_hat[i] = batch_data.std(axis=0, ddof=1)

        # avoid division by zero
        delta_hat = np.maximum(delta_hat, 1e-8)

        # --- 3. Empirical Bayes shrinkage (method of moments) ---
        gamma_star = np.zeros_like(gamma_hat)
        delta_star = np.zeros_like(delta_hat)

        for i, b in enumerate(batches):
            n_b = (batch == b).sum()

            # prior for gamma: Normal(g_bar, tau^2)
            g_bar = gamma_hat[i].mean()
            tau2 = gamma_hat[i].var()

            # single-sample batch: no shrinkage possible, use raw estimate
            if n_b <= 1 or tau2 == 0 or var_pooled.min() == 0:
                gamma_star[i] = gamma_hat[i]
                delta_star[i] = delta_hat[i]
                continue

            # prior for delta: InverseGamma(lambda, rate)
            # method of moments on log(delta^2)
            log_d2 = np.log(delta_hat[i] ** 2)
            d_bar = log_d2.mean()
            v = log_d2.var()

            # posterior mean for gamma (shrinkage toward g_bar)
            gamma_star[i] = (n_b * gamma_hat[i] / var_pooled + g_bar / tau2) / \
                            (n_b / var_pooled + 1 / tau2)

            # posterior for delta (shrinkage on log scale)
            v_star = 1.0 / (n_b / 2.0 + 1.0 / v) if v > 0 else 0
            d_star = (n_b / 2.0 * np.log(delta_hat[i] ** 2) + d_bar / v) * v_star if v > 0 else np.log(delta_hat[i] ** 2)
            delta_star[i] = np.exp(d_star / 2.0)

        # --- 4. Apply correction ---
        # ComBat formula:
        # y_corrected = (y - gamma_star) / delta_star * sqrt(var_pooled) + grand_mean
        df_out = df.values.copy().astype(float)
        for i, b in enumerate(batches):
            idx = batch == b
            idx_arr = np.where(idx.values)[0]
            for j in idx_arr:
                df_out[:, j] = (df.values[:, j] - gamma_star[i]) / delta_star[i] * np.sqrt(var_pooled) + grand_mean

        # replace any NaN/inf from edge cases
        df_out = np.nan_to_num(df_out, nan=0.0, posinf=0.0, neginf=0.0)
        return pd.DataFrame(df_out, index=df.index, columns=df.columns)

    @staticmethod
    def _filter_low_expression(
        df: pd.DataFrame,
        min_mean: float = 1.0,
        min_samples: int = 1,
    ) -> pd.DataFrame:
        """Remove genes with low mean expression."""
        mask = (df > min_mean).sum(axis=1) >= min_samples
        n_before = df.shape[0]
        df_filtered = df.loc[mask]
        logger.info(
            f"Filtered: {n_before} -> {df_filtered.shape[0]} genes "
            f"(min_mean={min_mean}, min_samples={min_samples})"
        )
        return df_filtered

    # ---- public properties ----

    @property
    def scaled(self) -> pd.DataFrame:
        """Z-scored matrix for clustering (genes x samples)."""
        if self._scaled is None:
            raise RuntimeError("No data loaded. Call load() or from_dataframe() first.")
        return self._scaled

    @property
    def raw(self) -> pd.DataFrame:
        """Filtered, log-transformed (if enabled), unscaled matrix."""
        if self._raw is None:
            raise RuntimeError("No data loaded. Call load() or from_dataframe() first.")
        return self._raw

    @property
    def gene_names(self) -> pd.Index:
        """Gene names (row index) of the processed matrix."""
        return self.scaled.index

    @property
    def sample_names(self) -> pd.Index:
        """Sample names (column index) of the processed matrix."""
        return self.scaled.columns


# ---------------------------------------------------------------------------
# Class 2: DistanceCalculator
# ---------------------------------------------------------------------------

class DistanceCalculator:
    """Compute distances and perform clustering.

    Parameters
    ----------
    metric : str
        Distance metric. One of:
        euclidean, maximum, manhattan, canberra, binary,
        minkowski, cosine.
    minkowski_p : float
        p parameter for minkowski distance (default 2.0).

    Supported metrics and behavior per clustering method:

    =============  ==============  ==============  ==============
    Metric         Hierarchical    KMeans          Leiden
    =============  ==============  ==============  ==============
    euclidean      native          native          native
    maximum        native          fallback*       fallback*
    manhattan      native          fallback*       native
    canberra       native          fallback*       fallback*
    binary         native          fallback*       fallback*
    minkowski      native          fallback*       fallback*
    cosine         native          fallback*       native
    =============  ==============  ==============  ==============

    * fallback = hierarchical (average linkage) is used instead.
    """

    SUPPORTED_METRICS = (
        "euclidean", "maximum", "manhattan", "canberra",
        "binary", "minkowski", "cosine",
    )

    # scipy pdist uses different names for some metrics
    _METRIC_ALIASES = {
        "manhattan": "cityblock",
    }

    def __init__(
        self,
        metric: str = "euclidean",
        minkowski_p: float = 2.0,
    ):
        if metric not in self.SUPPORTED_METRICS:
            raise ValueError(
                f"Unsupported metric '{metric}'. "
                f"Choose from: {', '.join(self.SUPPORTED_METRICS)}"
            )
        self.metric = metric
        self.minkowski_p = minkowski_p

    # ---- low-level distance ----

    def pairwise(self, df: pd.DataFrame) -> np.ndarray:
        """Compute condensed pairwise distance matrix.

        Parameters
        ----------
        df : pd.DataFrame
            Expression matrix (genes x samples), typically preprocessed.scaled.

        Returns
        -------
        np.ndarray
            Condensed distance vector (scipy pdist format).
        """
        pdist_metric = self._METRIC_ALIASES.get(self.metric, self.metric)
        kwargs = {}
        if self.metric == "minkowski":
            kwargs["p"] = self.minkowski_p
        return pdist(df.values, metric=pdist_metric, **kwargs)

    # ---- clustering methods ----

    def hierarchical(
        self,
        df: pd.DataFrame,
        n_clusters: int = 5,
        method: str = "ward",
    ) -> ClusterResult:
        """Hierarchical agglomerative clustering.

        Parameters
        ----------
        df : pd.DataFrame
            Scaled expression matrix (genes x samples).
        n_clusters : int
            Number of clusters.
        method : str
            Linkage method: ward, complete, average, single.
            Ward requires euclidean; if other metric is set,
            metric is overridden to euclidean with a warning.

        Returns
        -------
        ClusterResult
        """
        metric = self.metric
        if method == "ward" and metric != "euclidean":
            logger.warning(
                f"Ward linkage requires euclidean; overriding metric from '{metric}'"
            )
            metric = "euclidean"

        pdist_metric = self._METRIC_ALIASES.get(metric, metric)
        kwargs = {}
        if metric == "minkowski":
            kwargs["p"] = self.minkowski_p

        dist = pdist(df.values, metric=pdist_metric, **kwargs)
        Z = linkage(dist, method=method)
        labels = fcluster(Z, t=n_clusters, criterion="maxclust")

        logger.info(
            f"Hierarchical clustering: {n_clusters} clusters, "
            f"method={method}, metric={metric}"
        )
        return ClusterResult(
            labels=labels, Z=Z, method=f"hierarchical({method})",
            metric=metric, n_clusters=n_clusters,
        )

    def kmeans(
        self,
        df: pd.DataFrame,
        n_clusters: int = 5,
        random_state: int = 42,
        n_init: int = 10,
    ) -> ClusterResult:
        """K-means clustering.

        Note: sklearn KMeans only supports euclidean. For non-euclidean
        metrics, automatically falls back to hierarchical (average linkage).

        Parameters
        ----------
        df : pd.DataFrame
            Scaled expression matrix (genes x samples).
        n_clusters : int
            Number of clusters.
        random_state : int
            Random seed.
        n_init : int
            Number of initializations.

        Returns
        -------
        ClusterResult
        """
        if self.metric != "euclidean":
            logger.warning(
                f"KMeans only supports euclidean; '{self.metric}' "
                f"triggers hierarchical (average) fallback."
            )
            return self.hierarchical(df, n_clusters=n_clusters, method="average")

        from sklearn.cluster import KMeans

        km = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=n_init)
        labels = km.fit_predict(df.values)
        logger.info(f"K-means: {n_clusters} clusters, inertia={km.inertia_:.2f}")
        return ClusterResult(
            labels=labels, method="kmeans", metric="euclidean",
            n_clusters=n_clusters,
        )

    def leiden(
        self,
        df: pd.DataFrame,
        resolution: float = 1.0,
        n_neighbors: int = 15,
        random_state: int = 42,
    ) -> ClusterResult:
        """Leiden graph-based clustering.

        Requires: pip install scanpy leidenalg

        Parameters
        ----------
        df : pd.DataFrame
            Scaled expression matrix (genes x samples).
        resolution : float
            Leiden resolution (higher = more clusters).
        n_neighbors : int
            Number of neighbors for KNN graph.
        random_state : int
            Random seed.

        Returns
        -------
        ClusterResult
        """
        try:
            import scanpy as sc
        except ImportError:
            raise ImportError(
                "Leiden clustering requires scanpy: pip install scanpy leidenalg"
            )

        scanpy_metric = self._METRIC_ALIASES.get(self.metric, self.metric)

        adata = sc.AnnData(df.values)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X", metric=scanpy_metric)
        sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
        labels = adata.obs["leiden"].astype(int).values
        n_found = len(set(labels))

        logger.info(
            f"Leiden: {n_found} clusters, resolution={resolution}, metric={self.metric}"
        )
        return ClusterResult(
            labels=labels, method="leiden", metric=self.metric,
            n_clusters=n_found,
        )

    def dbscan(
        self,
        df: pd.DataFrame,
        eps: float = 0.5,
        min_samples: int = 5,
    ) -> ClusterResult:
        """DBSCAN density-based clustering.

        Good for non-convex clusters and automatic outlier detection.
        Noise points are labeled -1.

        Parameters
        ----------
        df : pd.DataFrame
            Scaled expression matrix (genes x samples).
        eps : float
            Maximum distance between two samples in the same neighborhood.
        min_samples : int
            Minimum number of samples in a neighborhood for a core point.

        Returns
        -------
        ClusterResult
        """
        from sklearn.cluster import DBSCAN

        pdist_metric = self._METRIC_ALIASES.get(self.metric, self.metric)
        db = DBSCAN(eps=eps, min_samples=min_samples, metric=pdist_metric)
        labels = db.fit_predict(df.values)
        n_clusters = len(set(labels) - {-1})
        n_noise = (labels == -1).sum()

        logger.info(
            f"DBSCAN: {n_clusters} clusters, {n_noise} noise points, "
            f"eps={eps}, min_samples={min_samples}, metric={self.metric}"
        )
        return ClusterResult(
            labels=labels, method="dbscan", metric=self.metric,
            n_clusters=n_clusters,
        )

    def spectral(
        self,
        df: pd.DataFrame,
        n_clusters: int = 5,
        n_neighbors: int = 10,
        random_state: int = 42,
    ) -> ClusterResult:
        """Spectral clustering using graph Laplacian.

        Effective for non-convex clusters and captures manifold structure.
        Uses affinity matrix built from nearest neighbors.

        Parameters
        ----------
        df : pd.DataFrame
            Scaled expression matrix (genes x samples).
        n_clusters : int
            Number of clusters.
        n_neighbors : int
            Number of neighbors for affinity matrix construction.
        random_state : int
            Random seed.

        Returns
        -------
        ClusterResult
        """
        from sklearn.cluster import SpectralClustering

        # Build affinity matrix from pairwise distances
        pdist_metric = self._METRIC_ALIASES.get(self.metric, self.metric)
        dist = pdist(df.values, metric=pdist_metric)
        from scipy.spatial.distance import squareform
        dist_matrix = squareform(dist)

        # Convert distance to affinity (Gaussian kernel)
        sigma = np.median(dist)
        affinity = np.exp(-dist_matrix ** 2 / (2 * sigma ** 2))

        sc = SpectralClustering(
            n_clusters=n_clusters,
            affinity="precomputed",
            random_state=random_state,
            n_neighbors=n_neighbors,
        )
        labels = sc.fit_predict(affinity)

        logger.info(
            f"Spectral: {n_clusters} clusters, n_neighbors={n_neighbors}, "
            f"metric={self.metric}"
        )
        return ClusterResult(
            labels=labels, method="spectral", metric=self.metric,
            n_clusters=n_clusters,
        )

    def gmm(
        self,
        df: pd.DataFrame,
        n_components: int = 5,
        covariance_type: str = "full",
        random_state: int = 42,
    ) -> ClusterResult:
        """Gaussian Mixture Model clustering.

        Soft clustering with probabilistic assignments. Good for
        overlapping clusters and provides uncertainty estimates.

        Parameters
        ----------
        df : pd.DataFrame
            Scaled expression matrix (genes x samples).
        n_components : int
            Number of mixture components.
        covariance_type : str
            Covariance type: full, tied, diag, spherical.
        random_state : int
            Random seed.

        Returns
        -------
        ClusterResult
        """
        from sklearn.mixture import GaussianMixture

        gmm = GaussianMixture(
            n_components=n_components,
            covariance_type=covariance_type,
            random_state=random_state,
        )
        labels = gmm.fit_predict(df.values)
        bic = gmm.bic(df.values)
        aic = gmm.aic(df.values)

        logger.info(
            f"GMM: {n_components} components, covariance={covariance_type}, "
            f"BIC={bic:.2f}, AIC={aic:.2f}"
        )
        return ClusterResult(
            labels=labels, method="gmm", metric=self.metric,
            n_clusters=n_components,
        )

    def elbow(
        self,
        df: pd.DataFrame,
        k_range: range = range(2, 11),
        random_state: int = 42,
    ) -> pd.DataFrame:
        """Compute inertia and silhouette for each k (for elbow analysis).

        Parameters
        ----------
        df : pd.DataFrame
            Scaled expression matrix.
        k_range : range
            Range of k values to evaluate.
        random_state : int
            Random seed.

        Returns
        -------
        pd.DataFrame
            Columns: k, inertia, silhouette.
        """
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score

        results = []
        for k in k_range:
            km = KMeans(n_clusters=k, random_state=random_state, n_init=10)
            labels = km.fit_predict(df.values)
            sil = silhouette_score(df.values, labels)
            results.append({"k": k, "inertia": km.inertia_, "silhouette": sil})
            logger.info(f"elbow k={k}: inertia={km.inertia_:.2f}, silhouette={sil:.4f}")
        return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Class 3: ClusterPlotter
# ---------------------------------------------------------------------------

class ClusterPlotter:
    """Visualize clustering results.

    All plot methods accept output path and produce publication-ready figures.
    """

    def heatmap(
        self,
        df_raw: pd.DataFrame,
        labels: np.ndarray,
        output: str | Path,
        cmap: str = "RdBu_r",
        figsize_per_gene: float = 0.04,
        figsize_per_sample: float = 0.5,
        vmin: float | None = None,
        vmax: float | None = None,
        max_genes: int = 2000,
    ) -> None:
        """Clustered heatmap with cluster color bar.

        Parameters
        ----------
        df_raw : pd.DataFrame
            Raw (unscaled) expression matrix (genes x samples).
        labels : np.ndarray
            Cluster label for each gene.
        output : str or Path
            Output file path (.png/.pdf).
        cmap : str
            Colormap.
        vmin, vmax : float, optional
            Color scale limits (auto from 2nd/98th percentile if None).
        max_genes : int
            Subsample for readability if more genes.
        """
        import seaborn as sns
        import matplotlib.pyplot as plt

        df = df_raw.copy()
        df["_cluster"] = labels
        df = df.sort_values("_cluster")
        cluster_labels = df.pop("_cluster").values

        if df.shape[0] > max_genes:
            step = df.shape[0] // max_genes
            df = df.iloc[::step]
            cluster_labels = cluster_labels[::step]
            logger.info(f"Heatmap subsampled to {df.shape[0]} genes")

        # z-score for display
        df_z = df.subtract(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0).fillna(0)
        if vmin is None:
            vmin = df_z.quantile(0.02).min()
        if vmax is None:
            vmax = df_z.quantile(0.98).max()

        n_genes, n_samples = df_z.shape
        fig_h = max(4, n_genes * figsize_per_gene)
        fig_w = max(6, n_samples * figsize_per_sample)

        fig, (ax_bar, ax_heat) = plt.subplots(
            1, 2, figsize=(fig_w + 1, fig_h),
            gridspec_kw={"width_ratios": [0.05, 1], "wspace": 0.01},
        )

        # cluster color bar
        unique = sorted(set(cluster_labels))
        palette = {c: plt.cm.tab20(i / max(len(unique), 1)) for i, c in enumerate(unique)}
        for i, cl in enumerate(cluster_labels):
            ax_bar.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=palette[cl], edgecolor="none"))
        ax_bar.set_xlim(0, 1)
        ax_bar.set_ylim(0, len(cluster_labels))
        ax_bar.invert_yaxis()
        ax_bar.set_xticks([])
        ax_bar.set_yticks([])
        ax_bar.set_ylabel("Cluster", fontsize=8)

        # heatmap
        sns.heatmap(
            df_z, ax=ax_heat, cmap=cmap, vmin=vmin, vmax=vmax,
            xticklabels=True, yticklabels=False, linewidths=0,
            cbar_kws={"shrink": 0.5, "label": "Z-score"},
        )
        ax_heat.set_xlabel("")
        ax_heat.set_ylabel("")
        ax_heat.tick_params(axis="x", labelsize=8, rotation=45)

        plt.savefig(output, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved heatmap to {output}")

    def umap(
        self,
        df_scaled: pd.DataFrame,
        labels: np.ndarray,
        output: str | Path,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        random_state: int = 42,
    ) -> None:
        """UMAP colored by cluster.

        Parameters
        ----------
        df_scaled : pd.DataFrame
            Scaled expression matrix.
        labels : np.ndarray
            Cluster labels.
        output : str or Path
            Output file path.
        n_neighbors : int
            UMAP n_neighbors.
        min_dist : float
            UMAP min_dist.
        random_state : int
            Random seed.

        Requires: pip install umap-learn
        """
        try:
            import umap
        except ImportError:
            raise ImportError("UMAP requires umap-learn: pip install umap-learn")

        import matplotlib.pyplot as plt

        reducer = umap.UMAP(
            n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state,
        )
        embedding = reducer.fit_transform(df_scaled.values)

        fig, ax = plt.subplots(figsize=(8, 6))
        scatter = ax.scatter(
            embedding[:, 0], embedding[:, 1],
            c=labels, cmap="tab20", s=5, alpha=0.7,
        )
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title("Gene clusters (UMAP)")
        plt.colorbar(scatter, ax=ax, label="Cluster")
        plt.savefig(output, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved UMAP to {output}")

    def dendrogram(
        self,
        Z: np.ndarray,
        output: str | Path,
        p: int = 30,
        figsize: tuple[int, int] = (12, 6),
    ) -> None:
        """Hierarchical clustering dendrogram.

        Parameters
        ----------
        Z : np.ndarray
            Linkage matrix (from ClusterResult.Z).
        output : str or Path
            Output file path.
        p : int
            Number of leaf nodes to show (truncation).
        figsize : tuple
            Figure size.
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)
        dendrogram(
            Z, truncate_mode="lastp", p=p,
            leaf_rotation=90, leaf_font_size=8, ax=ax,
        )
        ax.set_title("Hierarchical clustering dendrogram")
        ax.set_xlabel("Gene index (or cluster size)")
        ax.set_ylabel("Distance")
        plt.savefig(output, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved dendrogram to {output}")

    def elbow(
        self,
        elbow_df: pd.DataFrame,
        output: str | Path,
        figsize: tuple[int, int] = (10, 4),
    ) -> None:
        """Elbow (inertia) + silhouette analysis plot.

        Parameters
        ----------
        elbow_df : pd.DataFrame
            Output of DistanceCalculator.elbow().
        output : str or Path
            Output file path.
        figsize : tuple
            Figure size.
        """
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        ax1.plot(elbow_df["k"], elbow_df["inertia"], "bo-")
        ax1.set_xlabel("k")
        ax1.set_ylabel("Inertia")
        ax1.set_title("Elbow plot")

        ax2.plot(elbow_df["k"], elbow_df["silhouette"], "ro-")
        ax2.set_xlabel("k")
        ax2.set_ylabel("Silhouette score")
        ax2.set_title("Silhouette analysis")

        plt.tight_layout()
        plt.savefig(output, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved elbow plot to {output}")


# ---------------------------------------------------------------------------
# Pipeline orchestrator
# ---------------------------------------------------------------------------

def run_clustering(
    matrix: str | Path,
    output_dir: str | Path,
    method: Literal["hierarchical", "kmeans", "leiden", "dbscan", "spectral", "gmm"] = "hierarchical",
    n_clusters: int = 5,
    min_mean: float = 1.0,
    log_transform: bool = True,
    quantile_normalize: bool = False,
    zscore: bool = True,
    transpose: bool = False,
    linkage_method: str = "ward",
    metric: str = "euclidean",
    minkowski_p: float = 2.0,
    leiden_resolution: float = 1.0,
    random_state: int = 42,
    elbow: bool = False,
    umap: bool = False,
    heatmap: bool = True,
    dendrogram: bool = False,
    # DBSCAN parameters
    dbscan_eps: float = 0.5,
    dbscan_min_samples: int = 5,
    # Spectral parameters
    spectral_n_neighbors: int = 10,
    # GMM parameters
    gmm_covariance: str = "full",
) -> pd.DataFrame:
    """Run full clustering pipeline using the three classes.

    Returns
    -------
    pd.DataFrame
        Cluster assignments with columns: Gene, Cluster.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Preprocess
    prep = ExpressionPreprocessor(
        min_mean=min_mean, log_transform=log_transform,
        quantile_normalize=quantile_normalize, zscore=zscore,
    ).load(matrix, transpose=transpose)

    # 2. Distance + clustering
    calc = DistanceCalculator(metric=metric, minkowski_p=minkowski_p)

    if elbow:
        elbow_df = calc.elbow(prep.scaled, random_state=random_state)
        elbow_df.to_csv(output_dir / "elbow.tsv", sep="\t", index=False)
        ClusterPlotter().elbow(elbow_df, output_dir / "elbow.png")

    if method == "hierarchical":
        result = calc.hierarchical(prep.scaled, n_clusters=n_clusters, method=linkage_method)
    elif method == "kmeans":
        result = calc.kmeans(prep.scaled, n_clusters=n_clusters, random_state=random_state)
    elif method == "leiden":
        result = calc.leiden(prep.scaled, resolution=leiden_resolution, random_state=random_state)
    elif method == "dbscan":
        result = calc.dbscan(prep.scaled, eps=dbscan_eps, min_samples=dbscan_min_samples)
    elif method == "spectral":
        result = calc.spectral(prep.scaled, n_clusters=n_clusters, n_neighbors=spectral_n_neighbors, random_state=random_state)
    elif method == "gmm":
        result = calc.gmm(prep.scaled, n_components=n_clusters, covariance_type=gmm_covariance, random_state=random_state)
    else:
        raise ValueError(f"Unknown method: {method}")

    # save assignments
    assignments = result.to_dataframe(prep.gene_names)
    assignments.to_csv(output_dir / "cluster_assignments.tsv", sep="\t", index=False)
    summary = result.summary()
    summary.to_csv(output_dir / "cluster_summary.tsv", sep="\t", index=False)
    logger.info(f"Cluster sizes:\n{summary.to_string(index=False)}")

    # 3. Visualization
    plotter = ClusterPlotter()
    if heatmap:
        plotter.heatmap(prep.raw, result.labels, output_dir / "heatmap.png")
    if umap:
        plotter.umap(prep.scaled, result.labels, output_dir / "umap.png")
    if dendrogram and result.Z is not None:
        plotter.dendrogram(result.Z, output_dir / "dendrogram.png")

    return assignments


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cluster gene expression matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Hierarchical, cosine distance, complete linkage
  python expression_cluster.py matrix.tsv -o out/ -m hierarchical --metric cosine --linkage complete

  # K-means with elbow plot
  python expression_cluster.py matrix.tsv -o out/ -m kmeans --elbow

  # Manhattan distance + average linkage
  python expression_cluster.py matrix.tsv -o out/ --metric manhattan --linkage average
""",
    )
    parser.add_argument("matrix", help="Expression matrix (TSV, genes x samples)")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "-m", "--method",
        choices=["hierarchical", "kmeans", "leiden", "dbscan", "spectral", "gmm"],
        default="hierarchical",
        help="Clustering method (default: hierarchical)",
    )
    parser.add_argument("-k", "--n-clusters", type=int, default=5, help="Number of clusters")
    parser.add_argument("--min-mean", type=float, default=1.0, help="Min mean expression filter")
    parser.add_argument("--no-log", action="store_true", help="Skip log2(x+1) transform")
    parser.add_argument("--no-zscore", action="store_true", help="Skip z-score scaling")
    parser.add_argument("--qn", action="store_true", help="Quantile normalization (remove batch effects)")
    parser.add_argument("--transpose", action="store_true", help="Transpose (rows=samples)")
    parser.add_argument(
        "--linkage", default="ward",
        choices=["ward", "complete", "average", "single"],
        help="Linkage method for hierarchical (default: ward)",
    )
    parser.add_argument(
        "--metric", default="euclidean",
        choices=list(DistanceCalculator.SUPPORTED_METRICS),
        help="Distance metric (default: euclidean)",
    )
    parser.add_argument("--minkowski-p", type=float, default=2.0, help="Minkowski p")
    parser.add_argument("--leiden-resolution", type=float, default=1.0, help="Leiden resolution")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed")
    parser.add_argument("--elbow", action="store_true", help="Generate elbow plot")
    parser.add_argument("--umap", action="store_true", help="Generate UMAP plot")
    parser.add_argument("--no-heatmap", action="store_true", help="Skip heatmap")
    parser.add_argument("--dendrogram", action="store_true", help="Dendrogram (hierarchical only)")
    # DBSCAN parameters
    parser.add_argument("--dbscan-eps", type=float, default=0.5, help="DBSCAN eps (neighborhood size)")
    parser.add_argument("--dbscan-min-samples", type=int, default=5, help="DBSCAN min_samples")
    # Spectral parameters
    parser.add_argument("--spectral-n-neighbors", type=int, default=10, help="Spectral n_neighbors")
    # GMM parameters
    parser.add_argument("--gmm-covariance", default="full", choices=["full", "tied", "diag", "spherical"],
                        help="GMM covariance type")

    args = parser.parse_args(argv)
    if args.dendrogram and args.method != "hierarchical":
        parser.error("--dendrogram only works with hierarchical method")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    run_clustering(
        matrix=args.matrix,
        output_dir=args.output_dir,
        method=args.method,
        n_clusters=args.n_clusters,
        min_mean=args.min_mean,
        log_transform=not args.no_log,
        quantile_normalize=args.qn,
        zscore=not args.no_zscore,
        transpose=args.transpose,
        linkage_method=args.linkage,
        metric=args.metric,
        minkowski_p=args.minkowski_p,
        leiden_resolution=args.leiden_resolution,
        random_state=args.random_state,
        elbow=args.elbow,
        umap=args.umap,
        heatmap=not args.no_heatmap,
        dendrogram=args.dendrogram,
        dbscan_eps=args.dbscan_eps,
        dbscan_min_samples=args.dbscan_min_samples,
        spectral_n_neighbors=args.spectral_n_neighbors,
        gmm_covariance=args.gmm_covariance,
    )


if __name__ == "__main__":
    main()
