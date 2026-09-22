#!/usr/bin/env python3
"""Pseudo-bulk DEG pipeline: prepare → PyDESeq2 → plot.

Single entry point that:
  1. Aggregates single-cell counts into pseudo-bulk matrices
  2. Runs PyDESeq2 for pairwise DEG (replaces edgeR R subprocess)
  3. Generates PCA exploration, sample filtering, volcano, heatmap,
     GO dotplot, and overlap figures

Usage (CLI-only, no config file needed):
  python pseudobulk.py --out-dir ./results \\
    --h5ad ovaries=/path/to/ov.h5ad \\
    --h5ad uterus=/path/to/ut.h5ad \\
    --compare ovaries:luanchao-21310:luanchao-11238=Aging_vs_Youth \\
    --compare uterus:zigong-21310:ZIGONG-21224=Sham_vs_OV-POI \\
    --compare uterus:zigong-11238:zigong-1411200=Intact_vs_OV-only

--h5ad   tissue=/path/to.h5ad         (repeatable)
--compare tissue:ref:treat[=label]    (repeatable, label optional)
"""
import argparse, os, sys, warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

# ────────────────────────────────────────────────────────────
# 1. Config
# ────────────────────────────────────────────────────────────

def load_config(path):
    """Load a YAML config file into the internal config dict.

    Reads a YAML file with tissue-level comparison definitions and
    enriches each comparison entry with a ``key`` (``ref_vs_treat``)
    and a fallback ``label`` (defaults to the key if not specified).

    Args:
        path: Path to a YAML file structured as::

            <tissue>:
              h5ad: /path/to.h5ad
              comparisons:
                - ref:  <reference_sample_id>
                  treat: <treatment_sample_id>
                  label: <optional_short_name>

    Returns:
        dict: ``{tissue: {h5ad: str, comparisons: [{ref, treat, label, key}]}}``
    """
    import yaml
    with open(path) as f:
        raw = yaml.safe_load(f)
    for tissue, cfg in raw.items():
        for i, comp in enumerate(cfg["comparisons"]):
            comp["key"] = f"{comp['ref']}_vs_{comp['treat']}"
            if "label" not in comp:
                comp["label"] = comp["key"]
    return raw


def comp_map(comparisons):
    """Build a mapping from comparison key to human-readable label.

    Args:
        comparisons: List of comparison dicts with ``key`` and ``label``.

    Returns:
        dict: ``{ref_vs_treat: short_label, ...}``
    """
    return {c["key"]: c["label"] for c in comparisons}


# ────────────────────────────────────────────────────────────
# 2. Pseudobulk preparation
# ────────────────────────────────────────────────────────────

def _split_pseudo_reps(adata_sub, sample_id, raw_sub, var_names,
                       n_pseudo_reps=3, min_cells_per_rep=10, seed=42):
    """Split cells of one sample into pseudo-replicates for DESeq2.

    Randomly assigns cells to *n_pseudo_reps* groups and sums counts
    per group. Used when only 1 biological replicate exists per
    condition, giving DESeq2 enough ``replicates'' to estimate
    dispersion.

    Note: pseudo-replicates are NOT independent biological replicates.
    Dispersion will be underestimated; p-values are anti-conservative.

    Args:
        adata_sub: Subset AnnData for one sample (cells × genes).
        sample_id: Original sample ID (used as prefix for rep names).
        raw_sub: Raw count AnnData for the same subset.
        var_names: Gene names (for building the count dict).
        n_pseudo_reps: Number of pseudo-replicates to create.
        min_cells_per_rep: Minimum cells per pseudo-replicate.
            If the sample has fewer than ``n_pseudo_reps * min_cells_per_rep``
            cells, the number of replicates is reduced.
        seed: Random seed for reproducibility.

    Returns:
        list[dict]: One dict per pseudo-replicate with keys
        ``sample``, ``cell_type``, ``n_cells``, ``n_counts``, plus gene counts.
    """
    n_cells_total = adata_sub.shape[0]
    # Adjust n_pseudo_reps if too few cells
    actual_reps = min(n_pseudo_reps, n_cells_total // min_cells_per_rep)
    if actual_reps < 2:
        actual_reps = min(2, n_cells_total)  # at least try 2

    if actual_reps < 2 or n_cells_total < 4:
        # Can't split: return single aggregate
        counts = np.array(raw_sub.X.sum(axis=0)).flatten()
        return [{
            "sample": sample_id,
            "cell_type": adata_sub.obs["cell_type"].iloc[0] if "cell_type" in adata_sub.obs.columns else "unknown",
            "n_cells": n_cells_total,
            "n_counts": float(counts.sum()),
            **{gene: float(c) for gene, c in zip(var_names, counts)},
        }]

    # Equal distribution: each rep gets floor(N/M) or ceil(N/M) cells
    # then shuffle for randomness
    assignments = np.array([i % actual_reps for i in range(n_cells_total)])
    rng = np.random.RandomState(seed)
    rng.shuffle(assignments)

    rows = []
    cell_type = adata_sub.obs["cell_type"].iloc[0] if "cell_type" in adata_sub.obs.columns else "unknown"
    for rep_idx in range(actual_reps):
        rep_mask = assignments == rep_idx
        n_cells_rep = int(rep_mask.sum())
        if n_cells_rep == 0:
            continue
        counts = np.array(raw_sub[rep_mask].X.sum(axis=0)).flatten()
        rows.append({
            "sample": f"{sample_id}_rep{rep_idx}",
            "cell_type": cell_type,
            "n_cells": n_cells_rep,
            "n_counts": float(counts.sum()),
            **{gene: float(c) for gene, c in zip(var_names, counts)},
        })
    return rows


def aggregate_pseudobulk(adata, cell_type, sample_col="sample_id",
                         type_col="cell_type", min_cells=30,
                         min_counts=1000, n_pseudo_reps=0):
    """Sum raw UMI counts per sample for one cell type.

    Filters out samples with fewer than *min_cells* cells or
    *min_counts* total UMIs, then aggregates (sums) the raw count
    matrix across all cells belonging to each remaining sample.

    When *n_pseudo_reps* > 0 and only 1 sample passes the filter,
    that sample's cells are randomly split into *n_pseudo_reps*
    pseudo-replicates so PyDESeq2 can estimate dispersion.

    Args:
        adata: AnnData object with ``.raw`` holding the full count matrix.
        cell_type: Cell type to subset (must exist in ``adata.obs[type_col]``).
        sample_col: Column in ``adata.obs`` identifying biological samples.
        type_col: Column in ``adata.obs`` identifying cell type annotations.
        min_cells: Minimum number of cells per sample to retain it.
        min_counts: Minimum total UMI counts per sample to retain it.
        n_pseudo_reps: If > 0, split single-sample conditions into this
            many pseudo-replicates. 0 = disabled (default).

    Returns:
        tuple:
            - **df** (*pd.DataFrame | None*): Rows = samples (or pseudo-reps),
              columns = [cell_type, n_cells, n_counts, gene1, gene2, ...].
              Index is sample ID. ``None`` if fewer than 2 rows.
            - **samples_kept** (*list[str]*): Sample IDs that passed the filter.
    """
    sub = adata[adata.obs[type_col] == cell_type].copy()
    size_by_sample = sub.obs.groupby(sample_col).size()
    samples_by_cells = set(size_by_sample[size_by_sample >= min_cells].index.tolist())

    raw_sub = sub.raw.to_adata()

    # Compute total counts per sample
    counts_by_sample = {}
    for sample in samples_by_cells:
        mask = sub.obs[sample_col] == sample
        total = float(np.array(raw_sub[mask].X.sum()).flatten()[0])
        counts_by_sample[sample] = total

    # Apply min_counts filter
    samples_kept = [s for s in samples_by_cells
                    if counts_by_sample.get(s, 0) >= min_counts]

    if len(samples_kept) < 1:
        return None, []

    # Split into pseudo-replicates when n_pseudo_reps > 0
    # (needed when any comparison will have only 1 sample per condition)
    use_pseudo_reps = (n_pseudo_reps > 0)

    rows = []
    for i, sample in enumerate(samples_kept):
        mask = sub.obs[sample_col] == sample
        sub_cells = sub[mask]
        raw_cells = raw_sub[mask]

        if use_pseudo_reps:
            rep_rows = _split_pseudo_reps(
                sub_cells, sample, raw_cells, raw_sub.var_names,
                n_pseudo_reps=n_pseudo_reps, seed=42 + i,
            )
            rows.extend(rep_rows)
        else:
            counts = np.array(raw_cells.X.sum(axis=0)).flatten()
            rows.append({
                "sample": sample,
                "cell_type": cell_type,
                "n_cells": int(mask.sum()),
                "n_counts": counts_by_sample[sample],
                **{gene: float(c) for gene, c in zip(raw_sub.var_names, counts)},
            })

    if len(rows) < 2:
        return None, []

    df = pd.DataFrame(rows).set_index("sample")
    df.index.name = "sample"
    return df, samples_kept


def prepare_counts(h5ad_path, out_dir, label, min_cells=30, min_counts=1000,
                   n_pseudo_reps=0):
    """Prepare pseudo-bulk count CSVs for all cell types in an h5ad.

    Iterates over every unique cell type in the h5ad, aggregates counts
    per sample, and writes one CSV per cell type to *out_dir*.

    Args:
        h5ad_path: Path to an annotated h5ad file.
        out_dir: Directory for output count CSVs.
        label: Tissue label used as filename prefix (e.g. ``"uterus"``).
        min_cells: Minimum cells per sample to retain.
        min_counts: Minimum total UMIs per sample to retain.
        n_pseudo_reps: If > 0, split single-sample conditions into
            pseudo-replicates (see :func:`aggregate_pseudobulk`).

    Returns:
        tuple:
            - **written** (*list[dict]*): Metadata for each written file.
            - **adata** (*AnnData*): The loaded AnnData object (for downstream PCA, etc.).
    """
    import scanpy as sc
    os.makedirs(out_dir, exist_ok=True)
    print(f"[{label}] Loading {h5ad_path} ...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  shape: {adata.shape}, samples: {adata.obs['sample_id'].nunique()}")

    targets = sorted(adata.obs["cell_type"].unique().tolist())
    written = []
    for ct in targets:
        df, kept = aggregate_pseudobulk(adata, ct, min_cells=min_cells,
                                        min_counts=min_counts,
                                        n_pseudo_reps=n_pseudo_reps)
        if df is None:
            continue
        ct_safe = ct.replace(" ", "_").replace("/", "-")
        out_path = os.path.join(out_dir, f"{label}_{ct_safe}_counts.csv")
        df.to_csv(out_path)
        n_rows = len(df)
        pseudo_note = ""
        if n_pseudo_reps > 0 and n_rows > len(kept):
            pseudo_note = f" ({len(kept)} samples → {n_rows} pseudo-reps)"
        written.append({
            "label": label, "cell_type": ct,
            "n_samples": len(kept), "n_rows": n_rows,
            "n_genes": df.shape[1] - 3,
            "path": out_path,
        })
        print(f"  [{ct}] n_samples={len(kept)}, n_rows={n_rows}{pseudo_note} -> {out_path}")
    return written, adata


# ────────────────────────────────────────────────────────────
# 3. PyDESeq2 DEG analysis
# ────────────────────────────────────────────────────────────

def run_pydeseq2(counts_df, ref_id, treat_id):
    """Run PyDESeq2 on a pseudobulk count matrix for one comparison.

    Takes a pseudobulk DataFrame (rows=samples, columns=metadata+genes),
    subsets to the two groups, and runs DESeq2 analysis.

    Args:
        counts_df: DataFrame with ``cell_type``, ``n_cells``, ``n_counts``
            columns plus gene columns. Index = sample IDs.
        ref_id: Sample ID of the reference group.
        treat_id: Sample ID of the treatment group.

    Returns:
        pd.DataFrame: DEG results with columns ``gene``, ``cell_type``,
        ``comparison``, ``baseMean``, ``log2FC``, ``lfcSE``, ``stat``,
        ``pval``, ``pval_adj``, ``significant``.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("  ERROR: pydeseq2 not installed. pip install pydeseq2", file=sys.stderr)
        return pd.DataFrame()

    # Separate metadata from gene columns
    metadata_cols = ["cell_type", "n_cells", "n_counts"]
    gene_cols = [c for c in counts_df.columns if c not in metadata_cols]

    if len(gene_cols) < 100:
        return pd.DataFrame()

    # Subset to the two groups (handle pseudo-replicate IDs like sample_rep0)
    ref_samples = [s for s in counts_df.index if s == ref_id or s.startswith(f"{ref_id}_rep")]
    treat_samples = [s for s in counts_df.index if s == treat_id or s.startswith(f"{treat_id}_rep")]
    selected = ref_samples + treat_samples
    if len(selected) < 2 or len(ref_samples) < 1 or len(treat_samples) < 1:
        print(f"  [skip] not enough samples for {ref_id} vs {treat_id}")
        return pd.DataFrame()

    sub = counts_df.loc[selected]
    counts = sub[gene_cols].round().astype(int).copy()

    # Filter genes with zero total counts
    gene_totals = counts.sum(axis=0)
    counts = counts.loc[:, gene_totals > 0]
    if counts.shape[1] < 10:
        return pd.DataFrame()

    # Metadata: assign condition based on prefix matching
    metadata = pd.DataFrame({
        "condition": [ref_id if s in ref_samples else treat_id for s in counts.index]
    }, index=counts.index)
    metadata["condition"] = pd.Categorical(metadata["condition"],
                                           categories=[ref_id, treat_id])

    cell_type = counts_df["cell_type"].iloc[0] if "cell_type" in counts_df.columns else "unknown"

    # Fit DESeq2 (no fallback — pseudo-reps should provide enough replicates)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~condition",
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", treat_id, ref_id])
    stat_res.summary()
    results_df = stat_res.results_df.copy()

    # Format output
    results_df["gene"] = results_df.index
    results_df = results_df.rename(columns={
        "log2FoldChange": "log2FC",
        "pvalue": "pval",
        "padj": "pval_adj",
    })
    results_df["cell_type"] = cell_type
    results_df["comparison"] = f"{ref_id}_vs_{treat_id}"
    results_df["significant"] = (
        (results_df["pval_adj"] < 0.05) & (results_df["log2FC"].abs() > 1)
    )
    # Fill NaN significance as False
    results_df["significant"] = results_df["significant"].fillna(False)

    cols = ["gene", "cell_type", "comparison", "baseMean", "log2FC",
            "lfcSE", "stat", "pval", "pval_adj", "significant"]
    results_df = results_df[[c for c in cols if c in results_df.columns]]

    n_sig = int(results_df["significant"].sum())
    print(f"    {cell_type}: {ref_id} vs {treat_id} — "
          f"{len(results_df)} genes, {n_sig} significant")

    return results_df


# ────────────────────────────────────────────────────────────
# 4. Diagnostic plots (PyDESeq2)
# ────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False


def pydeseq2_diagnostic_plots(results_df, output_prefix, treat_id, ref_id):
    """Generate p-value histogram and MA plot for PyDESeq2 results.

    Args:
        results_df: DataFrame from :func:`run_pydeseq2` with ``log2FC``,
            ``pval``, ``baseMean`` columns.
        output_prefix: Path prefix for output PNGs (appended with suffix).
        treat_id: Treatment sample ID (for plot titles).
        ref_id: Reference sample ID (for plot titles).

    Returns:
        bool: True if plots were generated, False if too few genes.
    """
    if len(results_df) < 10:
        return False

    valid = results_df.dropna(subset=["pval"])
    if len(valid) < 10:
        return False

    # ── P-value histogram ──
    fig, ax = plt.subplots(figsize=(8, 6))
    pvals = valid["pval"].values
    n_bins = 50
    ax.hist(pvals, bins=n_bins, color="steelblue", edgecolor="white",
            alpha=0.85, label="Observed")
    expected = len(pvals) / n_bins
    ax.axhline(expected, color="red", linestyle="--", linewidth=1.5,
               label=f"Expected uniform ({expected:.0f}/bin)")

    # Enrichment near 0 annotation
    near_zero = np.sum(pvals < 0.05)
    total = len(pvals)
    ax.axvline(0.05, color="orange", linestyle=":", linewidth=1, alpha=0.7)
    ax.text(0.07, ax.get_ylim()[1] * 0.95,
            f"p < 0.05: {near_zero}/{total} ({100*near_zero/total:.1f}%)",
            fontsize=9, color="orange", va="top")

    ax.set_xlabel("P-value", fontsize=11)
    ax.set_ylabel("Frequency", fontsize=11)
    ax.set_title(f"P-value distribution\n{ref_id} vs {treat_id}",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_pval_hist.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close()

    # ── MA plot (log2FC vs baseMean) ──
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_df = valid.copy()
    plot_df["baseMean_log"] = np.log10(plot_df["baseMean"].clip(lower=1))
    not_sig = plot_df[~plot_df.get("significant", False)]
    sig = plot_df[plot_df.get("significant", False)]

    ax.scatter(not_sig["baseMean_log"], not_sig["log2FC"],
               s=2, c="#cccccc", alpha=0.3, edgecolor="none",
               label=f"NS ({len(not_sig)})", zorder=1)
    if len(sig) > 0:
        up = sig[sig["log2FC"] > 0]
        dn = sig[sig["log2FC"] < 0]
        if len(up) > 0:
            ax.scatter(up["baseMean_log"], up["log2FC"],
                       s=25, c="#d62728", alpha=0.9, edgecolor="darkred",
                       linewidth=0.3, label=f"Up ({len(up)})", zorder=4)
        if len(dn) > 0:
            ax.scatter(dn["baseMean_log"], dn["log2FC"],
                       s=25, c="#1f77b4", alpha=0.9, edgecolor="darkblue",
                       linewidth=0.3, label=f"Down ({len(dn)})", zorder=4)

    ax.axhline(0, color="black", linestyle="-", linewidth=0.5, alpha=0.3)
    ax.axhline(1, color="black", linestyle="--", linewidth=0.7, alpha=0.6)
    ax.axhline(-1, color="black", linestyle="--", linewidth=0.7, alpha=0.6)
    ax.set_xlabel("log10(baseMean)", fontsize=11)
    ax.set_ylabel("log2FC", fontsize=11)
    ax.set_title(f"MA plot: {ref_id} vs {treat_id}",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10, loc="upper right", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_ma_plot.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close()
    return True


# ────────────────────────────────────────────────────────────
# 5. QC / Exploration plots (best practices)
# ────────────────────────────────────────────────────────────

def pca_exploration_plot(adata, output_path, sample_col="sample_id",
                         condition_map=None):
    """PCA on pseudobulk samples to identify sources of variation.

    Aggregates counts per sample, normalizes, log-transforms, scales,
    computes PCA, and plots PC1 vs PC2 colored by condition and
    cell_type. Following sc-best-practices: explore before modeling.

    Args:
        adata: AnnData with raw counts in ``.raw``.
        output_path: Path to save the PNG figure.
        sample_col: Column identifying biological samples.
        condition_map: Optional dict ``{sample_id: condition}`` for coloring.

    Returns:
        bool: True if PCA was computed, False if too few samples.
    """
    import scanpy as sc

    samples = adata.obs[sample_col].unique().tolist()
    if len(samples) < 3:
        print(f"  PCA skipped: only {len(samples)} samples (need >= 3)")
        return False

    # Aggregate per sample
    raw = adata.raw.to_adata()
    rows = []
    for sample in samples:
        mask = adata.obs[sample_col] == sample
        counts = np.array(raw[mask].X.sum(axis=0)).flatten()
        cell_types = adata.obs.loc[mask, "cell_type"].value_counts()
        ct_label = cell_types.index[0] if len(cell_types) == 1 else "mixed"
        row = {
            "sample": sample,
            "cell_type": ct_label,
            "n_cells": int(mask.sum()),
            "n_counts": float(counts.sum()),
        }
        if condition_map:
            row["condition"] = condition_map.get(sample, "unknown")
        rows.append(row)

    meta = pd.DataFrame(rows).set_index("sample")
    count_rows = []
    for sample in samples:
        mask = adata.obs[sample_col] == sample
        count_rows.append(np.array(raw[mask].X.sum(axis=0)).flatten())
    counts_df = pd.DataFrame(count_rows, index=meta.index, columns=raw.var_names)

    # Build AnnData for PCA
    pb_adata = sc.AnnData(X=counts_df.values, obs=meta,
                          var=pd.DataFrame(index=counts_df.columns))
    sc.pp.normalize_total(pb_adata, target_sum=1e6)
    sc.pp.log1p(pb_adata)
    sc.pp.scale(pb_adata, max_value=10)
    sc.tl.pca(pb_adata)

    # Plot
    color_by = ["cell_type", "n_cells"]
    if condition_map:
        color_by.append("condition")

    n_panels = len(color_by)
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    pca_coords = pb_adata.obsm["X_pca"][:, :2]

    for i, col in enumerate(color_by):
        ax = axes[i]
        vals = pb_adata.obs[col]

        if pd.api.types.is_numeric_dtype(vals):
            # Continuous coloring (e.g., n_cells)
            sc = ax.scatter(pca_coords[:, 0], pca_coords[:, 1],
                            c=vals.values, cmap="viridis", s=100,
                            edgecolor="black", linewidth=0.5)
            plt.colorbar(sc, ax=ax, label=col)
        else:
            # Categorical coloring
            categories = sorted(vals.unique())
            palette = sns.color_palette("colorblind", len(categories))
            for j, cat in enumerate(categories):
                mask = (vals == cat).values
                ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1],
                           c=[palette[j]], label=cat, s=100,
                           edgecolor="black", linewidth=0.5)
            ax.legend(fontsize=8, loc="best")

        ax.set_xlabel("PC1", fontsize=10)
        ax.set_ylabel("PC2", fontsize=10)
        ax.set_title(f"PCA colored by {col}", fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Label sample names
        for k, sample_name in enumerate(meta.index):
            ax.annotate(sample_name, (pca_coords[k, 0], pca_coords[k, 1]),
                        fontsize=6, alpha=0.7, ha="left", va="bottom")

    plt.suptitle("Pseudobulk PCA exploration", fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  PCA plot: {output_path}")
    return True


def sample_filter_plot(adata, output_path, sample_col="sample_id",
                       type_col="cell_type", min_cells=10, min_counts=1000):
    """Visualize pseudobulk sample quality: n_cells vs n_counts.

    Similar to decoupler.pl.filter_samples. Each dot is one
    pseudobulk sample (cell type x donor). Dashed lines show filtering
    thresholds. Samples in the upper-right quadrant pass QC.

    Args:
        adata: AnnData with raw counts in ``.raw``.
        output_path: Path to save the PNG figure.
        sample_col: Column identifying biological samples.
        type_col: Column identifying cell type annotations.
        min_cells: Minimum cells threshold (vertical line).
        min_counts: Minimum total counts threshold (horizontal line).

    Returns:
        pd.DataFrame: Summary table with sample, cell_type, n_cells, n_counts, pass_qc.
    """
    raw = adata.raw.to_adata()
    records = []
    for (sample, ct), group_idx in adata.obs.groupby([sample_col, type_col]).groups.items():
        mask = adata.obs.index.isin(group_idx)
        n_cells = int(mask.sum())
        n_counts = float(np.array(raw[mask].X.sum()).flatten()[0])
        records.append({
            "sample": sample, "cell_type": ct,
            "n_cells": n_cells, "n_counts": n_counts,
            "pass_qc": n_cells >= min_cells and n_counts >= min_counts,
        })

    df = pd.DataFrame(records)
    if len(df) == 0:
        return df

    fig, ax = plt.subplots(figsize=(8, 6))

    passed = df[df["pass_qc"]]
    failed = df[~df["pass_qc"]]

    if len(failed) > 0:
        ax.scatter(failed["n_cells"], failed["n_counts"],
                   c="#cccccc", s=40, alpha=0.6, edgecolor="black",
                   linewidth=0.5, label=f"Fail ({len(failed)})", zorder=1)
    if len(passed) > 0:
        ax.scatter(passed["n_cells"], passed["n_counts"],
                   c="#2ca02c", s=40, alpha=0.8, edgecolor="black",
                   linewidth=0.5, label=f"Pass ({len(passed)})", zorder=2)

    ax.axvline(min_cells, color="red", linestyle="--", linewidth=1, alpha=0.7)
    ax.axhline(min_counts, color="red", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of cells (log)", fontsize=11)
    ax.set_ylabel("Total counts (log)", fontsize=11)
    ax.set_title(f"Pseudobulk sample filtering\n"
                 f"(min_cells={min_cells}, min_counts={min_counts})",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10, loc="lower right")

    # Annotate failed samples
    for _, row in failed.iterrows():
        ax.annotate(f"{row['sample']}\n({row['cell_type']})",
                    (row["n_cells"], row["n_counts"]),
                    fontsize=6, alpha=0.6, ha="left", va="bottom")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Filter plot: {output_path}")
    return df


def pseudobulk_summary_plot(written_records, output_path):
    """Bar chart of pseudobulk samples retained per cell type.

    Args:
        written_records: List of dicts from :func:`prepare_counts` with
            ``cell_type``, ``n_samples`` keys.
        output_path: Path to save the PNG figure.
    """
    if not written_records:
        return
    df = pd.DataFrame(written_records).sort_values("n_samples", ascending=True)
    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.35)))
    bars = ax.barh(df["cell_type"], df["n_samples"], color="#4c78a8",
                   edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, df["n_samples"]):
        ax.text(val + 0.2, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=9)
    ax.set_xlabel("Number of pseudobulk samples", fontsize=11)
    ax.set_title("Pseudobulk samples per cell type (after QC)",
                 fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Summary plot: {output_path}")


# ────────────────────────────────────────────────────────────
# 6. Result plotting (volcano, heatmap, GO, overlap)
# ────────────────────────────────────────────────────────────

def load_all_deg(deg_dir, label_prefix):
    """Load and concatenate all DEG CSVs matching a tissue prefix.

    Recursively scans *deg_dir* and its subdirectories for files
    named ``{label_prefix}*_DEG.csv`` and concatenates them.

    Args:
        deg_dir: Directory containing DEG CSV files (possibly in subdirectories).
        label_prefix: Filename prefix to filter by (e.g. ``"uterus_"``).

    Returns:
        pd.DataFrame: Combined DEG results, or empty DataFrame if none found.
    """
    rows = []
    if not os.path.isdir(deg_dir):
        return pd.DataFrame()
    for root, dirs, files in os.walk(deg_dir):
        for f in files:
            if f.startswith(label_prefix) and f.endswith("_DEG.csv"):
                df = pd.read_csv(os.path.join(root, f))
                rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def deg_count_heatmap(deg_df, cmap, title, output_path):
    """Heatmap of significant DEG counts per cell type × comparison.

    Pivots the DEG DataFrame to a cell-type × comparison matrix where
    each cell is the count of significant DEGs, and renders a heatmap.

    Args:
        deg_df: DataFrame with columns ``cell_type``, ``comparison``,
            ``significant``.
        cmap: ``{comparison_key: short_label}`` mapping for axis labels.
        title: Plot title.
        output_path: Path to save the PNG figure.
    """
    if len(deg_df) == 0:
        return
    deg_df = deg_df.copy()
    deg_df["comp_short"] = deg_df["comparison"].map(
        lambda x: cmap.get(x, x.replace("_vs_", " vs "))
    )
    pv = deg_df.pivot_table(
        index="cell_type", columns="comp_short", values="significant",
        aggfunc="sum", fill_value=0,
    )
    fig, ax = plt.subplots(figsize=(max(8, len(pv.columns) * 3), max(6, len(pv) * 0.4)))
    sns.heatmap(pv, annot=True, fmt=".0f", cmap="YlOrRd", ax=ax,
                cbar_kws={"label": "Significant DEG count"})
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha="right", fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def volcano_plot(deg_df, cell_type, comparison, output_path, title=None,
                 nlogp_clip=18, highlight_genes=None):
    """Volcano plot for one cell type × comparison.

    Plots log2FC vs -log10(pvalue) with up/down-regulated genes colored.
    All NS genes are shown as background (no filtering). Top 5 genes
    per direction by p-value are labeled with auto-repulsion (adjustText).
    Additional genes can be highlighted via *highlight_genes*. No log2FC
    clipping. Skips if no significant DEGs.

    Args:
        deg_df: DataFrame with ``gene``, ``cell_type``, ``comparison``,
            ``log2FC``, ``pval``, ``significant`` columns.
        cell_type: Cell type to plot.
        comparison: Comparison key (e.g. ``"ref_vs_treat"``).
        output_path: Path to save the PNG figure.
        title: Optional plot title (auto-generated if None).
        nlogp_clip: Clip -log10(p) display range to this value.
        highlight_genes: List of gene names to label on the plot regardless
            of significance. Labeled in red if up, blue if down, gray if NS.

    Returns:
        bool: True if the plot was generated, False if skipped.
    """
    from adjustText import adjust_text

    sub = deg_df[(deg_df["cell_type"] == cell_type) &
                 (deg_df["comparison"] == comparison)].copy()
    if len(sub) < 3:
        return False

    sub["nlogp"] = -np.log10(sub["pval"].clip(1e-300)).clip(upper=nlogp_clip)
    sub["sig"] = sub["significant"]

    sig_df = sub[sub["sig"]]
    up = sig_df[sig_df["log2FC"] > 0]
    dn = sig_df[sig_df["log2FC"] < 0]

    # Skip if no significant DEGs
    if len(sig_df) == 0:
        return False

    fig, ax = plt.subplots(figsize=(8, 6))
    not_sig = sub[~sub["sig"]]

    # Plot ALL NS genes as background
    if len(not_sig) > 0:
        ax.scatter(not_sig["log2FC"], not_sig["nlogp"],
                   s=2, c="#cccccc", alpha=0.3, edgecolor="none",
                   label=f"NS ({len(not_sig)})", zorder=1)
    if len(up) > 0:
        ax.scatter(up["log2FC"], up["nlogp"],
                   s=25, c="#d62728", alpha=0.9, edgecolor="darkred",
                   linewidth=0.3, label=f"Up ({len(up)})", zorder=4)
    if len(dn) > 0:
        ax.scatter(dn["log2FC"], dn["nlogp"],
                   s=25, c="#1f77b4", alpha=0.9, edgecolor="darkblue",
                   linewidth=0.3, label=f"Down ({len(dn)})", zorder=4)

    # Label top 5 genes per direction with auto-repulsion
    texts = []
    for _, r in pd.concat([up.nlargest(5, "nlogp"),
                           dn.nlargest(5, "nlogp")]).iterrows():
        texts.append(ax.text(r["log2FC"], r["nlogp"], r["gene"],
                             fontsize=7, color="black", zorder=5))

    # Highlight custom genes (e.g. aging markers, SASP factors)
    if highlight_genes:
        for gene_name in highlight_genes:
            gene_row = sub[sub["gene"] == gene_name]
            if len(gene_row) == 0:
                continue
            g = gene_row.iloc[0]
            color = "#d62728" if g["log2FC"] > 0 else "#1f77b4"
            ax.scatter(g["log2FC"], g["nlogp"], s=60, facecolors="none",
                       edgecolors=color, linewidths=1.5, zorder=6)
            texts.append(ax.text(g["log2FC"], g["nlogp"], gene_name,
                                 fontsize=8, fontweight="bold",
                                 color=color, zorder=7))

    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-",
                    color="gray", lw=0.5, alpha=0.6))

    # Threshold lines: p=0.05 and |log2FC|=1
    ax.axhline(-np.log10(0.05), color="black", linestyle="--",
               linewidth=0.7, alpha=0.6)
    ax.axvline(1, color="black", linestyle="--", linewidth=0.7, alpha=0.6)
    ax.axvline(-1, color="black", linestyle="--", linewidth=0.7, alpha=0.6)

    ax.set_xlabel("log2FC", fontsize=11)
    ax.set_ylabel("-log10(pvalue)", fontsize=11)
    if title is None:
        title = f'{cell_type}: {comparison.replace("_vs_", " vs ")}'
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend(fontsize=10, loc="upper right", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    return True


def top_deg_heatmap(adata, deg_df, cell_type, comparison, output_path, top_n=20):
    """Z-score heatmap of top DEGs across samples for one cell type.

    Selects the top *top_n* significant genes by p-value, computes
    mean expression per sample from the raw count matrix, converts to
    z-scores, and renders a clustered heatmap.

    Args:
        adata: AnnData with ``.raw`` and ``obs["sample_id"]`` / ``obs["cell_type"]``.
        deg_df: DataFrame with ``gene``, ``cell_type``, ``comparison``,
            ``pval``, ``significant`` columns.
        cell_type: Cell type to subset.
        comparison: Comparison key (e.g. ``"ref_vs_treat"``).
        output_path: Path to save the PNG figure.
        top_n: Number of top DEGs to display.

    Returns:
        bool: True if the plot was generated, False if too few genes or samples.
    """
    sub = deg_df[(deg_df["cell_type"] == cell_type) &
                 (deg_df["comparison"] == comparison) &
                 (deg_df["significant"])]
    if len(sub) < 3:
        return False
    sub = sub.sort_values("pval").head(top_n)
    genes = sub["gene"].tolist()

    sub_adata = adata[adata.obs["cell_type"] == cell_type]
    raw = sub_adata.raw.to_adata()
    gene_idx = [i for i, g in enumerate(raw.var_names) if g in genes]
    found_genes = [raw.var_names[i] for i in gene_idx]
    if not gene_idx:
        return False

    rows = []
    for sample in sorted(sub_adata.obs["sample_id"].unique()):
        mask = sub_adata.obs["sample_id"] == sample
        sample_cells = raw[mask]
        if sample_cells.X.shape[0] == 0:
            continue
        means = np.array(sample_cells.X[:, gene_idx].mean(axis=0)).flatten()
        rows.append({"sample": sample, **dict(zip(found_genes, means))})
    expr = pd.DataFrame(rows).set_index("sample")
    if expr.shape[0] < 2 or expr.shape[1] < 2:
        return False

    expr_z = (expr - expr.mean()) / (expr.std() + 1e-10)

    fig, ax = plt.subplots(figsize=(10, max(4, len(expr_z) * 0.5)))
    sns.heatmap(expr_z, cmap="RdBu_r", center=0, annot=False,
                linewidths=0.3, ax=ax, cbar_kws={"label": "Z-score"})
    ax.set_title(
        f'{cell_type} top {len(found_genes)} DEGs ({comparison.replace("_vs_", " vs ")})',
        fontsize=11, fontweight="bold",
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    return True


def go_dotplot(go_file, output_path, title, top_n=12):
    """Dotplot of GO enrichment results from clusterProfiler.

    Each dot represents one GO term. Dot size encodes -log10(pvalue)
    and color encodes the pvalue itself. X-axis is gene ratio
    (DE genes in term / total DE genes).

    Args:
        go_file: Path to a clusterProfiler GO CSV with columns
            ``GeneRatio``, ``pvalue``, ``Description``, ``Count``.
        output_path: Path to save the PNG figure.
        title: Plot title.
        top_n: Number of top GO terms (by p-value) to display.

    Returns:
        bool: True if the plot was generated, False if file missing or empty.
    """
    if not os.path.exists(go_file):
        return False
    df = pd.read_csv(go_file)
    if len(df) == 0:
        return False
    df = df.sort_values("pvalue").head(top_n).copy()

    # Parse GeneRatio "k/n" → float
    if "GeneRatio" in df.columns:
        parts = df["GeneRatio"].str.split("/", expand=True)
        df["gene_ratio_num"] = parts[0].astype(float) / parts[1].astype(float)
    else:
        df["gene_ratio_num"] = df["Count"] / 100
    df["nlogp"] = -np.log10(df["pvalue"].clip(1e-300))
    df["Description"] = df["Description"].astype(str).str[:60]

    df["_color"] = df["nlogp"]
    vmin, vmax = df["_color"].min(), df["_color"].max()
    if vmin == vmax:
        vmin, vmax = 0, max(1, vmax)

    fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.4)))
    sc_plot = ax.scatter(
        df["gene_ratio_num"], range(len(df)),
        s=df["nlogp"] * 25, c=df["_color"], cmap="YlOrRd",
        vmin=vmin, vmax=vmax, alpha=0.85,
        edgecolor="black", linewidth=0.5,
    )
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Description"], fontsize=9)
    ax.set_xlabel("Gene Ratio (DE in term / total DE)", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("-log10(pvalue)", fontsize=10)

    # Size legend for dot interpretation
    nlogp_legend_vals = [v for v in [2, 5, 10] if v <= df["nlogp"].max()]
    if nlogp_legend_vals:
        x_legend = df["gene_ratio_num"].max() * 1.15
        for i, nlp in enumerate(nlogp_legend_vals):
            ax.scatter([x_legend], [-(i + 0.5)],
                       s=nlp * 25, c="gray", alpha=0.5,
                       edgecolor="black", linewidth=0.5)
            ax.text(x_legend + df["gene_ratio_num"].max() * 0.05,
                    -(i + 0.5), f"-log10(p)={nlp:.0f}",
                    va="center", fontsize=8)
        ax.text(x_legend - df["gene_ratio_num"].max() * 0.02,
                0.5, "Size:", ha="right", va="center", fontsize=9, fontweight="bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    return True


def deg_overlap_bar(deg_df, sets_to_compare, output_path, title):
    """Horizontal bar plot of DEG set overlaps across comparisons.

    Computes Venn-diagram-style overlaps: genes unique to each comparison,
    shared by each pair, and shared by all (when >= 3 comparisons).

    Args:
        deg_df: DataFrame with ``comparison``, ``significant``, ``gene`` columns.
        sets_to_compare: List of ``(label, comparison_key)`` tuples.
        output_path: Path to save the PNG figure.
        title: Plot title.
    """
    from itertools import combinations

    sig_per_set = {}
    for label, comp in sets_to_compare:
        sub = deg_df[(deg_df["comparison"] == comp) & (deg_df["significant"])]
        sig_per_set[label] = set(sub["gene"])

    if not all(sig_per_set.values()):
        return

    keys = list(sig_per_set.keys())
    n = len(keys)
    overlap_data = {}

    # Genes shared by ALL comparisons (only meaningful when n >= 3)
    if n >= 3:
        common_all = set.intersection(*sig_per_set.values())
        overlap_data["All"] = len(common_all)

    # Pairwise: genes in exactly this pair, NOT in any other set
    for i, j in combinations(range(n), 2):
        label = f"{keys[i]} & {keys[j]}"
        shared = sig_per_set[keys[i]] & sig_per_set[keys[j]]
        for k in range(n):
            if k != i and k != j:
                shared -= sig_per_set[keys[k]]
        overlap_data[label] = len(shared)

    # Unique: genes found ONLY in this comparison, not in any other
    for i in range(n):
        label = f"Only {keys[i]}"
        only = sig_per_set[keys[i]].copy()
        for j in range(n):
            if j != i:
                only -= sig_per_set[keys[j]]
        overlap_data[label] = len(only)

    overlap_data = dict(sorted(overlap_data.items(), key=lambda x: -x[1]))

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = sns.color_palette("tab10", len(overlap_data))
    bars = ax.barh(list(overlap_data.keys()), list(overlap_data.values()),
                   color=colors)
    ax.set_xlabel("Number of significant DEGs", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    for bar, val in zip(bars, overlap_data.values()):
        ax.text(val + 20, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


# ────────────────────────────────────────────────────────────
# 7. Plot orchestration
# ────────────────────────────────────────────────────────────

def plot_all(config, out_dir, deg_dir, go_dir, fig_dir, top_cell_types=None,
             highlight_genes=None):
    """Generate all figures for all tissues defined in the config.

    For each tissue, produces:
      - DEG count heatmap (cell type × comparison)
      - Volcano plots for top cell types per comparison
      - Z-score expression heatmaps for top DEGs
      - GO dotplots (if GO results exist in *go_dir*)
      - DEG overlap bar plot (if >= 2 comparisons)

    Also saves a combined ``all_pseudobulk_DEG.csv``.

    Args:
        config: Tissue config dict (see :func:`load_config`).
        out_dir: Root output directory.
        deg_dir: Directory containing ``deg_{tissue}/`` subdirectories.
        go_dir: Directory containing GO enrichment CSVs.
        fig_dir: Directory for output figures.
        top_cell_types: Optional list of cell types to plot. If None,
            the top 3 by significant DEG count are used.
    """
    import scanpy as sc

    all_deg_frames = []

    for tissue, cfg in config.items():
        cmap = comp_map(cfg["comparisons"])
        prefix = f"{tissue}_"
        tissue_deg_dir = os.path.join(deg_dir, f"deg_{tissue}")

        # Load DEG results for this tissue
        t_deg = load_all_deg(tissue_deg_dir, prefix)
        if len(t_deg) == 0:
            print(f"  [{tissue}] No DEG results found, skipping plots")
            continue
        all_deg_frames.append(t_deg)
        print(f"  [{tissue}] DEG rows: {len(t_deg)}")

        # 1. Count heatmap (cross-comparison → fig_dir)
        os.makedirs(fig_dir, exist_ok=True)
        deg_count_heatmap(
            t_deg, cmap,
            f"{tissue.title()}: Significant DEG count per cell type",
            os.path.join(fig_dir, f"{tissue}_DEG_count_heatmap.png"),
        )

        # 2. Volcano + expression heatmap + GO dotplot (per comparison → comp_dir/ct_dir)
        adata = sc.read_h5ad(cfg["h5ad"])

        for comp in cfg["comparisons"]:
            key = comp["key"]
            label = comp["label"]
            comp_dir = os.path.join(tissue_deg_dir, f"{tissue}_{key}")

            # Pick cell types: user-specified, or top 3 by sig count,
            # or all cell types if no significant genes found
            if top_cell_types:
                cts = top_cell_types
            else:
                sig_count = t_deg[(t_deg["significant"]) &
                                  (t_deg["comparison"] == key)].groupby("cell_type").size()
                cts = sig_count.sort_values(ascending=False).head(3).index.tolist()
                if not cts:
                    cts = sorted(t_deg.loc[t_deg["comparison"] == key, "cell_type"].unique().tolist())

            for ct in cts:
                ct_safe = ct.replace(" ", "_")
                ct_dir = os.path.join(comp_dir, ct_safe)
                os.makedirs(ct_dir, exist_ok=True)
                comp_short = label.replace("_", " ")
                title = f"{ct} ({comp_short})"
                prefix = f"{tissue}_{key}_{ct_safe}"

                ok = volcano_plot(t_deg, ct, key,
                                  os.path.join(ct_dir, f"{prefix}_volcano.png"),
                                  title=title,
                                  highlight_genes=highlight_genes)
                if ok:
                    print(f"    volcano: {tissue}/{ct}/{label}")

                ok = top_deg_heatmap(adata, t_deg, ct, key,
                                     os.path.join(ct_dir, f"{prefix}_heatmap.png"))
                if ok:
                    print(f"    heatmap: {tissue}/{ct}/{label}")

            # GO dotplots → ct_dir
            if os.path.isdir(go_dir):
                for go_file in sorted(os.listdir(go_dir)):
                    if not go_file.endswith(".csv"):
                        continue
                    if "_go_up.csv" not in go_file and "_go_down.csv" not in go_file:
                        continue
                    if "_kegg_" in go_file:
                        continue
                    # Match: {tissue}_{key}_{cell_type}_go_{dir}.csv
                    expected_prefix = f"{tissue}_{key}_"
                    if not go_file.startswith(expected_prefix):
                        continue

                    base = go_file.replace(".csv", "")
                    direction = "UP" if base.endswith("_go_up") else "DOWN"
                    ct_part = base[len(expected_prefix):]
                    ct_name = ct_part.replace("_go_up", "").replace("_go_down", "")
                    ct_go_dir = os.path.join(comp_dir, ct_name)
                    os.makedirs(ct_go_dir, exist_ok=True)
                    ct_clean = ct_name.replace("_", " ")
                    title = f"{tissue.title()} {ct_clean} GO {direction}"
                    out = os.path.join(ct_go_dir, f"go_dotplot_{direction.lower()}.png")
                    ok = go_dotplot(os.path.join(go_dir, go_file), out, title, top_n=12)
                    if ok:
                        print(f"    GO: {base}")

        # 3. Overlap bar plot (cross-comparison → fig_dir)
        if len(cfg["comparisons"]) >= 2:
            sets = [(c["label"], c["key"]) for c in cfg["comparisons"]]
            deg_overlap_bar(
                t_deg, sets,
                os.path.join(fig_dir, f"{tissue}_DEG_overlap.png"),
                f"{tissue.title()}: DEG overlap across {len(sets)} comparisons",
            )

    # Save combined DEG table across all tissues
    if all_deg_frames:
        all_deg = pd.concat(all_deg_frames, ignore_index=True)
        all_deg.to_csv(os.path.join(deg_dir, "all_pseudobulk_DEG.csv"), index=False)
        print(f"\nCombined DEG: {len(all_deg)} rows -> {deg_dir}/all_pseudobulk_DEG.csv")


# ────────────────────────────────────────────────────────────
# 8. Main
# ────────────────────────────────────────────────────────────

def parse_cli_args(argv=None):
    """Parse CLI arguments and build the internal config dict.

    Supports two input modes:
      1. CLI-only: ``--h5ad tissue=path`` + ``--compare tissue:ref:treat[=label]``
      2. YAML file: ``--comp-config comparisons.yaml``

    In CLI mode, the label is optional and defaults to ``ref_vs_treat``.

    Args:
        argv: Argument list (defaults to ``sys.argv[1:]``).

    Returns:
        tuple:
            - **args** (*argparse.Namespace*): Parsed CLI arguments.
            - **config** (*dict*): Tissue config dict compatible with
              :func:`load_config` output format.
    """
    parser = argparse.ArgumentParser(
        description=__doc__ + "\n\n"
        "PARAMETER CONVENTION (READ FIRST):\n"
        "  --compare tissue:ref:treat[=label]   ← ref FIRST, treat SECOND\n"
        "  - File names: ref_vs_treat  (e.g. uterus_luanchao-21310_vs_luanchao-11238)\n"
        "  - Plot titles: 'ref vs treat'  (ref on left = baseline,\n"
        "    treat on right = perturbation)\n"
        "  - CSV 'comparison' column: ref_vs_treat  (e.g. 'luanchao-21310_vs_luanchao-11238')\n"
        "  - log2FC = log2(treat / ref); positive = higher in treat.\n"
        "  Convention follows 'reference (control) first, treatment second',\n"
        "  matching how biological contrasts are typically declared.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Ovaries — youth (ref) vs aged (treat), default label 'ref_vs_treat':
  python pseudobulk.py --out-dir ./results \\
    --h5ad ovaries=/path/to/ov.h5ad \\
    --compare ovaries:luanchao-21310:luanchao-11238

  # Uterus — sham (ref) vs OV-POI (treat), custom labels:
  python pseudobulk.py --out-dir ./results \\
    --h5ad ovaries=/path/to/ov.h5ad \\
    --h5ad uterus=/path/to/ut.h5ad \\
    --compare ovaries:luanchao-21310:luanchao-11238=Aging_vs_Youth \\
    --compare uterus:zigong-21310:ZIGONG-21224=Sham_vs_OV-POI

  # Or use a YAML config:
  python pseudobulk.py --out-dir ./results --comp-config comparisons.yaml
""")
    parser.add_argument("--comp-config",
                        help="YAML config (alternative to --h5ad/--compare)")
    parser.add_argument("--h5ad", action="append", metavar="tissue=path",
                        help="Tissue name = h5ad path (repeatable)")
    parser.add_argument("--compare", action="append", metavar="tissue:ref:treat[=label]",
                        help="Comparison (repeatable, label optional)")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-cells", type=int, default=30,
                        help="Min cells per pseudobulk sample (default: 30)")
    parser.add_argument("--min-counts", type=int, default=1000,
                        help="Min total UMIs per pseudobulk sample (default: 1000)")
    parser.add_argument("--n-pseudo-reps", type=int, default=3,
                        help="Split single-sample conditions into N pseudo-replicates for DESeq2 (default: 3, 0=disabled)")
    parser.add_argument("--top-cell-types", nargs="+", default=None)
    parser.add_argument("--skip-deseq2", action="store_true",
                        help="Skip PyDESeq2 DEG step")
    parser.add_argument("--skip-plot", action="store_true",
                        help="Skip plot generation step")
    parser.add_argument("--skip-prepare", action="store_true",
                        help="Skip pseudobulk preparation step")
    parser.add_argument("--skip-qc", action="store_true",
                        help="Skip QC/exploration plots (PCA, filter, summary)")
    parser.add_argument("--highlight-genes", nargs="+", default=None,
                        help="Genes to highlight on volcano plots (e.g. CDKN2A TP53 CCL2)")
    args = parser.parse_args(argv)

    # Load config from YAML or CLI args
    if args.comp_config:
        config = load_config(args.comp_config)
    elif args.h5ad and args.compare:
        config = {}
        # Parse --h5ad tissue=path
        for item in args.h5ad:
            if "=" not in item:
                parser.error(f"--h5ad format: tissue=path, got: {item}")
            tissue, path = item.split("=", 1)
            config[tissue] = {"h5ad": path, "comparisons": []}
        # Parse --compare tissue:ref:treat[=label]
        for item in args.compare:
            # Split on = for label (only the last =)
            if "=" in item:
                pair, label = item.rsplit("=", 1)
            else:
                pair, label = item, None
            parts = pair.split(":")
            if len(parts) != 3:
                parser.error(f"--compare format: tissue:ref:treat[=label], got: {item}")
            tissue, ref, treat = parts
            if tissue not in config:
                parser.error(f"Unknown tissue '{tissue}', define it with --h5ad first")
            key = f"{ref}_vs_{treat}"
            config[tissue]["comparisons"].append({
                "ref": ref, "treat": treat,
                "label": label or key, "key": key,
            })
        if not any(cfg["comparisons"] for cfg in config.values()):
            parser.error("No comparisons defined. Use --compare tissue:ref:treat")
    else:
        parser.error("Provide either --comp-config or both --h5ad and --compare")

    return args, config


def main():
    """Entry point: prepare counts → PyDESeq2 DEG → generate figures."""
    args, config = parse_cli_args()

    counts_dir = os.path.join(args.out_dir, "pseudobulk_counts")
    deg_dir = os.path.join(args.out_dir, "deg_pseudobulk")
    go_dir = os.path.join(args.out_dir, "go_macaque")
    fig_dir = os.path.join(args.out_dir, "figures_pseudobulk")

    # ── Step 1: Prepare counts ──
    if not args.skip_prepare:
        print("=" * 60)
        print("STEP 1: Prepare pseudo-bulk counts")
        print("=" * 60)
        all_written = []
        tissue_adata = {}
        for tissue, cfg in config.items():
            written, adata = prepare_counts(
                cfg["h5ad"], counts_dir, tissue,
                min_cells=args.min_cells, min_counts=args.min_counts,
                n_pseudo_reps=args.n_pseudo_reps,
            )
            all_written.extend(written)
            tissue_adata[tissue] = adata
        if all_written:
            pd.DataFrame(all_written).to_csv(
                os.path.join(counts_dir, "pseudobulk_summary.csv"), index=False)
            print(f"\nTotal count files: {len(all_written)}")

            # QC plots (PCA, sample filtering, summary)
            if not args.skip_qc:
                print("\n--- QC / Exploration plots ---")
                os.makedirs(fig_dir, exist_ok=True)
                pseudobulk_summary_plot(all_written,
                                        os.path.join(fig_dir, "pseudobulk_summary.png"))
                for tissue, cfg in config.items():
                    if tissue not in tissue_adata:
                        continue
                    adata = tissue_adata[tissue]
                    # Build condition map from comparisons
                    cond_map = {}
                    for comp in cfg["comparisons"]:
                        cond_map[comp["ref"]] = "reference"
                        cond_map[comp["treat"]] = "treatment"
                    pca_exploration_plot(
                        adata,
                        os.path.join(fig_dir, f"{tissue}_pca_exploration.png"),
                        condition_map=cond_map,
                    )
                    sample_filter_plot(
                        adata,
                        os.path.join(fig_dir, f"{tissue}_sample_filter.png"),
                        min_cells=args.min_cells, min_counts=args.min_counts,
                    )
    else:
        tissue_adata = {}

    # ── Step 2: PyDESeq2 DEG ──
    if not args.skip_deseq2:
        print("\n" + "=" * 60)
        print("STEP 2: PyDESeq2 DEG analysis")
        print("=" * 60)
        for tissue, cfg in config.items():
            t_deg_dir = os.path.join(deg_dir, f"deg_{tissue}")
            counts_dir_tissue = counts_dir

            # Load all counts CSVs for this tissue
            import glob
            pattern = os.path.join(counts_dir_tissue, f"{tissue}_*_counts.csv")
            counts_files = sorted(glob.glob(pattern))

            if not counts_files:
                print(f"  [{tissue}] No counts files found, skipping")
                continue

            print(f"\n[{tissue}] Found {len(counts_files)} counts files")

            for counts_file in counts_files:
                basename = os.path.basename(counts_file)
                ct_safe = basename.replace(f"{tissue}_", "").replace("_counts.csv", "")

                counts_df = pd.read_csv(counts_file, index_col=0)
                cell_type = counts_df["cell_type"].iloc[0] if "cell_type" in counts_df.columns else ct_safe

                for comp in cfg["comparisons"]:
                    ref_id, treat_id = comp["ref"], comp["treat"]
                    key = comp["key"]
                    label = comp["label"]

                    # Output directory: deg_{tissue}/{tissue}_{ref}_vs_{treat}/{cell_type}/
                    comp_dir = os.path.join(t_deg_dir, f"{tissue}_{key}")
                    ct_dir = os.path.join(comp_dir, ct_safe)
                    os.makedirs(ct_dir, exist_ok=True)

                    output_prefix = os.path.join(
                        ct_dir, f"{tissue}_{key}_{ct_safe}")

                    # Run PyDESeq2
                    try:
                        deg_results = run_pydeseq2(counts_df, ref_id, treat_id)
                    except Exception as e:
                        print(f"  ERROR: PyDESeq2 failed for {ct_safe} "
                              f"{ref_id} vs {treat_id}: {e}", file=sys.stderr)
                        continue
                    if len(deg_results) == 0:
                        continue

                    # Save DEG CSV
                    deg_results.to_csv(f"{output_prefix}_DEG.csv", index=False)

                    # Diagnostic plots
                    pydeseq2_diagnostic_plots(
                        deg_results, f"{output_prefix}_diag",
                        treat_id, ref_id,
                    )

    # ── Step 3: Plot ──
    if not args.skip_plot:
        print("\n" + "=" * 60)
        print("STEP 3: Generate figures")
        print("=" * 60)
        plot_all(config, args.out_dir, deg_dir, go_dir, fig_dir,
                 args.top_cell_types, args.highlight_genes)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
    if os.path.isdir(fig_dir):
        for root, dirs, files in os.walk(fig_dir):
            for f in sorted(files):
                if f.endswith(".png"):
                    fpath = os.path.join(root, f)
                    size = os.path.getsize(fpath)
                    rel = os.path.relpath(fpath, fig_dir)
                    print(f"  {rel}: {size/1024:.0f}KB")


if __name__ == "__main__":
    main()
