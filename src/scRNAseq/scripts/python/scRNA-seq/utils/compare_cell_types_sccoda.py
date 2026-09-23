#!/usr/bin/env python3
"""Differential abundance analysis using scCODA (sc-best-practices recommended).

Replaces compare_cell_types.py (Fisher's exact + chi2) with the Bayesian
Dirichlet-Multinomial model from Büttner et al. 2021 (Nat Comm).

Why scCODA over Fisher's exact:
  - Jointly models all cell-type proportions (compositional constraint).
  - Reference cell type auto-selected (least dispersion across samples).
  - Bayesian credible-effect output, designed for low-replicate data (n=2-5).
  - Supports covariate formula (batch correction, multi-factor designs).

Outputs per comparison:
  - {prefix}_{id1}_vs_{id2}_da.csv        credible effects + final_prob
  - {prefix}_effects_barplot.png         per-cell-type logFC bar
  - {prefix}_effects_umap.png            effects on UMAP
  - {prefix}_rel_abundance_disp.png      dispersion (for reference diagnostic)

Reference: https://www.sc-best-practices.org/conditions/compositional/
"""
import argparse
import logging
import warnings
from pathlib import Path
from typing import List, Tuple, Optional

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import pertpy as pt
import pertpy.tools  # noqa: F401  -- ensure jax deps register

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger(__name__)

def plot_effects_from_csv(csv_path: Path, out_path: Path, title: str) -> bool:
    """Recreate pertpy-style effects barplot from a saved DA CSV.

    Each cell type gets its own unique color (N colors for N cell types).
    Credible effects have black edges. Labels are horizontal for readability.

    Args:
        csv_path: Path to ``*_da.csv`` with columns
            ``cell_type, log2FC, is_credible, final_prob``.
        out_path: Path to save the PNG.
        title: Plot title.

    Returns:
        True if plot was generated, False if CSV is empty.
    """
    df = pd.read_csv(csv_path)
    if len(df) == 0:
        return False
    df = df.sort_values("log2FC").reset_index(drop=True)

    n = len(df)
    # Color only credible effects; non-credible = gray
    n_credible = int(df["is_credible"].sum())
    cmap = plt.cm.get_cmap("tab20" if n_credible <= 20 else "tab20b")
    cred_colors = [cmap(i / max(n_credible - 1, 1)) for i in range(n_credible)]
    ci = 0
    colors = []
    for _, row in df.iterrows():
        if row["is_credible"]:
            colors.append(cred_colors[ci])
            ci += 1
        else:
            colors.append("#cccccc")

    fig, ax = plt.subplots(figsize=(10, max(4, n * 0.45 + 1.5)))

    for i, (_, row) in enumerate(df.iterrows()):
        edgecolor = "black" if row["is_credible"] else "none"
        lw = 0.8 if row["is_credible"] else 0
        ax.barh(i, row["log2FC"], color=colors[i], edgecolor=edgecolor,
                linewidth=lw, height=0.7, zorder=2)

    ax.set_yticks(range(n))
    ax.set_yticklabels(df["cell_type"], fontsize=10)
    ax.set_xlabel("log2-fold change", fontsize=12)
    ax.set_title(title, fontsize=14, y=0.97)
    ax.axvline(0, color="black", linewidth=0.5, alpha=0.5, zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=10)

    fig.subplots_adjust(left=0.25, right=0.98, top=0.88, bottom=0.10)
    fig.savefig(out_path, dpi=300, facecolor="white")
    plt.close(fig)
    log.info("Saved: %s", out_path)
    return True


def run_sccoda(adata: sc.AnnData, ref: str, treat: str, tissue: str,
               formula: str = "condition",
               ref_celltype: str = "automatic",
               out_dir: Optional[Path] = None) -> pd.DataFrame:
    """Run scCODA for one ref vs treat comparison.

    Args:
        adata: AnnData with obs['sample_id'], obs['condition'], obs['cell_type'].
        ref: reference sample id.
        treat: treatment sample id.
        formula: design formula (default: 'condition').
        ref_celltype: 'automatic' or explicit cell type name.
        out_dir: directory to save plots and CSV.

    Returns:
        DataFrame with per-cell-type credible effect results.
    """
    # Subset to two conditions and rebuild count matrix.
    # Order matters for scCODA: the FIRST level of `condition` becomes the
    # base category (no [T.] marker in covariate names); the SECOND level
    # becomes the treatment ([T.<level>] in covariate names).
    # We force `categories=[ref, treat]` so the treatment is `treat` and
    # effects are reported as "ref vs treat" (matches our CLI convention).
    sub = adata[adata.obs["sample_id"].isin([ref, treat])].copy()
    sub.obs["condition"] = pd.Categorical(
        sub.obs["sample_id"].values,
        categories=[ref, treat],
    )
    log.info("Subsetting %d cells (ref=%s, treat=%s)",
             len(sub), ref, treat)

    sccoda = pt.tl.Sccoda()

    # Step 1: load — aggregate cell-level -> sample-level MuData
    log.info("Aggregating cell-level to sample-level (cell_type × sample)...")
    mdata = sccoda.load(
        sub,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier="cell_type",
        sample_identifier="sample_id",
        covariate_obs=["condition"],
    )
    log.info("  sample-level shape: %s × %s",
             mdata["coda"].shape[0], mdata["coda"].shape[1])

    # Step 2: prepare — design matrix + automatic reference selection
    sccoda.prepare(
        mdata,
        formula="condition",
        reference_cell_type=ref_celltype,
    )

    # Step 3: MCMC sampling
    log.info("Running scCODA NUTS sampler (num_samples=2000, num_warmup=500, ~20-40s)...")
    sccoda.run_nuts(
        mdata,
        num_samples=2000,
        num_warmup=500,
        rng_key=0,
    )

    # Step 4: extract credible effects
    summary = sccoda.get_effect_df(mdata)
    summary = summary.reset_index().rename(columns={
        "Cell Type": "cell_type",
        "log2-fold change": "log2FC",
        "Inclusion probability": "final_prob",
    })
    summary["final_prob"] = summary["final_prob"].astype(float)
    summary["log2FC"] = summary["log2FC"].astype(float)
    summary["is_credible"] = summary["final_prob"] >= 0.90

    # Log effect direction
    n_up = int(((summary["log2FC"] > 0) & summary["is_credible"]).sum())
    n_dn = int(((summary["log2FC"] < 0) & summary["is_credible"]).sum())
    n_total = int(summary["is_credible"].sum())
    ref_used = mdata["coda"].uns["scCODA_params"]["reference_cell_type"]
    log.info("  credible: %d (%d up / %d dn), reference=%s",
             n_total, n_up, n_dn, ref_used)
    log.info("  credible cell types:\n%s",
             summary[summary["is_credible"]][["cell_type", "log2FC", "final_prob"]]
             .to_string(index=False))

    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = f"{out_dir.parent.name}_{ref}_vs_{treat}"
        # CSV
        csv_path = out_dir / f"{prefix}_da.csv"
        summary.to_csv(csv_path, index=False)
        log.info("Saved: %s", csv_path)

        # Build human-readable title (full sample IDs only).
        # Convention: "A vs B" = A is reference (left), B is treatment (right).
        full_title = f"{tissue.upper()} DA: {ref} vs {treat}"

        # Effects barplot — use return_fig=True so we can set our own suptitle
        if n_total > 0:
            try:
                fg = sccoda.plot_effects_barplot(
                    mdata,
                    covariates=formula,
                    parameter="log2-fold change",
                    plot_facets=True,
                    plot_zero_covariate=False,
                    plot_zero_cell_type=False,
                    figsize=(10, 4),
                    return_fig=True,
                )
                if hasattr(fg, "fig"):
                    # Force the actual canvas size — pertpy leaves a tall FacetGrid
                    # even when there's only one facet, so resize manually.
                    fg.fig.set_size_inches(10, 5.5)
                    fg.fig.suptitle(full_title, fontsize=14, y=0.97)
                    for ax in fg.axes.flat:
                        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
                        ax.set_title("")
                        ax.set_ylabel(ax.get_ylabel(), fontsize=12)
                    fg.fig.subplots_adjust(
                        left=0.08, right=0.98, top=0.85, bottom=0.38,
                    )
                    fg.fig.savefig(
                        out_dir / f"{prefix}_effects_barplot.png",
                        dpi=300, facecolor="white",
                    )
                    plt.close(fg.fig)
                    log.info("Saved: %s", out_dir / f"{prefix}_effects_barplot.png")
            except Exception as e:
                log.warning("plot_effects_barplot failed: %s", e)
                plt.close("all")
        else:
            log.info("  skip effects_barplot: no credible effects to plot")

        # Stacked barplot of cell-type proportions per sample
        try:
            fig = sccoda.plot_stacked_barplot(
                mdata, feature_name="samples", figsize=(10, 4.5), return_fig=True,
            )
            if fig is None:
                # pertpy sometimes still leaves the figure open — grab current
                fig = plt.gcf()
            fig.set_size_inches(10, 5.0)
            fig.suptitle(full_title, fontsize=14, y=0.97)
            fig.subplots_adjust(left=0.08, right=0.78, top=0.86, bottom=0.22)
            fig.savefig(
                out_dir / f"{prefix}_stacked_barplot.png",
                dpi=300, facecolor="white",
            )
            plt.close(fig)
            log.info("Saved: %s", out_dir / f"{prefix}_stacked_barplot.png")
        except Exception as e:
            log.warning("plot_stacked_barplot failed: %s", e)
            plt.close("all")

    return summary


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "scCODA-based differential abundance analysis.\n"
            "\n"
            "PARAMETER CONVENTION (READ FIRST):\n"
            "  --compare tissue:ref:treat   ← ref FIRST, treat SECOND\n"
            "  In every output:\n"
            "    - File names: ref_vs_treat  (e.g. luanchao-21310_vs_luanchao-11238)\n"
            "    - Plot titles: 'ref vs treat'  (ref on the LEFT, treat on the RIGHT)\n"
            "  In the CSV covariate column, treat is tagged with [T.]: [T.treat]\n"
            "  so log2FC = log2(treat / ref).  ref = baseline / control / unperturbed;\n"
            "  treat = experimental / treated / perturbed.  Convention follows\n"
            "  'reference appears first, treatment second'."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Ovaries — youth (ref) vs aged (treat):
  python compare_cell_types_sccoda.py \\
    --h5ad ovaries=/path/ovaries.h5ad \\
    --compare ovaries:luanchao-21310:luanchao-11238 \\
    --out-dir ./da_sccoda

  # Uterus — sham (ref) vs POI (treat):
  python compare_cell_types_sccoda.py \\
    --h5ad uterus=/path/uterus.h5ad \\
    --compare uterus:zigong-21310:ZIGONG-21224 \\
    --out-dir ./da_sccoda
""")
    p.add_argument("--h5ad", action="append", required=True,
                   metavar="tissue=path",
                   help="Tissue name = h5ad path (repeatable)")
    p.add_argument("--compare", action="append", required=True,
                   metavar="tissue:ref:treat",
                   help="Comparison ref vs treat (repeatable)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--ref-celltype", default="automatic",
                   help="Reference cell type ('automatic' or explicit name)")
    p.add_argument("--formula", default="condition",
                   help="Design formula (default: condition)")
    p.add_argument("--credible-threshold", type=float, default=0.90,
                   help="Posterior probability threshold for credible effect")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Parse --h5ad
    tissues = {}
    for item in args.h5ad:
        if "=" not in item:
            raise SystemExit(f"--h5ad format: tissue=path, got: {item}")
        tissue, path = item.split("=", 1)
        tissues[tissue] = Path(path)

    # Parse --compare
    comparisons = []
    for item in args.compare:
        if ":" not in item:
            raise SystemExit(f"--compare format: tissue:ref:treat, got: {item}")
        tissue, ref, treat = item.split(":")
        if tissue not in tissues:
            raise SystemExit(f"Unknown tissue '{tissue}', define --h5ad first")
        comparisons.append((tissue, ref, treat))

    # Run
    for tissue, ref, treat in comparisons:
        log.info("=" * 60)
        log.info("[%s] %s vs %s", tissue, ref, treat)
        log.info("=" * 60)
        adata = sc.read_h5ad(tissues[tissue])
        log.info("Loaded %s: shape=%s", tissue, adata.shape)

        tissue_out = out_dir / tissue
        run_sccoda(adata, ref, treat, tissue=tissue,
                   formula=args.formula,
                   ref_celltype=args.ref_celltype,
                   out_dir=tissue_out)

    log.info("DONE")


if __name__ == "__main__":
    main()
