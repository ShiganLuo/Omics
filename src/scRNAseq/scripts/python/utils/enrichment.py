#!/usr/bin/env python3
"""Post-DEG enrichment analysis: GO/KEGG for gene DEGs + TE subfamily for TE DEGs.

Iterates over DEG CSV files from pseudobulk.py output:
  - Gene mode: GO/KEGG enrichment via clusterProfiler (Apptainer)
  - TE mode: TE subfamily-level aggregation + plots

Usage:
    python enrichment.py \\
        --deg-dir results/DEG/deg_pseudobulk_gene \\
        --te-dir results/DEG/deg_pseudobulk_te \\
        --te-bed /path/to/rheMac10_rmsk_TE.bed \\
        --go-kegg-script workflow/Omics/modules/function/bin/go-kegg.r \\
        --container /path/to/function.sif \\
        [--lfc-cut 1] [--p-cut 0.05] [--skip-existing]
"""
import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


# ── TE BED parser ─────────────────────────────────────────────────────

def parse_te_bed(bed_path: str) -> dict[str, tuple[str, str]]:
    """Parse TE BED: gene_id -> (family_id, class_id)."""
    result: dict[str, tuple[str, str]] = {}
    with open(bed_path, encoding="utf-8") as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            gene_id = parts[4]
            if gene_id not in result:
                result[gene_id] = (parts[5], parts[6])
    return result


# ── GO/KEGG enrichment ───────────────────────────────────────────────

def run_go_kegg_for_deg(deg_csv: Path, go_kegg_script: str,
                        container: str, lfc_cut: float, p_cut: float,
                        skip_existing: bool) -> bool:
    """Run GO/KEGG enrichment for one gene DEG CSV. Returns True if ran."""
    ct_dir = deg_csv.parent
    go_out = ct_dir / "go_kegg"

    if skip_existing and (go_out / "go_up.csv").exists() and \
       (go_out / "go_down.csv").exists() and \
       (go_out / "kegg_up.csv").exists() and (go_out / "kegg_down.csv").exists():
        logger.info("  SKIP (exists): %s", go_out)
        return False

    go_out.mkdir(parents=True, exist_ok=True)

    cmd = [
        "apptainer", "exec", "--cleanenv", "--no-home",
        "-B", f"{os.getcwd()}:{os.getcwd()}",
        container, "Rscript", go_kegg_script,
        "-i", str(deg_csv),
        "-o", str(go_out),
        "-s", "macaque",
        "--gene-col", "gene",
        "--value-col", "log2FC",
        "--p-col", "pval_adj",
        "--lfc-cut", str(lfc_cut),
        "--p-cut", str(p_cut),
    ]

    logger.info("  GO/KEGG: %s", go_out)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("  GO/KEGG FAILED: %s\n%s", deg_csv, result.stderr[-500:])
        return False
    return True


def run_go_kegg_all(deg_dir: Path, go_kegg_script: str,
                    container: str, lfc_cut: float, p_cut: float,
                    skip_existing: bool) -> int:
    """Run GO/KEGG for all gene DEG CSVs. Returns count of processed files."""
    deg_files = sorted(deg_dir.rglob("*_DEG.csv"))
    logger.info("Found %d gene DEG files in %s", len(deg_files), deg_dir)

    count = 0
    for deg_csv in deg_files:
        if run_go_kegg_for_deg(deg_csv, go_kegg_script, container,
                               lfc_cut, p_cut, skip_existing):
            count += 1
    logger.info("GO/KEGG processed: %d / %d", count, len(deg_files))
    return count


# ── TE subfamily analysis ─────────────────────────────────────────────

def te_subfamily_analysis(deg_csv: Path, te_map: dict[str, tuple[str, str]],
                          lfc_cut: float, p_cut: float,
                          skip_existing: bool) -> bool:
    """Aggregate TE DEGs by subfamily, generate plots. Returns True if ran."""
    ct_dir = deg_csv.parent
    out_dir = ct_dir / "te_subfamily"
    summary_csv = out_dir / "te_subfamily_summary.csv"

    if skip_existing and summary_csv.exists():
        logger.info("  SKIP (exists): %s", out_dir)
        return False

    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(deg_csv)
    if "gene" not in df.columns or "log2FC" not in df.columns:
        logger.warning("  Invalid DEG CSV: %s", deg_csv)
        return False

    # Map TE genes to family/class
    df["family_id"] = df["gene"].map(lambda g: te_map.get(g, ("", ""))[0])
    df["class_id"] = df["gene"].map(lambda g: te_map.get(g, ("", ""))[1])

    # Filter out unmapped genes
    mapped = df[df["family_id"] != ""].copy()
    if len(mapped) == 0:
        logger.info("  No TE genes mapped: %s", deg_csv)
        return False

    # Aggregate by family
    family_stats = []
    for fam, grp in mapped.groupby("family_id"):
        n_total = len(grp)
        n_sig = int(grp["significant"].sum()) if "significant" in grp.columns else 0
        mean_lfc = grp["log2FC"].mean()
        max_lfc = grp["log2FC"].max()
        min_lfc = grp["log2FC"].min()
        mean_base = grp["baseMean"].mean() if "baseMean" in grp.columns else 0
        cls = grp["class_id"].iloc[0]
        family_stats.append({
            "family": fam, "class_id": cls,
            "n_genes": n_total, "n_sig": n_sig,
            "mean_log2FC": mean_lfc, "max_log2FC": max_lfc, "min_log2FC": min_lfc,
            "mean_baseMean": mean_base,
        })

    fam_df = pd.DataFrame(family_stats).sort_values("mean_log2FC", ascending=False)
    fam_df.to_csv(summary_csv, index=False)
    logger.info("  TE subfamily summary: %d families, %d with sig genes",
                len(fam_df), (fam_df.n_sig > 0).sum())

    # Filter for plotting: families with sig genes or high mean_baseMean
    plot_df = fam_df[(fam_df.n_sig > 0) | (fam_df.mean_baseMean > 0.01)].copy()
    if len(plot_df) < 3:
        plot_df = fam_df.nlargest(16, "mean_baseMean").copy()

    plot_df = plot_df.sort_values("mean_log2FC")

    # ── Plot 1: family-level log2FC bar chart ──
    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(plot_df))))
    colors = ["#d62728" if v > 0 else "#1f77b4" for v in plot_df["mean_log2FC"]]
    ax.barh(plot_df["family"], plot_df["mean_log2FC"], color=colors)
    ax.axvline(0, color="black", lw=0.5)

    # Mark families with sig genes
    for _, row in plot_df[plot_df.n_sig > 0].iterrows():
        idx = list(plot_df.family).index(row.family)
        ax.text(row.mean_log2FC, idx, f" {row.n_sig}sig",
                va="center", fontsize=7, color="gray")

    ax.set_xlabel("mean log2FC")
    ax.set_title(f"TE Subfamily Changes ({ct_dir.parent.parent.name})")
    fig.tight_layout()
    fig.savefig(out_dir / "te_subfamily_log2FC.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ── Plot 2: class-level summary ──
    class_stats = mapped.groupby("class_id").agg(
        n_genes=("gene", "count"),
        mean_lfc=("log2FC", "mean"),
    ).reset_index().sort_values("mean_lfc", ascending=True)

    if len(class_stats) > 1:
        fig, ax = plt.subplots(figsize=(8, max(3, 0.5 * len(class_stats))))
        colors = ["#d62728" if v > 0 else "#1f77b4" for v in class_stats["mean_lfc"]]
        ax.barh(class_stats["class_id"], class_stats["mean_lfc"], color=colors)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlabel("mean log2FC")
        ax.set_title(f"TE Class Changes ({ct_dir.parent.parent.name})")
        fig.tight_layout()
        fig.savefig(out_dir / "te_class_log2FC.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Save full DEG with family/class annotation
    mapped.to_csv(out_dir / "te_deg_annotated.csv", index=False)

    return True


def run_te_subfamily_all(te_dir: Path, te_bed: str,
                         lfc_cut: float, p_cut: float,
                         skip_existing: bool) -> int:
    """Run TE subfamily analysis for all TE DEG CSVs."""
    te_map = parse_te_bed(te_bed)
    logger.info("TE BED: %d unique TE subfamilies", len(te_map))

    deg_files = sorted(te_dir.rglob("*_DEG.csv"))
    logger.info("Found %d TE DEG files in %s", len(deg_files), te_dir)

    count = 0
    for deg_csv in deg_files:
        if te_subfamily_analysis(deg_csv, te_map, lfc_cut, p_cut, skip_existing):
            count += 1
    logger.info("TE subfamily processed: %d / %d", count, len(deg_files))
    return count


# ── CLI ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Post-DEG enrichment analysis")
    sub = parser.add_subparsers(dest="command")

    # go-kegg subcommand
    p_go = sub.add_parser("go-kegg", help="GO/KEGG enrichment for gene DEGs")
    p_go.add_argument("--deg-dir", type=Path, required=True,
                       help="Gene DEG dir (deg_pseudobulk_gene)")
    p_go.add_argument("--go-kegg-script",
                       default="workflow/Omics/modules/function/bin/go-kegg.r")
    p_go.add_argument("--container",
                       default="/home/luosg/Database/env/function/function.sif")
    p_go.add_argument("--lfc-cut", type=float, default=1.0)
    p_go.add_argument("--p-cut", type=float, default=0.05)
    p_go.add_argument("--skip-existing", action="store_true")

    # te-subfamily subcommand
    p_te = sub.add_parser("te-subfamily", help="TE subfamily aggregation for TE DEGs")
    p_te.add_argument("--te-dir", type=Path, required=True,
                       help="TE DEG dir (deg_pseudobulk_te)")
    p_te.add_argument("--te-bed", type=Path, required=True,
                       help="TE BED file with family_id/class_id")
    p_te.add_argument("--lfc-cut", type=float, default=1.0)
    p_te.add_argument("--p-cut", type=float, default=0.05)
    p_te.add_argument("--skip-existing", action="store_true")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == "go-kegg":
        if not args.deg_dir.exists():
            logger.error("Gene DEG dir not found: %s", args.deg_dir)
            sys.exit(1)
        run_go_kegg_all(args.deg_dir, args.go_kegg_script, args.container,
                        args.lfc_cut, args.p_cut, args.skip_existing)

    elif args.command == "te-subfamily":
        if not args.te_dir.exists():
            logger.error("TE DEG dir not found: %s", args.te_dir)
            sys.exit(1)
        run_te_subfamily_all(args.te_dir, str(args.te_bed),
                             args.lfc_cut, args.p_cut, args.skip_existing)


if __name__ == "__main__":
    main()
