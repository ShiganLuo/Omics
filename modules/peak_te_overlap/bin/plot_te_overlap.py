#!/usr/bin/env python3
"""Generate a grouped bar chart of TE family overlap ratios.

Reads per-sample TE overlap count TSVs and produces:
  1. A grouped bar chart (PNG) comparing IP vs Input samples per TE class
  2. A combined TSV with all samples' overlap data
"""
import argparse
import logging
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("plot_te_overlap")

DPI = 300


def load_tsvs(tsv_paths):
    """Load and concatenate per-sample overlap TSVs."""
    dfs = []
    for path in tsv_paths:
        df = pd.read_csv(path, sep="\t")
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def make_pairs(sample_ip_input_map):
    """Parse IP:Input pairs from semicolon-separated string."""
    pairs = {}
    for pair in sample_ip_input_map.split(";"):
        if ":" in pair:
            ip, inp = pair.split(":", 1)
            pairs[ip.strip()] = inp.strip()
    return pairs


def plot_grouped_bar(df, pairs, output_path):
    """
    Create a grouped bar chart:
    - X-axis: TE families (sorted by overlap ratio)
    - Y-axis: overlap ratio (peaks overlapping TE / total peaks)
    - Bars: grouped by IP (solid) vs Input (hatched) for each protein
    """
    # Get all TE classes sorted by mean overlap ratio
    te_classes = (
        df.groupby("te_class")["overlap_ratio"]
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    # Limit to top 15 TE classes for readability
    te_classes = te_classes[:15]

    # Filter data
    df_plot = df[df["te_class"].isin(te_classes)].copy()

    # Build protein groups: each IP + its matched Input
    protein_groups = []
    for ip, inp in pairs.items():
        protein_name = ip.replace("IP", "").replace("Input", "")
        protein_groups.append({
            "name": protein_name,
            "ip": ip,
            "input": inp,
        })

    if not protein_groups:
        # No pairs, just plot all samples
        samples = df_plot["sample_id"].unique().tolist()
        protein_groups = [{"name": s, "ip": s, "input": None} for s in samples]

    n_classes = len(te_classes)
    n_groups = len(protein_groups)
    bar_width = 0.8 / (n_groups * 2)

    fig, ax = plt.subplots(figsize=(max(10, n_classes * 0.8), 6))

    cmap = matplotlib.colormaps["Set2"]
    colors = cmap(np.linspace(0, 1, n_groups))

    for i, group in enumerate(protein_groups):
        ip_data = df_plot[df_plot["sample_id"] == group["ip"]]
        ip_ratios = []
        for tc in te_classes:
            row = ip_data[ip_data["te_class"] == tc]
            ip_ratios.append(row["overlap_ratio"].values[0] if len(row) > 0 else 0.0)

        x_base = np.arange(n_classes)
        offset_ip = (i * 2) * bar_width - (n_groups * bar_width) + bar_width / 2

        bars_ip = ax.bar(
            x_base + offset_ip,
            ip_ratios,
            bar_width,
            label=f"{group['name']} IP",
            color=colors[i],
            edgecolor="black",
            linewidth=0.5,
        )

        if group["input"]:
            inp_data = df_plot[df_plot["sample_id"] == group["input"]]
            inp_ratios = []
            for tc in te_classes:
                row = inp_data[inp_data["te_class"] == tc]
                inp_ratios.append(row["overlap_ratio"].values[0] if len(row) > 0 else 0.0)

            offset_inp = (i * 2 + 1) * bar_width - (n_groups * bar_width) + bar_width / 2
            ax.bar(
                x_base + offset_inp,
                inp_ratios,
                bar_width,
                label=f"{group['name']} Input",
                color=colors[i],
                edgecolor="black",
                linewidth=0.5,
                alpha=0.5,
                hatch="//",
            )

    ax.set_xlabel("TE Family", fontsize=12)
    ax.set_ylabel("Peak Overlap Ratio", fontsize=12)
    ax.set_title("Peak-TE Family Overlap: IP vs Input", fontsize=14)
    ax.set_xticks(np.arange(n_classes))
    ax.set_xticklabels(te_classes, rotation=45, ha="right", fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.set_ylim(0, 1.0)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Chart saved to {output_path}")


def main():
    p = argparse.ArgumentParser("Plot TE family overlap chart")
    p.add_argument("--tsvs", required=True, help="Comma-separated per-sample TSV paths")
    p.add_argument("--output", required=True, help="Output chart PNG path")
    p.add_argument("--combined-tsv", required=True, help="Output combined TSV path")
    p.add_argument("--ip-input-pairs", default="", help="IP:Input pairs, semicolon-separated")
    args = p.parse_args()

    tsv_paths = [p.strip() for p in args.tsvs.split(",") if p.strip()]
    logger.info(f"Loading {len(tsv_paths)} TSV files")
    df = load_tsvs(tsv_paths)

    # Save combined TSV
    df.to_csv(args.combined_tsv, sep="\t", index=False)
    logger.info(f"Combined TSV saved to {args.combined_tsv}")

    # Parse IP:Input pairs
    pairs = {}
    if args.ip_input_pairs:
        pairs = make_pairs(args.ip_input_pairs)
    logger.info(f"IP:Input pairs: {pairs}")

    # Plot
    plot_grouped_bar(df, pairs, args.output)


if __name__ == "__main__":
    main()
