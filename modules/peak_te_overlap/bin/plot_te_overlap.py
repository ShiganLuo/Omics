#!/usr/bin/env python3
"""Generate TE subfamily enrichment plot.

Reads per-sample TE subfamily locus-level TSVs and produces:
  - Per IP:Input pair: separate enrichment bar charts (default)
  - Or combined: single figure with all pairs (with --combine flag)
  - A combined TSV with all samples' subfamily overlap data

Each figure has:
  - Top: line plot of mean TE length per subfamily
  - Bottom: bar chart of log2(IP / Input) enrichment per subfamily
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
    """Load and concatenate per-sample subfamily overlap TSVs."""
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


def assign_condition(sample_id, pairs):
    """Assign 'IP' or 'Input' condition based on sample_ip_input_map."""
    if sample_id in pairs:
        return "IP"
    for ip, inp in pairs.items():
        if sample_id == inp:
            return "Input"
    return "Unknown"


def get_boxplot_column(method):
    """Return the column name to use for enrichment based on method."""
    if method == "interval":
        return "interval_overlap_frac"
    elif method == "count":
        return "overlap_peak_count"
    elif method == "reads":
        return "ip_reads"
    else:
        raise ValueError(f"Unknown method: {method}")


def compute_enrichment_single_pair(df_pair, method):
    """Compute enrichment (IP vs Input) for a single pair's data.

    method='reads': sum ip_reads / sum input_reads per subfamily, then log2
    method='count': normalized peak count ratio
    method='interval': mean interval overlap fraction ratio

    Returns a Series of log2(IP/Input) per subfamily.
    """
    has_reads = "ip_reads" in df_pair.columns and "input_reads" in df_pair.columns

    if method == "reads":
        if not has_reads:
            raise ValueError("method='reads' requires ip_reads/input_reads columns in data. "
                             "Ensure --ip-bam and --input-bam are provided to intersect_te.py")
        # Sum reads per subfamily directly (ip_reads and input_reads are in the same row)
        agg = df_pair.groupby("te_subfamily").agg(
            ip_total=("ip_reads", "sum"),
            input_total=("input_reads", "sum"),
        )
        pseudocount = max(1, min(agg["ip_total"].min(), agg["input_total"].min()) * 0.01) if len(agg) > 0 else 1
        log2fc = np.log2((agg["ip_total"] + pseudocount) / (agg["input_total"] + pseudocount))
        return log2fc

    if method == "count":
        sample_sf = df_pair.groupby(["sample_id", "te_subfamily"]).agg(
            peak_count=("overlap_peak_count", "max"),
            condition=("condition", "first"),
        ).reset_index()
        total_per_sample = sample_sf.groupby("sample_id")["peak_count"].sum()
        sample_sf["norm_count"] = sample_sf.apply(
            lambda r: r["peak_count"] / total_per_sample[r["sample_id"]] if total_per_sample[r["sample_id"]] > 0 else 0,
            axis=1,
        )
        cond_mean = sample_sf.groupby(["condition", "te_subfamily"])["norm_count"].mean().unstack(level=0)

    elif method == "interval":
        sample_sf = df_pair.groupby(["sample_id", "te_subfamily"]).agg(
            mean_frac=("interval_overlap_frac", "mean"),
            condition=("condition", "first"),
        ).reset_index()
        cond_mean = sample_sf.groupby(["condition", "te_subfamily"])["mean_frac"].mean().unstack(level=0)

    else:
        raise ValueError(f"Unknown method: {method}")

    if "IP" not in cond_mean.columns or "Input" not in cond_mean.columns:
        logger.warning("Missing IP or Input conditions, cannot compute enrichment")
        return pd.Series(dtype=float)

    pseudocount = cond_mean[["IP", "Input"]].min().min()
    if pseudocount <= 0:
        pseudocount = 1e-6

    ip_vals = cond_mean["IP"].fillna(0) + pseudocount
    input_vals = cond_mean["Input"].fillna(0) + pseudocount

    log2fc = np.log2(ip_vals / input_vals)
    return log2fc.sort_values(ascending=False)


def select_top_subfamilies(log2fc, df_pair, top_n):
    """Select top N subfamilies: filter by locus count >= 3, then top by enrichment."""
    locus_counts = df_pair.groupby("te_subfamily").size()
    valid = locus_counts[locus_counts >= 3].index
    log2fc = log2fc[log2fc.index.isin(valid)]
    return log2fc.sort_values(ascending=False).index.tolist()[:top_n]


def plot_single_pair(df_pair, pairs, ip_name, input_name, output_path,
                     method="count", sort_by="te_length", top_n=30):
    """Plot enrichment for a single IP:Input pair."""
    # Compute enrichment
    log2fc = compute_enrichment_single_pair(df_pair, method)
    if log2fc.empty:
        logger.warning(f"No enrichment data for {ip_name} vs {input_name}")
        return

    # Select top N
    sorted_subfamilies = select_top_subfamilies(log2fc, df_pair, top_n)
    log2fc = log2fc.loc[sorted_subfamilies]

    # Sort display
    mean_length = df_pair.groupby("te_subfamily")["te_length"].mean()
    if sort_by == "te_length":
        common = log2fc.index.intersection(mean_length.index)
        sorted_subfamilies = mean_length.loc[common].sort_values().index.tolist()
        log2fc = log2fc.reindex(sorted_subfamilies)
    # else: already sorted by enrichment descending

    mean_length = mean_length.reindex(sorted_subfamilies)

    # Build figure
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1,
        figsize=(max(12, len(sorted_subfamilies) * 0.5), 8),
        gridspec_kw={"height_ratios": [1, 2]},
        sharex=True,
    )

    x_pos = np.arange(len(sorted_subfamilies))

    # Top: mean TE length
    lengths = [mean_length.get(sf, 0) for sf in sorted_subfamilies]
    ax_top.plot(x_pos, lengths, color="black", linewidth=1.5)
    ax_top.set_ylabel("Mean TE length (bp)", fontsize=11)
    ax_top.set_xlim(-0.5, len(sorted_subfamilies) - 0.5)
    ax_top.grid(axis="y", alpha=0.3)

    # Bottom: bar chart log2(IP/Input)
    values = [log2fc.get(sf, 0) for sf in sorted_subfamilies]
    colors = ["#DD8452" if v > 0 else "#4C72B0" for v in values]
    ax_bot.bar(x_pos, values, color=colors, width=0.6, edgecolor="black", linewidth=0.3)
    ax_bot.axhline(0, color="black", linewidth=0.8)

    ax_bot.set_xlabel("TE subfamily", fontsize=11)
    ax_bot.set_ylabel("log2(IP / Input)", fontsize=11)
    ax_bot.set_xticks(x_pos)
    ax_bot.set_xticklabels(sorted_subfamilies, rotation=45, ha="right", fontsize=8)
    ax_bot.grid(axis="y", alpha=0.3)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#DD8452", edgecolor="black", label="IP enriched"),
        Patch(facecolor="#4C72B0", edgecolor="black", label="Input enriched"),
    ]
    ax_bot.legend(handles=legend_elements, fontsize=9, loc="upper right")

    fig.suptitle(f"{ip_name} vs {input_name}", fontsize=13, fontweight="bold", y=0.98)
    plt.tight_layout()
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Chart saved: {output_path}")


def plot_combined(df_all, pairs, output_path, method="count", sort_by="te_length", top_n=30):
    """Plot enrichment for all pairs combined in one figure."""
    n_pairs = len(pairs)
    if n_pairs == 0:
        logger.warning("No pairs to plot")
        return

    fig, axes = plt.subplots(
        n_pairs * 2, 1,
        figsize=(max(14, top_n * 0.5), 4 * n_pairs * 2),
        gridspec_kw={"height_ratios": [1, 2] * n_pairs},
    )
    if n_pairs == 1:
        axes = [axes]

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#DD8452", edgecolor="black", label="IP enriched"),
        Patch(facecolor="#4C72B0", edgecolor="black", label="Input enriched"),
    ]

    for idx, (ip_name, input_name) in enumerate(pairs.items()):
        ax_top = axes[idx * 2]
        ax_bot = axes[idx * 2 + 1]

        # Filter data for this pair
        pair_samples = [ip_name, input_name]
        df_pair = df_all[df_all["sample_id"].isin(pair_samples)].copy()
        df_pair["condition"] = df_pair["sample_id"].apply(
            lambda x: "IP" if x == ip_name else "Input"
        )

        log2fc = compute_enrichment_single_pair(df_pair, method)
        if log2fc.empty:
            continue

        sorted_subfamilies = select_top_subfamilies(log2fc, df_pair, top_n)
        log2fc = log2fc.loc[sorted_subfamilies]

        mean_length = df_pair.groupby("te_subfamily")["te_length"].mean()
        if sort_by == "te_length":
            common = log2fc.index.intersection(mean_length.index)
            sorted_subfamilies = mean_length.loc[common].sort_values().index.tolist()
            log2fc = log2fc.reindex(sorted_subfamilies)

        mean_length = mean_length.reindex(sorted_subfamilies)
        x_pos = np.arange(len(sorted_subfamilies))

        # Top
        lengths = [mean_length.get(sf, 0) for sf in sorted_subfamilies]
        ax_top.plot(x_pos, lengths, color="black", linewidth=1.5)
        ax_top.set_ylabel("Mean TE length (bp)", fontsize=9)
        ax_top.set_xlim(-0.5, len(sorted_subfamilies) - 0.5)
        ax_top.grid(axis="y", alpha=0.3)
        ax_top.set_title(f"{ip_name} vs {input_name}", fontsize=11, fontweight="bold")

        # Bottom
        values = [log2fc.get(sf, 0) for sf in sorted_subfamilies]
        colors = ["#DD8452" if v > 0 else "#4C72B0" for v in values]
        ax_bot.bar(x_pos, values, color=colors, width=0.6, edgecolor="black", linewidth=0.3)
        ax_bot.axhline(0, color="black", linewidth=0.8)
        ax_bot.set_ylabel("log2(IP / Input)", fontsize=9)
        ax_bot.set_xticks(x_pos)
        ax_bot.set_xticklabels(sorted_subfamilies, rotation=45, ha="right", fontsize=7)
        ax_bot.grid(axis="y", alpha=0.3)
        if idx == 0:
            ax_bot.legend(handles=legend_elements, fontsize=8, loc="upper right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Combined chart saved: {output_path}")


def main():
    p = argparse.ArgumentParser("Plot TE subfamily enrichment chart")
    p.add_argument("--tsvs", required=True, help="Comma-separated per-sample TSV paths")
    p.add_argument("--output", required=True, help="Output directory (for separate mode) or file (for combined mode)")
    p.add_argument("--combined-tsv", required=True, help="Output combined TSV path")
    p.add_argument("--ip-input-pairs", default="", help="IP:Input pairs, semicolon-separated")
    p.add_argument("--method", default="reads", choices=["interval", "count", "reads"],
                   help="Enrichment method: 'reads' (IP/Input read ratio, default), 'count' (peak count), 'interval' (overlap fraction)")
    p.add_argument("--sort-by", default="te_length", choices=["enrichment", "te_length"],
                   help="Sort subfamilies by 'te_length' (default) or 'enrichment'")
    p.add_argument("--top-n", type=int, default=30, help="Number of top subfamilies to show (default: 30)")
    p.add_argument("--combine", action="store_true",
                   help="Generate a single combined plot instead of separate per-pair plots")
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

    if args.combine:
        # Single combined plot
        plot_combined(df, pairs, args.output,
                      method=args.method, sort_by=args.sort_by, top_n=args.top_n)
    else:
        # Separate plot per pair
        output_dir = args.output
        os.makedirs(output_dir, exist_ok=True)
        for ip_name, input_name in pairs.items():
            pair_samples = [ip_name, input_name]
            df_pair = df[df["sample_id"].isin(pair_samples)].copy()
            df_pair["condition"] = df_pair["sample_id"].apply(
                lambda x: "IP" if x == ip_name else "Input"
            )
            out_path = os.path.join(output_dir, f"{ip_name}_vs_{input_name}_enrichment.png")
            plot_single_pair(df_pair, pairs, ip_name, input_name, out_path,
                             method=args.method, sort_by=args.sort_by, top_n=args.top_n)


if __name__ == "__main__":
    main()
