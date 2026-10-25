#!/usr/bin/env python3
"""Statistics and visualization for Decoil ecDNA reconstruction results.

Reads per-sample Decoil outputs (``summary.txt`` plus
``reconstruct.ecDNA.filtered.bed``) and produces:

- per-sample: cycle size distribution histogram, topology composition bar,
  ecDNA region table (TSV)
- cross-sample (``--mode compare``): circle/ecDNA counts bar, size
  distribution overlay, topology stacked bar, shared/unique ecDNA region
  table (TSV)

All matplotlib text is English (CJK font compatibility, project convention).

CLI conventions follow the project's plotting tools: ``LogUtil``-style
logger is NOT required here (standalone bin/ script), logging goes to
stderr; image formats are repeatable via ``-f``.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Dict, List, Optional

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
except ImportError as e:  # pragma: no cover
    sys.stderr.write(
        f"Missing dependency: {e}\nInstall matplotlib and pandas in the "
        "execution environment (they are provided by the decoil/track conda env).\n"
    )
    raise

PlotFormat = str  # png/pdf/svg/...; validated by matplotlib at save time

TOP_COLORS = {
    "simple_circle": "#3498DB",
    "simple_sv": "#9B59B6",
    "multi_region_intra_chr": "#E67E22",
    "foldbacks": "#E74C3C",
    "multi_region_inter_chr": "#16A085",
}

DPI = 300


def setup_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("decoil_plot")


logger = setup_logging()


def load_summary(path: str, sample_id: str) -> pd.DataFrame:
    """Load one decoil summary.txt into a DataFrame with sample column.

    Parameters
    ----------
    path : str
        Path to ``summary.txt`` (tab-separated, has header).
    sample_id : str
        Sample label attached to every row.

    Returns
    -------
    pandas.DataFrame
        Columns: sample_id, circ_id, chr_origin, size_mb, label, topology.

    Raises
    ------
    ValueError
        If the file is missing required columns.
    """
    try:
        df = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        # Empty summary.txt: no ecDNA reconstructed for this sample
        return pd.DataFrame(columns=["sample_id", "circ_id", "chr_origin",
                                     "size_mb", "label", "topology"])
    df.columns = [str(c).strip() for c in df.columns]
    required = {"circ_id", "chr_origin", "size(MB)", "topology_name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    out = pd.DataFrame({
        "sample_id": sample_id,
        "circ_id": df["circ_id"],
        "chr_origin": df["chr_origin"],
        "size_mb": pd.to_numeric(df["size(MB)"], errors="coerce"),
        "label": df.get("label", "").fillna(""),
        "topology": df["topology_name"],
    })
    return out.dropna(subset=["size_mb"])


def load_ecdna_bed(path: str, sample_id: str) -> pd.DataFrame:
    """Load reconstruct.ecDNA.filtered.bed (comment-line tolerant).

    Parameters
    ----------
    path : str
        Path to the filtered ecDNA BED.
    sample_id : str
        Sample label attached to every row.

    Returns
    -------
    pandas.DataFrame
        Columns: sample_id, chrom, start, end, circ_id, coverage, proportion.
        Empty frame when the file has no data lines.
    """
    try:
        df = pd.read_csv(path, sep="\t", comment="#", header=None,
                         dtype=str)
    except pd.errors.EmptyDataError:
        # All-comment / empty file: no ecDNA passed the filter
        df = pd.DataFrame()
    if df.empty:
        return pd.DataFrame(columns=["sample_id", "chrom", "start", "end",
                                     "circ_id", "coverage", "proportion"])
    # Expected: chr start end circ_id fragment_id strand coverage proportions
    df = df.iloc[:, [0, 1, 2, 3, 6, 7]] if df.shape[1] >= 8 else df
    df.columns = ["chrom", "start", "end", "circ_id", "coverage", "proportion"]
    df.insert(0, "sample_id", sample_id)
    for col in ("start", "end"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def merge_ecdna_fragments(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse BED fragment rows belonging to one circ_id into one region.

    Parameters
    ----------
    df : pandas.DataFrame
        Fragment-level ecDNA BED rows (see :func:`load_ecdna_bed`).

    Returns
    -------
    pandas.DataFrame
        One row per (sample_id, circ_id) with min(start), max(end) and the
        max coverage/proportion across fragments, plus full genomic
        coordinates in ``chr:start-end`` format (full precision, project
        convention — never abbreviated).
    """
    if df.empty:
        return df
    agg = (df.groupby(["sample_id", "circ_id"])
             .agg(chrom=("chrom", "first"),
                  start=("start", "min"),
                  end=("end", "max"),
                  coverage=("coverage", "max"),
                  proportion=("proportion", "max"))
             .reset_index())
    agg["region"] = (
        agg["chrom"] + ":" + agg["start"].astype("Int64").astype(str)
        + "-" + agg["end"].astype("Int64").astype(str)
    )
    return agg


def save_fig(fig: "plt.Figure", out_prefix: str, image_formats: List[str]) -> None:
    """Save a figure in every requested format."""
    for fmt in image_formats:
        path = f"{out_prefix}.{fmt}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        logger.info(f"Saved {path}")


def plot_per_sample(summary: pd.DataFrame, ecdna: pd.DataFrame,
                    sample_id: str, out_prefix: str,
                    image_formats: List[str]) -> None:
    """Size histogram + topology composition + ecDNA table for one sample.

    Parameters
    ----------
    summary : pandas.DataFrame
        All reconstructed cycles (see :func:`load_summary`).
    ecdna : pandas.DataFrame
        Merged ecDNA regions (see :func:`merge_ecdna_fragments`).
    sample_id : str
        Sample label used in titles.
    out_prefix : str
        Output path prefix; ``.<fmt>`` appended per format.
    image_formats : list of str
        Image formats (png/pdf/...).
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    sizes_kb = summary["size_mb"] * 1000.0
    import numpy as np
    bins = 10.0 ** np.linspace(0, 3, 25)  # 1 kb .. 1 Mb log bins
    ax.hist(sizes_kb.clip(lower=0.5), bins=bins, color="#3498DB",
            edgecolor="white", linewidth=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("Circle size (kb, log scale)")
    ax.set_ylabel("Number of cycles")
    ax.set_title(f"{sample_id}: cycle size distribution (n={len(summary)})")

    ax = axes[1]
    topo_counts = summary["topology"].value_counts()
    colors = [TOP_COLORS.get(t, "#95A5A6") for t in topo_counts.index]
    ax.bar(range(len(topo_counts)), topo_counts.values, color=colors)
    ax.set_xticks(range(len(topo_counts)))
    ax.set_xticklabels(topo_counts.index, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Number of cycles")
    ax.set_title(f"{sample_id}: topology composition")

    fig.tight_layout()
    save_fig(fig, out_prefix, image_formats)
    plt.close(fig)

    # ecDNA region table (full coordinates)
    table_path = f"{out_prefix}_ecdna_regions.tsv"
    ecdna_out = ecdna[["sample_id", "region", "circ_id", "coverage",
                       "proportion"]] if not ecdna.empty else ecdna
    ecdna_out.to_csv(table_path, sep="\t", index=False)
    logger.info(f"Saved {table_path} ({len(ecdna_out)} ecDNA regions)")


def plot_compare(all_summary: pd.DataFrame, all_ecdna: pd.DataFrame,
                 sample_order: List[str], out_prefix: str,
                 image_formats: List[str]) -> None:
    """Cross-sample comparison panel.

    Parameters
    ----------
    all_summary : pandas.DataFrame
        Concatenated per-sample summaries.
    all_ecdna : pandas.DataFrame
        Concatenated merged ecDNA regions.
    sample_order : list of str
        Sample display order (x-axis / legend order).
    out_prefix : str
        Output path prefix.
    image_formats : list of str
        Image formats (png/pdf/...).
    """
    import numpy as np

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    import numpy as np
    x = np.arange(len(sample_order))
    width = 0.38

    # Panel 1: cycle / ecDNA counts
    ax = axes[0]
    total = [int((all_summary["sample_id"] == s).sum()) for s in sample_order]
    ecdna_n = [int((all_ecdna["sample_id"] == s).sum()) for s in sample_order]
    ax.bar(x - width / 2, total, width, label="all cycles", color="#3498DB")
    ax.bar(x + width / 2, ecdna_n, width, label="ecDNA (filtered)",
           color="#E74C3C")
    for xi, v in zip(x - width / 2, total):
        ax.text(xi, v, str(v), ha="center", va="bottom", fontsize=8)
    for xi, v in zip(x + width / 2, ecdna_n):
        ax.text(xi, v, str(v), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order)
    ax.set_ylabel("Count")
    ax.set_title("Reconstructed cycles per sample")
    ax.legend(fontsize=8)

    # Panel 2: size distribution overlay
    ax = axes[1]
    bins = 10.0 ** np.linspace(0, 3, 25)
    cmap = plt.cm.tab10
    for i, s in enumerate(sample_order):
        sizes_kb = all_summary.loc[all_summary["sample_id"] == s,
                                   "size_mb"] * 1000.0
        ax.hist(sizes_kb.clip(lower=0.5), bins=bins, histtype="step",
                linewidth=1.6, color=cmap(i / max(len(sample_order) - 1, 1)),
                label=s)
    ax.set_xscale("log")
    ax.set_xlabel("Circle size (kb, log scale)")
    ax.set_ylabel("Number of cycles")
    ax.set_title("Cycle size distribution")
    ax.legend(fontsize=8)

    # Panel 3: topology composition (stacked, per sample)
    ax = axes[2]
    topologies = sorted(all_summary["topology"].unique())
    bottom = np.zeros(len(sample_order))
    for topo in topologies:
        vals = np.array([
            int(((all_summary["sample_id"] == s)
                 & (all_summary["topology"] == topo)).sum())
            for s in sample_order
        ])
        ax.bar(x, vals, bottom=bottom, label=topo,
               color=TOP_COLORS.get(topo, "#95A5A6"))
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order)
    ax.set_ylabel("Number of cycles")
    ax.set_title("Topology composition")
    ax.legend(fontsize=7, loc="upper right")

    fig.tight_layout()
    save_fig(fig, out_prefix, image_formats)
    plt.close(fig)

    # Shared / unique ecDNA regions (keyed by full coordinate region string)
    regions = all_ecdna[["sample_id", "region", "circ_id", "coverage",
                         "proportion"]]
    table_path = f"{out_prefix}_ecdna_regions.tsv"
    regions.to_csv(table_path, sep="\t", index=False)
    logger.info(f"Saved {table_path}")

    if len(sample_order) >= 2:
        sets = {s: set(all_ecdna.loc[all_ecdna["sample_id"] == s, "region"])
                for s in sample_order}
        shared = set.intersection(*sets.values()) if sets else set()
        union = set.union(*sets.values()) if sets else set()
        rows = []
        for region in sorted(union):
            status = "shared" if region in shared else "unique"
            per_sample = {
                s: int(region in sets[s]) for s in sample_order
            }
            row = {"region": region, "status": status}
            row.update(per_sample)
            rows.append(row)
        shared_path = f"{out_prefix}_shared_unique.tsv"
        pd.DataFrame(rows).to_csv(shared_path, sep="\t", index=False)
        logger.info(f"Saved {shared_path} ({len(shared)} shared / "
                    f"{len(union)} total ecDNA regions)")


def parse_sample_args(pairs: List[str]) -> List:
    """Parse repeated ``--sample SAMPLE:PATH`` args.

    Parameters
    ----------
    pairs : list of str
        Each entry ``SAMPLE:PATH`` where PATH may itself contain colons
        (Windows-style paths are tolerated via ``split(":", 1)``).

    Returns
    -------
    list of tuple(str, str)
        (sample_id, summary_path) pairs.

    Raises
    ------
    SystemExit
        If any entry lacks the ``SAMPLE:PATH`` separator.
    """
    parsed = []
    for item in pairs:
        if ":" not in item:
            sys.exit(f"Invalid --sample entry '{item}', expected SAMPLE:PATH")
        name, path = item.split(":", 1)
        if not os.path.isfile(path):
            sys.exit(f"Summary file not found: {path}")
        parsed.append((name, path))
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Statistics and plots for Decoil ecDNA results.")
    parser.add_argument("--mode", choices=("single", "compare"),
                        default="single",
                        help="single: one sample panel; compare: cross-sample")
    parser.add_argument("--sample", action="append", required=True,
                        metavar="SAMPLE:SUMMARY_TSV",
                        help="sample id + path to summary.txt; repeatable")
    parser.add_argument("--ecdna-bed", action="append", default=None,
                        metavar="SAMPLE:BED",
                        help="sample id + reconstruct.ecDNA.filtered.bed; "
                             "repeatable, same order as --sample")
    parser.add_argument("-o", "--out-prefix", required=True,
                        help="output path prefix (extensions appended)")
    parser.add_argument("-f", "--format", action="append", dest="formats",
                        metavar="FMT", help="image format, repeatable; "
                                            "default: png")
    args = parser.parse_args()

    image_formats = args.formats or ["png"]

    pairs = parse_sample_args(args.sample)
    bed_pairs = []
    if args.ecdna_bed:
        bed_pairs = parse_sample_args(args.ecdna_bed)
        if len(bed_pairs) != len(pairs):
            sys.exit("--ecdna-bed count must match --sample count")

    summaries, ecdnas = [], []
    for i, (sid, summary_path) in enumerate(pairs):
        summaries.append(load_summary(summary_path, sid))
        bed_path = bed_pairs[i][1] if bed_pairs else None
        if bed_path:
            ecdnas.append(merge_ecdna_fragments(load_ecdna_bed(bed_path, sid)))
        else:
            ecdnas.append(pd.DataFrame())
    all_summary = pd.concat(summaries, ignore_index=True)
    all_ecdna = (pd.concat(ecdnas, ignore_index=True)
                 if ecdnas and not all(d.empty for d in ecdnas)
                 else pd.DataFrame(columns=["sample_id", "region", "circ_id",
                                            "coverage", "proportion"]))

    sample_order = [s for s, _ in pairs]
    if args.mode == "single":
        sid, _ = pairs[0]
        plot_per_sample(all_summary, all_ecdna, sid, args.out_prefix,
                        image_formats)
    else:
        plot_compare(all_summary, all_ecdna, sample_order, args.out_prefix,
                     image_formats)


if __name__ == "__main__":
    main()
