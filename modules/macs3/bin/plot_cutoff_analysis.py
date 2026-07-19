#!/usr/bin/env python3
"""Plot peak count trends across different p/q-score cutoffs for each sample."""

import argparse
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_cutoff_table(path: str) -> dict[str, list[float]]:
    """Parse cutoff_analysis.txt into column dict."""
    pscores, qscores, npeaks = [], [], []
    with open(path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            pscores.append(float(parts[0]))
            qscores.append(float(parts[1]))
            npeaks.append(float(parts[2]))
    return {"pscore": pscores, "qscore": qscores, "npeaks": npeaks}


def main():
    parser = argparse.ArgumentParser(description="Plot peak count vs score cutoffs")
    parser.add_argument("--input-files", action="append", required=True,
                        help="cutoff_analysis.txt paths (one per sample)")
    parser.add_argument("--sample-names", action="append", required=True,
                        help="Sample names corresponding to input files")
    parser.add_argument("--output", required=True, help="Output figure path")
    args = parser.parse_args()

    if len(args.input_files) != len(args.sample_names):
        print("ERROR: --input-files and --sample-names count mismatch", file=sys.stderr)
        sys.exit(1)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

    colors = ["#E64B35", "#4DBBD5", "#00A087", "#F39B7F", "#8491B4"]

    for i, (fpath, name) in enumerate(zip(args.input_files, args.sample_names)):
        data = load_cutoff_table(fpath)
        c = colors[i % len(colors)]

        axes[0].plot(data["pscore"], data["npeaks"], marker="o", ms=3,
                     label=name, color=c, linewidth=1.5)
        axes[1].plot(data["qscore"], data["npeaks"], marker="o", ms=3,
                     label=name, color=c, linewidth=1.5)

    axes[0].set_xlabel("-log10(p-value) cutoff", fontsize=12)
    axes[0].set_ylabel("Number of peaks", fontsize=12)
    axes[0].set_title("Peak count vs p-score cutoff", fontsize=13)
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=10)
    axes[0].invert_xaxis()
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("-log10(q-value) cutoff", fontsize=12)
    axes[1].set_ylabel("Number of peaks", fontsize=12)
    axes[1].set_title("Peak count vs q-score cutoff", fontsize=13)
    axes[1].set_yscale("log")
    axes[1].legend(fontsize=10)
    axes[1].invert_xaxis()
    axes[1].grid(True, alpha=0.3)

    fig.suptitle("MACS3 Peak Count Trends Across Score Cutoffs", fontsize=14, fontweight="bold")
    fig.savefig(args.output, dpi=200)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
