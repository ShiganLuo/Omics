#!/usr/bin/env python3
"""fig9-style plot (title-free): epigenetic drift from the computed windows.

Port of conclusion_figures/plot_fig9_replot.py to a reusable module script:
group labels come from the summary JSON instead of hardcoded aged/youth.

Left: histogram of log2 variance ratio (B / A) of the per-fiber m5CpG
fraction across windows, with the median marked.  Right: drift vs mean
shift scatter.

Usage:
  python epi_drift_plot.py \
      --windows epigenetic_drift_windows.tsv \
      --summary epigenetic_drift.json \
      --out-base out/epigenetic_drift --format png --dpi 200
"""
import argparse
import csv
import json
import math
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out-base", required=True)
    ap.add_argument("--format", action="append", default=None)
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()
    fmts = args.format or ["png"]

    meta = json.load(open(args.summary))
    label_a, label_b = meta["groups"]
    n, signs = meta["n"], meta["signs"]
    k = min(signs, n - signs)
    p = 2 * sum(math.comb(n, i) for i in range(0, k + 1)) / 2 ** n
    p = min(1.0, p)

    rows = list(csv.DictReader(open(args.windows), delimiter="\t"))
    var_ratios = [(float(r[f"var_{label_b}"]) + 1e-12)
                  / (float(r[f"var_{label_a}"]) + 1e-12) for r in rows]
    lr = [math.log2(v) for v in var_ratios]

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    axes[0].hist(lr, bins=40, color="#2E86C1", alpha=0.75)
    axes[0].axvline(0, color="black", lw=1)
    axes[0].axvline(statistics.median(lr), color="#C0392B", lw=1.4, ls="--",
                    label=f"median = {2 ** statistics.median(lr):.3f}")
    axes[0].set_xlabel(f"log2 variance ratio ({label_b} / {label_a})")
    axes[0].set_ylabel("windows")
    axes[0].legend(fontsize=8)
    md = [float(r[f"mean_{label_b}"]) - float(r[f"mean_{label_a}"])
          for r in rows]
    axes[1].scatter(md, lr, s=7, alpha=0.4, color="#7F8C8D")
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].axvline(0, color="black", lw=0.8)
    axes[1].set_xlabel(f"mean difference ({label_b} - {label_a})")
    axes[1].set_ylabel("log2 variance ratio")
    axes[1].locator_params(axis="x", nbins=5)
    axes[1].ticklabel_format(axis="x", style="sci", scilimits=(0, 0),
                             useOffset=False)
    fig.tight_layout()
    for fmt in fmts:
        out = f"{args.out_base}.{fmt}"
        fig.savefig(out, dpi=args.dpi)
        print("saved", out)
    print(f"n={n} {label_b}>{label_a} in {signs}/{n} windows "
          f"(sign-test p={p:.3g}) median_var_ratio="
          f"{statistics.median(var_ratios):.3f}")


if __name__ == "__main__":
    main()
