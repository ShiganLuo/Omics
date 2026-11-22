#!/usr/bin/env python3
"""fig2-style plot: telomere/centromere region classes, depth-normalized.

Port of conclusion_figures/plot_fig2_telocen_v5.py to a reusable module
script: group labels come from --group instead of hardcoded youth/aged.

Input: region_normalized.tsv (from telocen_normalized.py).
Two panels: m5CpG per aligned Mb (endogenous methylation) and m6A per
aligned Mb (chromatin accessibility footprint).  Paired Wilcoxon over
regions, bracket stars; stats honesty note: n=1 individual per group ->
direction-consistency across regions, NOT population significance.

Usage:
  python plot_telocen_norm.py \
      --input region_normalized.tsv \
      --group youth=youth --group aged=aged \
      --out-base out/telocen_normalized --format png --dpi 200
"""
import argparse
import csv
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

CLASSES = ["cen_core", "cen_flank", "subtelomere", "telomere_core"]
# first group = baseline (youth colour), second = compared (aged colour)
COLORS = ["#2E86C1", "#C0392B"]


def parse_groups(items):
    groups = {}
    for it in items or []:
        if "=" not in it:
            raise ValueError(f"--group expects LABEL=SAMPLES, got: {it!r}")
        label, s = it.split("=", 1)
        groups[label] = [x for x in s.split(",") if x]
    if len(groups) != 2:
        raise ValueError(f"exactly 2 --group entries required (got {len(groups)})")
    return groups


def stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="region_normalized.tsv")
    ap.add_argument("--group", action="append", required=True,
                    help="LABEL=S1,S2 (repeatable; exactly 2 groups)")
    ap.add_argument("--out-base", required=True, help="extension-less figure base")
    ap.add_argument("--format", action="append", default=None,
                    help="output format (repeatable, default png)")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()
    fmts = args.format or ["png"]

    groups = parse_groups(args.group)
    labels = list(groups)
    label_a, label_b = labels
    sample2label = {}
    for label, slist in groups.items():
        for s in slist:
            sample2label[s] = label

    # (metric, class) -> {label: {region_key: [values]}}
    per = {}
    for r in csv.DictReader(open(args.input), delimiter="\t"):
        key = (r["chrom"], int(r["start"]), int(r["end"]))
        label = sample2label.get(r["sample"])
        if label is None:
            continue
        for metric in ("cpg_per_mb", "m6a_per_mb"):
            if r[metric] == "NA":
                continue
            per.setdefault(metric, {}).setdefault(r["region_class"], {}) \
              .setdefault(label, {}).setdefault(key, []).append(float(r[metric]))

    # keep known class order, append any extra classes found in the data
    classes = [c for c in CLASSES if any(c in per.get(m, {}) for m in per)]
    for m in per:
        for c in per[m]:
            if c not in classes:
                classes.append(c)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.2))
    panels = [("cpg_per_mb", "m5CpG per aligned Mb (endogenous methylation)",
               "m5CpG density, depth-normalized"),
              ("m6a_per_mb", "Hia5-m6A per aligned Mb (= accessibility footprint)",
               "chromatin accessibility, depth-normalized")]

    print_stats = []
    for ax, (metric, ylab, title) in zip(axes, panels):
        ymax = 0.0
        for i, cls in enumerate(classes):
            d = per.get(metric, {}).get(cls, {})
            pa = {k: statistics.mean(v) for k, v in d.get(label_a, {}).items()}
            pb = {k: statistics.mean(v) for k, v in d.get(label_b, {}).items()}
            pairs = [(pa[k], pb[k]) for k in sorted(set(pa) & set(pb))]
            if not pairs:
                continue
            av = [p[0] for p in pairs]
            bv = [p[1] for p in pairs]
            ma, mb = statistics.mean(av), statistics.mean(bv)
            sa = statistics.stdev(av) if len(av) > 1 else 0.0
            sb = statistics.stdev(bv) if len(bv) > 1 else 0.0
            try:
                p = float(wilcoxon(av, bv).pvalue)
            except ValueError:
                p = 1.0
            print_stats.append((metric, cls, len(pairs), ma, mb, mb / ma, p))
            ymax = max(ymax, ma + sa, mb + sb, max(av), max(bv))
            for a1, b1 in pairs:
                ax.plot([i - 0.19, i + 0.19], [a1, b1], color="#B8B8B8",
                        lw=0.5, alpha=0.55, zorder=2)
            ax.scatter([i - 0.19] * len(av), av, s=13, color=COLORS[0],
                       alpha=0.85, zorder=3, linewidths=0)
            ax.scatter([i + 0.19] * len(bv), bv, s=13, color=COLORS[1],
                       alpha=0.85, zorder=3, linewidths=0)
            ax.bar(i - 0.19, ma, 0.34, color=COLORS[0], alpha=0.45, zorder=1,
                   yerr=sa, capsize=3, error_kw=dict(lw=0.9, ecolor="#333333",
                                                     zorder=4))
            ax.bar(i + 0.19, mb, 0.34, color=COLORS[1], alpha=0.45, zorder=1,
                   yerr=sb, capsize=3, error_kw=dict(lw=0.9, ecolor="#333333",
                                                     zorder=4))
            for xoff, val, col in ((i - 0.19, ma + sa, "#1F4E79"),
                                   (i + 0.19, mb + sb, "#7B241C")):
                ax.annotate(f"{val:,.0f}", (xoff, val), ha="center",
                            va="bottom", fontsize=6.5, color=col,
                            xytext=(0, 3), textcoords="offset points",
                            zorder=7,
                            bbox=dict(facecolor="white", edgecolor="none",
                                      alpha=0.75, pad=0.6))
            ax.annotate(f"\u00d7{mb / ma:.2f}", (i, ymax * 0.02), ha="center",
                        va="bottom", fontsize=7.5, color="#333333", zorder=6,
                        bbox=dict(facecolor="white", edgecolor="none",
                                  alpha=0.85, pad=0.8))
        ybar = ymax * 1.15
        for i, cls in enumerate(classes):
            d = per.get(metric, {}).get(cls, {})
            pa = {k for k in d.get(label_a, {})}
            pb = {k for k in d.get(label_b, {})}
            pairs_n = len(pa & pb)
            if not pairs_n:
                continue
            av = [statistics.mean(v) for k, v in sorted(d.get(label_a, {}).items())
                  if k in pb]
            bv = [statistics.mean(v) for k, v in sorted(d.get(label_b, {}).items())
                  if k in pa]
            try:
                p = float(wilcoxon(av, bv).pvalue)
            except ValueError:
                p = 1.0
            ax.plot([i - 0.19, i - 0.19, i + 0.19, i + 0.19],
                    [ybar * 0.97, ybar, ybar, ybar * 0.97], color="#333333",
                    lw=1.0, zorder=5)
            ax.annotate(stars(p), (i, ybar), ha="center", va="bottom",
                        fontsize=10, fontweight="bold", xytext=(0, 2),
                        textcoords="offset points", zorder=5)
        ax.set_ylim(0, ymax * 1.32)
        ax.set_xticks(range(len(classes)))
        ax.set_xticklabels(
            [f"{c}\n(n={len(set(per.get(metric, {}).get(c, {}).get(label_a, {})) & set(per.get(metric, {}).get(c, {}).get(label_b, {})))} regions)"
             for c in classes], fontsize=8)
        ax.set_ylabel(ylab)
        ax.grid(alpha=0.2, axis="y")
        ax.spines[["top", "right"]].set_visible(False)

    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=COLORS[0], alpha=0.45),
                        plt.Rectangle((0, 0), 1, 1, fc=COLORS[1], alpha=0.45)],
               labels=[label_a, label_b], fontsize=8.5, frameon=False,
               loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=2)
    fig.tight_layout(rect=(0, 0.01, 1, 0.955))
    for fmt in fmts:
        out = f"{args.out_base}.{fmt}"
        fig.savefig(out, dpi=args.dpi)
        print("saved", out)
    for m, cls, n, ma, mb, ratio, p in print_stats:
        print(f"  {m:10s} {cls:14s} n={n:2d}  {label_a} {ma:.1f}  "
              f"{label_b} {mb:.1f}  ratio {ratio:.3f}  wilcoxon_p={p:.3g}")


if __name__ == "__main__":
    main()
