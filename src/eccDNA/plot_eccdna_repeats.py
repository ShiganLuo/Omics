#!/usr/bin/env python3
"""Summary figure for eccDNA repeat-content analyses.

Reads the TSV outputs of `breakpoint_repeat.py`, `ecdna_repeat_enrichment.py`
and `junction_microhomology.py` from --out-dir and renders a 2x2 panel figure:

- A  breakpoint repeat-hit rate vs window (obs vs three nulls)
- B  repeat families flanking breakpoints (enrichment vs local null)
- C  repeat families on eccDNA bodies (matched-background enrichment)
- D  junction microhomology thresholds (obs vs null)

Parameters
----------
--out-dir : directory containing the analysis TSVs
--min-family-bp : minimum eccDNA bp for families drawn in panel C (default 5000)
--top-n : families drawn per bar panel (default 10)
--format : png or pdf (default png)
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import sys

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("ecdna_plot_repeats")

NAVY = "#1F3864"
GOLD = "#C9A227"
SLATE = "#6E7B8B"
RUST = "#A64B2A"
TEAL = "#2E6E6B"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--out-dir", required=True,
                        help="Directory with the analysis TSVs")
    parser.add_argument("--min-family-bp", type=int, default=5000,
                        help="Min eccDNA bp for panel C families (default: 5000)")
    parser.add_argument("--top-n", type=int, default=10,
                        help="Families per bar panel (default: 10)")
    parser.add_argument("--format", choices=("png", "pdf"), default="png",
                        help="Figure format (default: png)")
    return parser.parse_args(argv)


def read_tsv(path: str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv: list[str] | None = None) -> int:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        logger.error(f"matplotlib unavailable ({exc}); install it and rerun")
        return 1
    args = parse_args(argv)
    d = args.out_dir

    win = read_tsv(os.path.join(d, "breakpoint_window_summary.tsv"))
    fam_bp = read_tsv(os.path.join(d, "breakpoint_family_composition.tsv"))
    enrich = read_tsv(os.path.join(d, "ecdna_repeat_enrichment_family.tsv"))
    mh = read_tsv(os.path.join(d, "microhomology_summary.tsv"))

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # ---- A: breakpoint repeat-hit rate vs window ------------------------ #
    any_rows = [r for r in win if r["key"] == "any"]
    windows = sorted({int(r["window_bp"]) for r in any_rows})
    nulls = ["local", "circle_interior", "genome"]
    labels = {"local": "local null (circle span ±50kb)",
              "circle_interior": "circle interior null",
              "genome": "genome null"}
    colors = {"local": GOLD, "circle_interior": SLATE, "genome": TEAL}
    obs = {(int(r["window_bp"])): float(r["obs_frac"]) for r in any_rows
           if r["null"] == "local"}
    ax_a.plot(windows, [obs[w] for w in windows], "-o", color=NAVY,
              lw=2, label="observed")
    for mode in nulls:
        vals = {int(r["window_bp"]): float(r["null_mean_frac"])
                for r in any_rows if r["null"] == mode}
        ax_a.plot(windows, [vals[w] for w in windows], "--s", color=colors[mode],
                  lw=1.5, ms=4, label=labels[mode])
    ax_a.set_xscale("symlog", linthresh=25)
    ax_a.set_xticks(windows)
    ax_a.set_xticklabels([str(w) for w in windows])
    ax_a.set_xlabel("window around breakpoint (bp)")
    ax_a.set_ylabel("fraction of breakpoints with a repeat")
    ax_a.set_title("A  repeat content at eccDNA breakpoints")
    ax_a.legend(fontsize=8)
    ax_a.grid(alpha=0.3)

    # ---- B: families flanking breakpoints ------------------------------- #
    loc_rows = [r for r in fam_bp if r["null"] == "local"
                and float(r["obs_locus_frac"]) >= 0.10]
    loc_rows.sort(key=lambda r: -float(r["enrichment"]))
    loc_rows = loc_rows[: args.top_n]
    loc_rows.sort(key=lambda r: float(r["enrichment"]))
    names = [r["family"] for r in loc_rows]
    enrs = [float(r["enrichment"]) for r in loc_rows]
    bar_colors = [RUST if float(r["p_fdr"]) < 0.05 else NAVY for r in loc_rows]
    ax_b.barh(names, enrs, color=bar_colors)
    ax_b.axvline(1.0, color="black", lw=0.8, ls="--")
    for y, r in enumerate(loc_rows):
        star = ("**" if float(r["p_fdr"]) < 0.05 else
                "*" if float(r["p_emp"]) < 0.05 else "")
        ax_b.text(float(r["enrichment"]) + 0.05, y,
                  f"{float(r['p_emp']):.3f}{star}", va="center", fontsize=7)
    ax_b.set_xlim(0, max(enrs) * 1.45 if enrs else 1.0)
    ax_b.set_xlabel("enrichment vs local null (locus fraction)")
    ax_b.set_title(f"B  repeat families within ±500 bp of breakpoints")
    ax_b.grid(alpha=0.3, axis="x")

    # ---- C: families on eccDNA bodies ---------------------------------- #
    body = [r for r in enrich
            if int(r["ecdna_bp"]) >= args.min_family_bp
            and float(r["enrichment"]) > 0]
    body.sort(key=lambda r: -float(r["enrichment"]))
    body = body[: args.top_n]
    body.sort(key=lambda r: float(r["enrichment"]))
    names = [r["key"] for r in body]
    enrs = [float(r["enrichment"]) for r in body]
    bar_colors = [GOLD if float(r["p_emp"]) < 0.05 else NAVY for r in body]
    ax_c.barh(names, enrs, color=bar_colors)
    ax_c.axvline(1.0, color="black", lw=0.8, ls="--")
    for y, r in enumerate(body):
        star = ("**" if float(r["p_fdr"]) < 0.05 else
                "*" if float(r["p_emp"]) < 0.05 else "")
        ax_c.text(float(r["enrichment"]) + 0.05, y,
                  f"{float(r['p_emp']):.3f}{star}", va="center", fontsize=7)
    ax_c.set_xlim(0, max(enrs) * 1.45 if enrs else 1.0)
    ax_c.set_xlabel("enrichment vs matched random intervals")
    ax_c.set_title(f"C  repeat families on eccDNA bodies "
                   f"(≥{args.min_family_bp // 1000} kb)")
    ax_c.grid(alpha=0.3, axis="x")

    # ---- D: microhomology thresholds ----------------------------------- #
    th = [r["threshold"] for r in mh]
    obs_n = [float(r["obs_n"]) for r in mh]
    null_n = [float(r["null_mean_n"]) for r in mh]
    x = range(len(th))
    w = 0.38
    ax_d.bar([i - w / 2 for i in x], obs_n, width=w, color=NAVY, label="observed")
    ax_d.bar([i + w / 2 for i in x], null_n, width=w, color=GOLD, label="null mean")
    for i, r in enumerate(mh):
        if float(r["p_emp"]) < 0.05:
            ax_d.text(i - w / 2, float(r["obs_n"]) + 0.1,
                      f"p={float(r['p_emp']):.3f}", ha="center", fontsize=7)
    ax_d.set_xticks(list(x))
    ax_d.set_xticklabels(th)
    ax_d.set_xlabel("microhomology length threshold")
    ax_d.set_ylabel("n junctions")
    ax_d.set_title("D  junction microhomology")
    ax_d.legend(fontsize=8)
    ax_d.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    out = os.path.join(d, f"eccdna_repeats_summary.{args.format}")
    fig.savefig(out, dpi=200)
    logger.info(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
