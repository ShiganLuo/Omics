#!/usr/bin/env python3
"""Multi-sample combined telomere figures (all samples in one plot).

Per-method combined figure: per-read value distribution (log-scale density
histogram + violin/box) across all samples in ONE figure. With ``--design``
(the flattened ``comparisons`` module config key), samples are ordered
control-first and annotated with their comparison group + role, colored by
group. Also renders a cross-method overview (per-method per-sample summary).

Usage:
    # per-method combined figure (length column from each sample's TSV)
    python plot_combined.py -m read_density -s youth=youth_read_telomere.tsv \
        -s aged=aged_read_telomere.tsv --col tel_bp_estimate -o out/ \
        --design '{"youth_vs_aged": {"control_group_name": "youth",
        "experimental_group_name": "aged", "control_samples": ["s1"],
        "experimental_samples": ["s2"]}}'
    # cross-method overview (each method's summary TSV -> one bar/summary fig)
    python plot_combined.py -m overview -s read_density=...summary.tsv \
        -s assembly_scan=...stats.tsv -o out/ -f png
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    from .LogUtil import setup_logger
    from .plot_util import normalize_formats, save_fig
except ImportError:
    from LogUtil import setup_logger
    from plot_util import normalize_formats, save_fig

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

logger = setup_logger("plot_combined", stream=sys.stdout)

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#000000"]

_RC: Dict[str, Any] = {
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9.5,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "-",
    "axes.axisbelow": True,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
}


def _load_column(path: Path, col: str) -> List[float]:
    """Read numeric values of `col` from a TSV (header row)."""
    vals: List[float] = []
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or col not in header:
            return vals
        j = header.index(col)
        for row in reader:
            if j < len(row):
                try:
                    vals.append(float(row[j]))
                except ValueError:
                    pass
    return vals


def _design_sample_map(design: Optional[Dict]) -> Tuple[List[str], Dict[str, str], Dict[str, str]]:
    """Flattened comparisons design -> (control-first sample order,
    sample->group name, sample->role) maps; first occurrence wins."""
    order: List[str] = []
    group_of: Dict[str, str] = {}
    role_of: Dict[str, str] = {}
    for comp in (design or {}).values():
        if not isinstance(comp, dict):
            continue
        cg = comp.get("control_group_name") or "control"
        eg = comp.get("experimental_group_name") or "experiment"
        for s in comp.get("control_samples") or []:
            if s not in group_of:
                group_of[s], role_of[s] = cg, "control"
                order.append(s)
        for s in comp.get("experimental_samples") or []:
            if s not in group_of:
                group_of[s], role_of[s] = eg, "experiment"
                order.append(s)
    return order, group_of, role_of


def _order_by_design(pairs: List[Tuple[str, str]], design: Optional[Dict]
                     ) -> List[Tuple[str, str]]:
    """Stable-sort (label, tsv) pairs: design control samples first (leftmost),
    then experiment samples, then samples absent from the design."""
    if not design:
        return pairs
    order, _, _ = _design_sample_map(design)
    rank = {s: i for i, s in enumerate(order)}
    return [p for _, p in sorted(
        enumerate(pairs), key=lambda t: (rank.get(t[1][0], len(rank)), t[0]))]


def _group_colors(labels: Sequence[str], group_of: Dict[str, str]
                  ) -> Tuple[List[str], Dict[str, str]]:
    """One palette color per design group (fallback: per sample)."""
    colors: List[str] = []
    color_of: Dict[str, str] = {}
    for lab in labels:
        key = group_of.get(lab, lab)
        if key not in color_of:
            color_of[key] = PALETTE[len(color_of) % len(PALETTE)]
        colors.append(color_of[key])
    return colors, color_of


def _style(ax, ylog: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", visible=False)
    if ylog:
        ax.set_yscale("log")
        ax.grid(axis="y", which="both", alpha=0.2)


def plot_per_method(pairs: List[Tuple[str, str]], col: str, method: str,
                    out_base: str, formats: Optional[Sequence[str]],
                    design: Optional[Dict] = None) -> List[str]:
    """Combined figure: density histogram + violin/box of `col` across samples."""
    pairs = _order_by_design(pairs, design)
    order, group_of, role_of = _design_sample_map(design)
    data: List[List[float]] = []
    labels: List[str] = []
    for label, tsv in pairs:
        vals = _load_column(Path(tsv), col)
        if vals:
            data.append(vals)
            labels.append(label)
    if not data:
        logger.warning(f"No numeric values for column '{col}'; nothing to plot.")
        return []
    colors, _ = _group_colors(labels, group_of)

    all_pos = all(min(v) > 0 for v in data)
    xlog = all_pos
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    with plt.rc_context(_RC):
        ax = axes[0]
        for vals, lab, color in zip(data, labels, colors):
            if xlog:
                lo = math.log10(min(min(v) for v in data))
                hi = math.log10(max(max(v) for v in data))
                bins = np.logspace(lo, hi, 60)
            else:
                bins = 60
            ax.hist(vals, bins=bins, density=True, histtype="stepfilled",
                    alpha=0.35, color=color, edgecolor=color, linewidth=1.2)
        _style(ax, ylog=False)
        if xlog:
            ax.set_xscale("log")
            ax.grid(axis="x", which="both", alpha=0.2)
        ax.set_xlabel(f"{method} · {col} per read" + (" (log scale)" if xlog else ""))
        ax.set_ylabel("density")

        ax = axes[1]
        vp = ax.violinplot(data, showmedians=False, showextrema=False, widths=0.85)
        for body, color in zip(vp["bodies"], colors):
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.55)
        ax.boxplot(data, widths=0.16, showfliers=False, zorder=3,
                   boxprops=dict(color="0.25", linewidth=1.0),
                   whiskerprops=dict(color="0.25", linewidth=1.0),
                   capprops=dict(color="0.25", linewidth=1.0),
                   medianprops=dict(color="black", linewidth=1.4))
        for i, (vals, color) in enumerate(zip(data, colors), start=1):
            med = sorted(vals)[len(vals) // 2]
            ax.annotate(f"{med:,.0f}", (i, med), xytext=(0, 6),
                        textcoords="offset points", ha="center",
                        fontsize=9, color=color, fontweight="bold")
        _style(ax, ylog=all_pos)
        ax.set_ylabel(f"{col} per read" + (" (log scale)" if all_pos else ""))
        ax.set_xticks(range(1, len(labels) + 1))
        ticklabels = [lab if lab not in group_of
                      else f"{lab}\n{group_of[lab]} · {role_of[lab]}"
                      for lab in labels]
        ax.set_xticklabels(ticklabels)
        for tick, color in zip(ax.get_xticklabels(), colors):
            tick.set_color(color)

        legend_labels = [lab if lab not in group_of else f"{lab} ({group_of[lab]})"
                         for lab in labels]
        fig.legend(handles=[Patch(facecolor=c, edgecolor=c, alpha=0.55, label=l)
                            for c, l in zip(colors, legend_labels)],
                   loc="upper center", ncol=max(1, len(labels)), frameon=False,
                   bbox_to_anchor=(0.5, 1.0))
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def plot_overview(pairs: List[Tuple[str, str]], out_base: str,
                  formats: Optional[Sequence[str]]) -> List[str]:
    """Cross-method overview: one summary value per method per sample."""
    # Read each method's summary TSV: first numeric column's mean as the summary.
    methods = []
    vals = []
    for method, tsv in pairs:
        nums: List[float] = []
        with Path(tsv).open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            if not header:
                continue
            # pick the first numeric column
            j = None
            for k in range(len(header)):
                col = [r[k] for r in reader if k < len(r)]
                fh.seek(0); next(csv.reader(fh, delimiter="\t"))  # reset
                if col and sum(_isnum(x) for x in col) / len(col) > 0.8:
                    j = k
                    break
            if j is None:
                continue
            fh.seek(0); reader = csv.reader(fh, delimiter="\t"); next(reader, None)
            for row in reader:
                if j < len(row) and _isnum(row[j]):
                    nums.append(float(row[j]))
        if nums:
            methods.append(method)
            vals.append(sum(nums) / len(nums))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    with plt.rc_context(_RC):
        bars = ax.bar(methods, vals,
                      color=[PALETTE[i % len(PALETTE)] for i in range(len(vals))],
                      width=0.55)
        for bar, v in zip(bars, vals):
            ax.annotate(f"{v:,.1f}", (bar.get_x() + bar.get_width() / 2, v),
                        xytext=(0, 4), textcoords="offset points",
                        ha="center", fontsize=9)
        _style(ax)
        ax.set_ylabel("mean value")
        ax.set_ylim(0, max(vals) * 1.15 if vals else 1)
        fig.tight_layout()
        written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def _isnum(x: str) -> bool:
    try:
        float(x); return True
    except ValueError:
        return False


def main() -> None:
    """CLI: render per-method combined figure or cross-method overview."""
    p = argparse.ArgumentParser(description="Multi-sample combined telomere figures.")
    p.add_argument("-m", "--method", required=True,
                   help="method name (or 'overview')")
    p.add_argument("-s", "--sample", action="append", required=True,
                   metavar="LABEL=TSV", help="label + per-sample data TSV")
    p.add_argument("-c", "--col", default="tel_bp_estimate",
                   help="numeric column to plot (per-method mode)")
    p.add_argument("-D", "--design", default="",
                   help="flattened comparisons design JSON "
                        "({name: {control_group_name, experimental_group_name, "
                        "control_samples, experimental_samples}}); orders samples "
                        "control-first and annotates groups")
    p.add_argument("-o", "--outdir", default="results")
    p.add_argument("-b", "--out_base", default="")
    p.add_argument("-F", "--format", action="append", dest="formats", metavar="FMT",
                   help="output figure format(s), repeatable (default png)")
    p.add_argument("-l", "--log", default="")
    args = p.parse_args()

    log = logger
    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        log = setup_logger("plot_combined", log_file=args.log, stream=sys.stdout)

    design = None
    if args.design:
        try:
            design = json.loads(args.design)
        except json.JSONDecodeError as e:
            raise SystemExit(
                f"--design is not valid JSON ({e}); pass the flattened "
                f"comparisons dict as produced by the subworkflow.")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    pairs = []
    for spec in args.sample:
        label, _, tsv = spec.partition("=")
        pairs.append((label, tsv))
    base = args.out_base or (f"combined_{args.method}" if args.method != "overview"
                             else "combined_overview")
    out_base = str(out / base)
    if args.method == "overview":
        written = plot_overview(pairs, out_base, args.formats)
    else:
        written = plot_per_method(pairs, args.col, args.method, out_base,
                                  args.formats, design)
    for w in written:
        log.info(f"Wrote: {w}")


if __name__ == "__main__":
    main()
