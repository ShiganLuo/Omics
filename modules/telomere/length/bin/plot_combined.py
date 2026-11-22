#!/usr/bin/env python3
"""Multi-sample combined telomere figures (all samples in one plot).

Per-method combined figure: overlaid length distribution (histogram + violin)
across all samples, so youth/aged (etc.) compare in ONE plot. Also renders a
cross-method overview (per-method per-sample summary) when the user runs all
methods.

Usage:
    # per-method combined figure (length column from each sample's TSV)
    python plot_combined.py -m read_density -s youth=youth_read_telomere.tsv \
        -s aged=aged_read_telomere.tsv --col tel_bp_estimate -o out/ -f png pdf
    # cross-method overview (each method's summary TSV -> one bar/summary fig)
    python plot_combined.py -m overview -s read_density=...summary.tsv \
        -s assembly_scan=...stats.tsv -o out/ -f png
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

try:
    from .LogUtil import setup_logger
    from .plot_util import normalize_formats, save_fig
except ImportError:
    from LogUtil import setup_logger
    from plot_util import normalize_formats, save_fig

import matplotlib.pyplot as plt

logger = setup_logger("plot_combined", stream=sys.stdout)


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


def plot_per_method(pairs: List[Tuple[str, str]], col: str, method: str,
                    out_base: str, formats: Optional[Sequence[str]]) -> List[str]:
    """Combined figure: overlaid histogram + violin of `col` across samples."""
    data = []
    labels = []
    for label, tsv in pairs:
        vals = _load_column(Path(tsv), col)
        if vals:
            data.append(vals)
            labels.append(f"{label} (n={len(vals)})")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    for vals, lab in zip(data, labels):
        ax.hist(vals, bins=60, alpha=0.55, label=lab, edgecolor="none")
    ax.set_xlabel(col)
    ax.set_ylabel("count")
    ax.set_title(f"{method}: {col} across samples")
    ax.legend(fontsize=8)
    if data:
        axes[1].violinplot(data, showmedians=True)
        axes[1].set_xticks(range(1, len(labels) + 1))
        axes[1].set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
        axes[1].set_ylabel(col)
        axes[1].set_title(f"{method}: distribution")
    fig.suptitle(f"Telomere {method} — all samples", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
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
    ax.bar(methods, vals)
    ax.set_ylabel("mean value")
    ax.set_title("Cross-method telomere overview")
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
        written = plot_per_method(pairs, args.col, args.method, out_base, args.formats)
    for w in written:
        log.info(f"Wrote: {w}")


if __name__ == "__main__":
    main()
