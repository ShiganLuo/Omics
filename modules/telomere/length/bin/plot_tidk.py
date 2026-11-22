#!/usr/bin/env python3
"""Plot telomere repeat locations/lengths from a tidk ``search`` result (Figure 3).

Reads tidk's ``*_telomeric_repeat.tsv`` and renders a per-contig telomere
picture. The exact tidk column layout varies across versions, so this parser is
ADAPTIVE: it detects the contig/sequence column and numeric length/size/count
columns by name + dtype, then plots the most length-like numeric column (with a
per-contig panel). If the layout can't be parsed it logs a warning and writes a
"no parsable data" placeholder instead of crashing.

Usage:
    python plot_tidk.py -i s1_tidk_search_telomeric_repeat.tsv \
        -s s1 -o out/s1/tidk -F png -F svg
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Tuple

try:
    from .LogUtil import setup_logger
    from .plot_util import normalize_formats, save_fig
except ImportError:
    from LogUtil import setup_logger
    from plot_util import normalize_formats, save_fig

import matplotlib.pyplot as plt

logger = setup_logger("plot_tidk", stream=sys.stdout)

_CONTIG_HINTS = ("contig", "sequence", "seq", "chrom", "chr", "name", "id", "scaffold")
_LEN_HINTS = ("length", "len", "size", "span", "bp")
_COUNT_HINTS = ("count", "repeat", "copies", "n_")


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def parse_tidk_tsv(path: Path) -> Tuple[Optional[str], Optional[str], dict]:
    """Adaptively parse a tidk TSV.

    Returns (contig_col, value_col, rows) where rows maps contig -> list of
    numeric values on value_col. Picks the length-like numeric column when
    present, else the first numeric column.
    """
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows_raw = [r for r in reader if r and any(c.strip() for c in r)]
    if not rows_raw:
        return None, None, {}
    header = [c.strip() for c in rows_raw[0]]
    data = rows_raw[1:]
    ncol = len(header)

    # column classification
    numeric_cols = []
    for j in range(ncol):
        vals = [r[j] for r in data if j < len(r)]
        if vals and sum(_is_number(v) for v in vals) / len(vals) > 0.8:
            numeric_cols.append(j)
    contig_col = None
    for j in range(ncol):
        if any(h in header[j].lower() for h in _CONTIG_HINTS):
            contig_col = j
            break
    if contig_col is None and not numeric_cols:
        contig_col = 0
    if contig_col is None:
        contig_col = 0

    # value column: prefer length-like, then count-like, else first numeric
    value_col = None
    for j in numeric_cols:
        if any(h in header[j].lower() for h in _LEN_HINTS):
            value_col = j
            break
    if value_col is None:
        for j in numeric_cols:
            if any(h in header[j].lower() for h in _COUNT_HINTS):
                value_col = j
                break
    if value_col is None and numeric_cols:
        value_col = numeric_cols[0]
    if value_col is None:
        return header[contig_col] if contig_col < ncol else None, None, {}

    rows: dict = defaultdict(list)
    for r in data:
        if contig_col < len(r) and value_col < len(r) and _is_number(r[value_col]):
            rows[r[contig_col]].append(float(r[value_col]))
    return header[contig_col], header[value_col], dict(rows)


def plot_tidk(contig_col: Optional[str], value_col: Optional[str], rows: dict,
              sample_name: str, out_base: str,
              formats: Optional[List[str]] = None) -> List[str]:
    """Figure 3: telomere repeat length distribution + per-contig (multi-format)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    all_vals = [v for vs in rows.values() for v in vs]

    ax = axes[0]
    if all_vals:
        ax.hist(all_vals, bins=40, color="#3C6E9F", edgecolor="none")
    else:
        ax.text(0.5, 0.5, "no parsable data", ha="center", va="center")
    ax.set_xlabel(value_col or "value")
    ax.set_ylabel("Count")
    ax.set_title(f"tidk telomere repeat — {value_col or 'value'}")

    ax = axes[1]
    contigs = list(rows.keys())
    means = [sum(rows[c]) / len(rows[c]) if rows[c] else 0 for c in contigs]
    if contigs:
        ax.bar(range(len(contigs)), means, color="#4C9F70")
        ax.set_xlabel("Contig index")
        ax.set_ylabel(f"mean {value_col or 'value'}")
        ax.set_title(f"Per-contig (n={len(contigs)})")
    else:
        ax.text(0.5, 0.5, "no contig data", ha="center", va="center")

    fig.suptitle(f"tidk telomere repeats — {sample_name}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def main() -> None:
    """Parse CLI, read tidk TSV, render Figure 3."""
    parser = argparse.ArgumentParser(
        description="Plot tidk telomere repeat result (Figure 3).")
    parser.add_argument("-i", "--input", required=True,
                        help="tidk *_telomeric_repeat.tsv")
    parser.add_argument("-s", "--sample_name", required=True,
                        help="Sample name (for titles)")
    parser.add_argument("-o", "--output_dir", default="results",
                        help="Output directory")
    parser.add_argument("-b", "--out_base", default="",
                        help="Output figure base name (default: <sample>_tidk_telomeres)")
    parser.add_argument("-F", "--format", action="append", dest="formats", metavar="FMT",
                        help="Output figure format(s), repeatable (default: png)")
    parser.add_argument("-l", "--log", default="",
                        help="Optional log file path (default: log to stdout only)")
    args = parser.parse_args()

    log = logger
    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        log = setup_logger("plot_tidk", log_file=args.log, stream=sys.stdout)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    contig_col, value_col, rows = parse_tidk_tsv(Path(args.input))
    if not rows:
        log.warning(f"No parsable rows in {args.input}; writing placeholder figure.")
    base = args.out_base or f"{args.sample_name}_tidk_telomeres"
    out_base = str(out / base)
    written = plot_tidk(contig_col, value_col, rows, args.sample_name,
                        out_base, args.formats)
    for p in written:
        log.info(f"Wrote: {p}")


if __name__ == "__main__":
    main()
