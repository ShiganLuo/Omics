#!/usr/bin/env python3
"""Scan assembly contig ends for telomeric repeats.

Approach A: For each contig in a hifiasm assembly, scan the first/last N bp
for contiguous telomeric repeat (CCCTAA/TTAGGG) stretches. Reports per-contig
telomere length at each end.

This is the recommended approach for mouse telomeres (30-150kb) where HiFi
reads (~15-25kb) cannot span the full telomere.

Usage:
    python scan_assembly_telomere.py --fasta assembly.fa --output results/
    python scan_assembly_telomere.py --fasta assembly.fa --output results/ --scan_length 100000
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
try:
    from .LogUtil import setup_logger
    from .plot_util import normalize_formats, save_fig
except ImportError:
    from LogUtil import setup_logger
    from plot_util import normalize_formats, save_fig

import matplotlib.pyplot as plt

logger = setup_logger("scan_assembly_telomere", stream=sys.stdout)



TEL_MOTIF_FWD = "CCCTAA"
TEL_MOTIF_REV = "TTAGGG"
DEFAULT_SCAN_LENGTH = 50_000
DEFAULT_MIN_TEL_FRAC = 0.5
DEFAULT_WINDOW_SIZE = 200
DEFAULT_STEP_SIZE = 50


@dataclass
class ContigTelomere:
    """Telomere measurement for one end of a contig."""
    contig: str
    contig_length: int
    end: str
    tel_start: int
    tel_end: int
    tel_length: int
    tel_motif: str
    supporting_windows: int
    total_windows_scanned: int


def find_telomeric_region(seq: str, scan_from_start: bool,
                          scan_length: int, window_size: int,
                          step_size: int, min_tel_frac: float) -> tuple[int, int, str, int, int]:
    """Find the longest contiguous telomeric region in a sequence segment.

    Args:
        seq: DNA sequence (uppercased)
        scan_from_start: if True, scan from 5' end; if False, scan from 3' end
        scan_length: bp to scan
        window_size: sliding window size
        step_size: step size for sliding window
        min_tel_frac: minimum fraction of telomeric bases in a window

    Returns:
        (tel_start, tel_end, motif, supporting_windows, total_windows)
        Positions are relative to the scanned segment.
    """
    n = len(seq)
    scan_len = min(scan_length, n)

    if scan_from_start:
        segment = seq[:scan_len]
    else:
        segment = seq[n - scan_len:]

    seg_len = len(segment)

    tel_positions = []
    for m in re.finditer(r"(CCCTAA|TTAGGG)", segment):
        tel_positions.append((m.start(), m.end()))

    if not tel_positions:
        return (0, 0, "", 0, 0)

    clusters = []
    current_cluster = [tel_positions[0]]
    for pos in tel_positions[1:]:
        if pos[0] - current_cluster[-1][1] < 30:
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)

    best_cluster = None
    best_tel_bases = 0
    for cluster in clusters:
        tel_bases = sum(e - s for s, e in cluster)
        if tel_bases > best_tel_bases:
            best_tel_bases = tel_bases
            best_cluster = cluster

    if best_cluster is None:
        return (0, 0, "", 0, 0)

    cluster_start = best_cluster[0][0]
    cluster_end = best_cluster[-1][1]
    cluster_len = cluster_end - cluster_start

    ccctaa_count = 0
    ttaggg_count = 0
    for m in re.finditer(r"CCCTAA", segment[cluster_start:cluster_end]):
        ccctaa_count += 1
    for m in re.finditer(r"TTAGGG", segment[cluster_start:cluster_end]):
        ttaggg_count += 1
    motif = TEL_MOTIF_FWD if ccctaa_count >= ttaggg_count else TEL_MOTIF_REV

    supporting = 0
    total_windows = 0
    for start in range(0, seg_len - window_size + 1, step_size):
        window = segment[start:start + window_size]
        tel_count = window.count("CCCTAA") + window.count("TTAGGG")
        tel_frac = (tel_count * 6) / window_size
        total_windows += 1
        if tel_frac >= min_tel_frac:
            supporting += 1

    if not scan_from_start:
        tel_start = n - scan_len + cluster_start
        tel_end = n - scan_len + cluster_end
    else:
        tel_start = cluster_start
        tel_end = cluster_end

    return (tel_start, tel_end, motif, supporting, total_windows)


def scan_assembly(fasta_path: str, scan_length: int,
                  window_size: int, step_size: int,
                  min_tel_frac: float) -> list[ContigTelomere]:
    """Scan all contigs in an assembly for telomeric repeats at ends."""
    results = []

    contigs = {}
    current_name = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    contigs[current_name] = "".join(current_seq).upper()
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name is not None:
        contigs[current_name] = "".join(current_seq).upper()

    logger.info(f"Scanning {len(contigs)} contigs for telomeric repeats "
          f"(scan_length={scan_length}bp, window={window_size}bp, step={step_size}bp)")

    for name, seq in contigs.items():
        contig_len = len(seq)

        tel5_start, tel5_end, motif5, sup5, tot5 = find_telomeric_region(
            seq, scan_from_start=True, scan_length=scan_length,
            window_size=window_size, step_size=step_size, min_tel_frac=min_tel_frac
        )
        tel5_len = tel5_end - tel5_start

        tel3_start, tel3_end, motif3, sup3, tot3 = find_telomeric_region(
            seq, scan_from_start=False, scan_length=scan_length,
            window_size=window_size, step_size=step_size, min_tel_frac=min_tel_frac
        )
        tel3_len = tel3_end - tel3_start

        if tel5_len > 100 or sup5 >= 3:
            results.append(ContigTelomere(
                contig=name, contig_length=contig_len, end="5prime",
                tel_start=tel5_start, tel_end=tel5_end, tel_length=tel5_len,
                tel_motif=motif5, supporting_windows=sup5, total_windows_scanned=tot5
            ))

        if tel3_len > 100 or sup3 >= 3:
            results.append(ContigTelomere(
                contig=name, contig_length=contig_len, end="3prime",
                tel_start=tel3_start, tel_end=tel3_end, tel_length=tel3_len,
                tel_motif=motif3, supporting_windows=sup3, total_windows_scanned=tot3
            ))

    logger.info(f"Found {len(results)} contig ends with telomeric signal")
    return results


def plot_telomere_length(results: list[ContigTelomere], sample_name: str,
                         out_base: str,
                         formats: Optional[List[str]] = None) -> List[str]:
    """Figure 2: assembly contig-end telomere length distribution (multi-format).

    Panels: (1) histogram of telomere length across contig ends; (2) per-contig
    end length scatter (5'/3' distinguished). Assembly/reference telomere length
    (not sample length — see README).
    """
    if not results:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "no telomeric regions detected at contig ends",
                ha="center", va="center")
        ax.set_axis_off()
        fig.suptitle(f"Assembly telomere length — {sample_name}", fontsize=12)
        written = save_fig(fig, out_base, formats)
        plt.close(fig)
        return written
    sig = [r for r in results if r.tel_length >= 1000]
    show = sig if sig else results
    lengths = [r.tel_length for r in show]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    ax = axes[0]
    if lengths:
        ax.hist(lengths, bins=40, color="#4C9F70", edgecolor="none")
    ax.set_xlabel("Telomere length (bp)")
    ax.set_ylabel("Contig ends")
    ax.set_title(f"Telomere length distribution\n(>=1kb, n={len(sig)})")

    ax = axes[1]
    colors = {"5prime": "#3C6E9F", "3prime": "#C0504D"}
    marks = {"5prime": "o", "3prime": "^"}
    for end in ("5prime", "3prime"):
        xs = [i for i, r in enumerate(show) if r.end == end]
        ys = [r.tel_length for r in show if r.end == end]
        if xs:
            ax.scatter(xs, ys, s=12, c=colors.get(end, "gray"),
                       marker=marks.get(end, "o"), label=end, alpha=0.7)
    ax.set_xlabel("Contig end index")
    ax.set_ylabel("Telomere length (bp)")
    ax.set_title("Per-contig-end telomere length\n(5' vs 3')")
    ax.legend()

    fig.suptitle(f"Assembly telomere length — {sample_name}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    written = save_fig(fig, out_base, formats)
    plt.close(fig)
    return written


def write_results(results: list[ContigTelomere], sample_name: str,
                  output_dir: Path, scan_length: int,
                  formats: Optional[List[str]] = None):
    """Write results to TSV and summary files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    detail_path = output_dir / f"{sample_name}_assembly_telomere.tsv"
    with open(detail_path, "w") as f:
        f.write("contig\tcontig_length\tend\ttel_start\ttel_end\t"
                "tel_length\ttel_motif\tsupporting_windows\ttotal_windows\n")
        for r in results:
            f.write(f"{r.contig}\t{r.contig_length}\t{r.end}\t"
                    f"{r.tel_start}\t{r.tel_end}\t{r.tel_length}\t"
                    f"{r.tel_motif}\t{r.supporting_windows}\t{r.total_windows_scanned}\n")
    logger.info(f"Wrote: {detail_path}")

    stats_path = output_dir / f"{sample_name}_assembly_telomere_stats.txt"
    with open(stats_path, "w") as f:
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Method: Assembly contig end scanning (scan_length={scan_length}bp)\n\n")

        f.write(f"Total contig ends with telomeric signal: {len(results)}\n\n")

        sig_results = [r for r in results if r.tel_length >= 1000]
        f.write(f"Contig ends with telomere >= 1kb: {len(sig_results)}\n\n")

        if sig_results:
            import numpy as np
            lengths = np.array([r.tel_length for r in sig_results])

            f.write("=== Telomere Length Statistics (>= 1kb) ===\n")
            f.write(f"  Count:   {len(lengths)}\n")
            f.write(f"  Mean:    {np.mean(lengths):.0f} bp\n")
            f.write(f"  Median:  {np.median(lengths):.0f} bp\n")
            f.write(f"  Std:     {np.std(lengths):.0f} bp\n")
            f.write(f"  Min:     {np.min(lengths):.0f} bp\n")
            f.write(f"  Max:     {np.max(lengths):.0f} bp\n")
            f.write(f"  Q25:     {np.percentile(lengths, 25):.0f} bp\n")
            f.write(f"  Q75:     {np.percentile(lengths, 75):.0f} bp\n\n")

            prime5 = [r for r in sig_results if r.end == "5prime"]
            prime3 = [r for r in sig_results if r.end == "3prime"]

            if prime5:
                arr5 = np.array([r.tel_length for r in prime5])
                f.write("=== 5' End Telomere Length ===\n")
                f.write(f"  Count:   {len(arr5)}\n")
                f.write(f"  Mean:    {np.mean(arr5):.0f} bp\n")
                f.write(f"  Median:  {np.median(arr5):.0f} bp\n\n")

            if prime3:
                arr3 = np.array([r.tel_length for r in prime3])
                f.write("=== 3' End Telomere Length ===\n")
                f.write(f"  Count:   {len(arr3)}\n")
                f.write(f"  Mean:    {np.mean(arr3):.0f} bp\n")
                f.write(f"  Median:  {np.median(arr3):.0f} bp\n\n")

            f.write("=== Length Distribution ===\n")
            thresholds = [0, 1000, 5000, 10000, 20000, 50000, 100000, 200000]
            for i in range(len(thresholds) - 1):
                lo, hi = thresholds[i], thresholds[i + 1]
                count = sum(1 for r in sig_results if lo <= r.tel_length < hi)
                f.write(f"  {lo:>8,} - {hi:>8,} bp:  {count}\n")
            count = sum(1 for r in sig_results if r.tel_length >= thresholds[-1])
            f.write(f"  >{thresholds[-1]:>7,} bp:  {count}\n")

            f.write("\n=== Top 20 Longest Telomeres ===\n")
            f.write(f"{'Contig':<20} {'End':<8} {'Length':>10} {'Motif':<8} {'Windows':>8}\n")
            for r in sorted(sig_results, key=lambda x: x.tel_length, reverse=True)[:20]:
                f.write(f"{r.contig:<20} {r.end:<8} {r.tel_length:>10,} "
                        f"{r.tel_motif:<8} {r.supporting_windows:>5}/{r.total_windows_scanned}\n")
        else:
            f.write("No significant telomeric regions found at contig ends.\n")
            f.write("This may indicate:\n")
            f.write("  - Assembly did not span telomeric regions\n")
            f.write("  - Telomeres are too short for detection\n")
            f.write("  - Contigs are internal (no chromosome ends)\n")

    logger.info(f"Wrote: {stats_path}")

    with open(stats_path) as f:
        logger.info(f.read())

    plot_base = str(output_dir / f"{sample_name}_assembly_telomere_length")
    for _p in plot_telomere_length(results, sample_name, plot_base, formats):
        logger.info(f"Wrote: {_p}")


def main():
    global logger
    parser = argparse.ArgumentParser(
        description="Scan assembly contig ends for telomeric repeats")
    parser.add_argument("-f", "--fasta", required=True,
                        help="Input assembly FASTA file")
    parser.add_argument("-s", "--sample_name", required=True,
                        help="Sample name for output files")
    parser.add_argument("-o", "--output_dir", default="results",
                        help="Output directory")
    parser.add_argument("-L", "--scan_length", type=int, default=DEFAULT_SCAN_LENGTH,
                        help=f"bp to scan from each contig end (default: {DEFAULT_SCAN_LENGTH})")
    parser.add_argument("-w", "--window_size", type=int, default=DEFAULT_WINDOW_SIZE,
                        help=f"Sliding window size (default: {DEFAULT_WINDOW_SIZE})")
    parser.add_argument("-e", "--step_size", type=int, default=DEFAULT_STEP_SIZE,
                        help=f"Step size for sliding window (default: {DEFAULT_STEP_SIZE})")
    parser.add_argument("-t", "--min_tel_frac", type=float, default=DEFAULT_MIN_TEL_FRAC,
                        help=f"Minimum telomeric fraction in window (default: {DEFAULT_MIN_TEL_FRAC})")
    parser.add_argument("-l", "--log", default="",
                        help="Optional log file path (default: log to stdout only)")
    parser.add_argument("-F", "--format", action="append", dest="formats", metavar="FMT",
                        help="Output figure format(s), repeatable (default: png)")
    args = parser.parse_args()

    if args.log:
        Path(args.log).parent.mkdir(parents=True, exist_ok=True)
        logger = setup_logger("scan_assembly_telomere", log_file=args.log, stream=sys.stdout)

    results = scan_assembly(
        args.fasta,
        scan_length=args.scan_length,
        window_size=args.window_size,
        step_size=args.step_size,
        min_tel_frac=args.min_tel_frac
    )
    write_results(results, args.sample_name, Path(args.output_dir),
                  args.scan_length, args.formats)


if __name__ == "__main__":
    main()