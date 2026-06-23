#!/usr/bin/env python3
"""Estimate genome-wide average telomere length from read-level k-mer density.

Approach B: Count telomeric k-mers (CCCTAA/TTAGGG) across all HiFi reads.
Divide total telomeric bp by estimated number of chromosome arms to get
average telomere length per arm.

This gives a genome-wide average, not per-chromosome-arm measurements.

Usage:
    python read_density_telomere.py --bam input.bam --sample_name SAMPLE --output results/
    python read_density_telomere.py --bam input.bam --sample_name SAMPLE --output results/ \
        --n_chrom_arms 40 --min_read_length 5000
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TEL_KMER_FWD = "CCCTAA"
TEL_KMER_REV = "TTAGGG"
# Mouse: 20 chromosomes, diploid = 40 chromosome arms
DEFAULT_N_CHROM_ARMS = 40
DEFAULT_MIN_READ_LENGTH = 5000
# Minimum telomeric k-mer count to consider a read as having telomeric content
DEFAULT_MIN_TEL_KMERS = 16


@dataclass
class ReadTelStats:
    """Telomere statistics for a single read."""
    read_name: str
    read_length: int
    tel_kmer_count: int       # total telomeric k-mers found
    tel_bp_estimate: int      # estimated telomeric bp (kmer_count * 6)
    tel_fraction: float       # fraction of read that is telomeric
    has_terminal_tel: bool    # whether telomeric signal is at read end
    terminal_tel_len: int     # length of terminal telomeric region (bp)


def count_tel_kmers(seq: str) -> tuple[int, list[tuple[int, int]]]:
    """Count telomeric k-mers and find their positions.

    Returns:
        (total_kmer_count, list of (start, end) positions)
    """
    positions = []
    seq_upper = seq.upper()
    for m in re.finditer(r"(CCCTAA|TTAGGG)", seq_upper):
        positions.append((m.start(), m.end()))
    return len(positions), positions


def find_terminal_telomere(seq: str, positions: list[tuple[int, int]],
                           max_gap: int = 30) -> int:
    """Find the length of terminal telomeric region.

    Clusters telomeric k-mers near the read end (within 500bp) and
    returns the total length of the terminal telomeric stretch.
    """
    if not positions:
        return 0

    read_len = len(seq)
    # Find clusters near the 3' end (last 500bp)
    terminal_positions = [p for p in positions if p[0] > read_len - 500]
    if not terminal_positions:
        # Check 5' end
        terminal_positions = [p for p in positions if p[1] < 500]
    if not terminal_positions:
        return 0

    # Cluster nearby hits
    clusters = []
    current_cluster = [terminal_positions[0]]
    for pos in terminal_positions[1:]:
        if pos[0] - current_cluster[-1][1] < max_gap:
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)

    # Find the largest cluster
    best_len = 0
    for cluster in clusters:
        start = cluster[0][0]
        end = cluster[-1][1]
        tel_len = end - start
        best_len = max(best_len, tel_len)

    return best_len


def analyze_reads(bam_path: str, min_read_length: int,
                  min_tel_kmers: int) -> list[ReadTelStats]:
    """Analyze all reads in a BAM file for telomeric content."""
    results = []
    n_total = 0
    n_tel = 0

    print(f"Reading BAM: {bam_path}", file=sys.stderr)
    print(f"  min_read_length={min_read_length}, min_tel_kmers={min_tel_kmers}",
          file=sys.stderr)

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped and read.query_sequence is None:
                continue
            n_total += 1
            if n_total % 500000 == 0:
                print(f"  Processed {n_total:,} reads, {n_tel:,} with telomere...",
                      file=sys.stderr)

            seq = read.query_sequence
            if seq is None or len(seq) < min_read_length:
                continue

            kmer_count, positions = count_tel_kmers(seq)
            if kmer_count < min_tel_kmers:
                continue

            n_tel += 1
            read_len = len(seq)
            tel_bp = kmer_count * 6  # each k-mer is 6bp
            tel_frac = tel_bp / read_len if read_len > 0 else 0
            terminal_len = find_terminal_telomere(seq, positions)
            has_terminal = terminal_len > 0

            results.append(ReadTelStats(
                read_name=read.query_name,
                read_length=read_len,
                tel_kmer_count=kmer_count,
                tel_bp_estimate=tel_bp,
                tel_fraction=tel_frac,
                has_terminal_tel=has_terminal,
                terminal_tel_len=terminal_len
            ))

    print(f"Done. Total reads: {n_total:,}, with telomere: {n_tel:,}",
          file=sys.stderr)
    return results


def write_results(results: list[ReadTelStats], sample_name: str,
                  output_dir: Path, n_chrom_arms: int):
    """Write results to TSV and summary files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Per-read detail TSV
    detail_path = output_dir / f"{sample_name}_read_telomere.tsv"
    with open(detail_path, "w") as f:
        f.write("read_name\tread_length\ttel_kmer_count\ttel_bp_estimate\t"
                "tel_fraction\thas_terminal_tel\tterminal_tel_len\n")
        for r in results:
            f.write(f"{r.read_name}\t{r.read_length}\t{r.tel_kmer_count}\t"
                    f"{r.tel_bp_estimate}\t{r.tel_fraction:.4f}\t"
                    f"{r.has_terminal_tel}\t{r.terminal_tel_len}\n")
    print(f"Wrote: {detail_path}", file=sys.stderr)

    # 2. Statistical summary
    stats_path = output_dir / f"{sample_name}_read_telomere_stats.txt"
    with open(stats_path, "w") as f:
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Method: Read-level k-mer density estimation\n")
        f.write(f"N_chromosome_arms: {n_chrom_arms}\n\n")

        f.write(f"Reads with telomeric signal: {len(results)}\n\n")

        if not results:
            f.write("No reads with telomeric signal found.\n")
            print(f"Wrote: {stats_path}", file=sys.stderr)
            with open(stats_path) as ff:
                print(ff.read())
            return

        total_tel_bp = sum(r.tel_bp_estimate for r in results)
        total_read_bp = sum(r.read_length for r in results)

        f.write("=== Genome-wide Telomere Content ===\n")
        f.write(f"  Total telomeric bp:     {total_tel_bp:>15,} bp\n")
        f.write(f"  Total read bp:          {total_read_bp:>15,} bp\n")
        f.write(f"  Telomeric fraction:     {total_tel_bp/total_read_bp:.6f}\n")
        f.write(f"  Avg telomere per arm:   {total_tel_bp / n_chrom_arms:>15,.0f} bp\n\n")

        # Per-read statistics
        tel_bps = np.array([r.tel_bp_estimate for r in results])
        tel_fracs = np.array([r.tel_fraction for r in results])

        f.write("=== Per-Read Telomeric Content ===\n")
        f.write(f"  Telomeric bp per read:\n")
        f.write(f"    Mean:    {np.mean(tel_bps):.0f} bp\n")
        f.write(f"    Median:  {np.median(tel_bps):.0f} bp\n")
        f.write(f"    Max:     {np.max(tel_bps):.0f} bp\n")
        f.write(f"  Telomeric fraction per read:\n")
        f.write(f"    Mean:    {np.mean(tel_fracs):.4f}\n")
        f.write(f"    Median:  {np.median(tel_fracs):.4f}\n\n")

        # Terminal telomere statistics
        terminal_reads = [r for r in results if r.has_terminal_tel and r.terminal_tel_len > 0]
        f.write(f"Reads with terminal telomere: {len(terminal_reads)}\n\n")

        if terminal_reads:
            terminal_lens = np.array([r.terminal_tel_len for r in terminal_reads])
            f.write("=== Terminal Telomere Length (reads with telomere at end) ===\n")
            f.write(f"  Count:   {len(terminal_lens)}\n")
            f.write(f"  Mean:    {np.mean(terminal_lens):.0f} bp\n")
            f.write(f"  Median:  {np.median(terminal_lens):.0f} bp\n")
            f.write(f"  Std:     {np.std(terminal_lens):.0f} bp\n")
            f.write(f"  Min:     {np.min(terminal_lens):.0f} bp\n")
            f.write(f"  Max:     {np.max(terminal_lens):.0f} bp\n")
            f.write(f"  Q25:     {np.percentile(terminal_lens, 25):.0f} bp\n")
            f.write(f"  Q75:     {np.percentile(terminal_lens, 75):.0f} bp\n\n")

            # Length distribution
            f.write("=== Terminal Telomere Length Distribution ===\n")
            thresholds = [0, 500, 1000, 5000, 10000, 20000, 50000, 100000]
            for i in range(len(thresholds) - 1):
                lo, hi = thresholds[i], thresholds[i + 1]
                count = sum(1 for x in terminal_lens if lo <= x < hi)
                f.write(f"  {lo:>8,} - {hi:>8,} bp:  {count}\n")
            count = sum(1 for x in terminal_lens if x >= thresholds[-1])
            f.write(f"  >{thresholds[-1]:>7,} bp:  {count}\n")

    print(f"Wrote: {stats_path}", file=sys.stderr)

    # Print summary to stdout
    with open(stats_path) as f:
        print(f.read())


def main():
    parser = argparse.ArgumentParser(
        description="Estimate genome-wide telomere length from read-level k-mer density")
    parser.add_argument("--bam", required=True,
                        help="Input BAM file")
    parser.add_argument("--sample_name", required=True,
                        help="Sample name for output files")
    parser.add_argument("--output_dir", default="results",
                        help="Output directory")
    parser.add_argument("--n_chrom_arms", type=int, default=DEFAULT_N_CHROM_ARMS,
                        help=f"Number of chromosome arms (default: {DEFAULT_N_CHROM_ARMS})")
    parser.add_argument("--min_read_length", type=int, default=DEFAULT_MIN_READ_LENGTH,
                        help=f"Minimum read length (default: {DEFAULT_MIN_READ_LENGTH})")
    parser.add_argument("--min_tel_kmers", type=int, default=DEFAULT_MIN_TEL_KMERS,
                        help=f"Minimum telomeric k-mers per read (default: {DEFAULT_MIN_TEL_KMERS})")
    args = parser.parse_args()

    results = analyze_reads(
        args.bam,
        min_read_length=args.min_read_length,
        min_tel_kmers=args.min_tel_kmers
    )
    write_results(results, args.sample_name, Path(args.output_dir), args.n_chrom_arms)


if __name__ == "__main__":
    main()
