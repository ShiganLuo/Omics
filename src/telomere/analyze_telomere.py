#!/usr/bin/env python3
"""Analyze telomere length from PacBio HiFi BAM files.

Extracts reads containing telomeric repeats (TTAGGG)n and measures
the length of telomeric sequence at read termini.

Usage:
    python analyze_telomere.py --bam input.bam --sample_name DMSO_P20 --output_dir results/
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TEL_PATTERN = re.compile(r"(TTAGGG|CCCTAA){5,}", re.IGNORECASE)
# For mouse, telomeres are long; use a sliding window to find contiguous
# telomeric stretches allowing small interruptions
WINDOW_SIZE = 100
STEP_SIZE = 20
MIN_TEL_CONTENT = 0.6  # 60% of window must be telomeric


@dataclass
class TelomereHit:
    """A telomere repeat region found in a read."""
    read_name: str
    read_length: int
    strand: str  # "G-rich" (TTAGGG) or "C-rich" (CCCTAA)
    start: int
    end: int
    tel_length: int
    position: str  # "5prime", "3prime", "internal"


@dataclass
class ReadTelomere:
    """Summary of telomere content in a single read."""
    read_name: str
    read_length: int
    has_telomere: bool
    prime5_tel_len: int  # telomere length at 5' end
    prime3_tel_len: int  # telomere length at 3' end
    internal_tel_len: int  # internal telomeric sequence
    hits: list[TelomereHit] = field(default_factory=list)


def find_telomeric_regions(seq: str) -> list[tuple[int, int, str]]:
    """Find all telomeric repeat regions in a sequence.

    Uses regex for exact repeats, then merges nearby hits.
    Returns list of (start, end, strand_type).
    """
    hits = []
    for m in TEL_PATTERN.finditer(seq):
        strand = "G-rich" if m.group().upper().startswith("TTAGGG") else "C-rich"
        hits.append((m.start(), m.end(), strand))

    # Merge overlapping/nearby hits (within 20bp)
    if not hits:
        return []

    merged = [hits[0]]
    for start, end, strand in hits[1:]:
        prev_start, prev_end, prev_strand = merged[-1]
        if start <= prev_end + 20 and strand == prev_strand:
            merged[-1] = (prev_start, max(end, prev_end), prev_strand)
        else:
            merged.append((start, end, strand))
    return merged


def find_telomere_with_gaps(seq: str, min_pct: float = 0.55) -> list[tuple[int, int, str]]:
    """Find telomeric regions allowing gaps (for real HiFi data with errors).

    Slides a 100bp window; if >55% of kmers are telomeric, marks as telomeric.
    This catches telomeric regions with occasional non-TTAGGG motifs.
    """
    results = []
    seq_upper = seq.upper()
    n = len(seq_upper)

    # Quick check: if no TTAGGG/CCCTAA at all, skip
    if "TTAGGG" not in seq_upper and "CCCTAA" not in seq_upper:
        return []

    # Find regions by scanning for TTAGGG/CCCTAA density
    tel_positions = []
    for m in re.finditer(r"(TTAGGG|CCCTAA)", seq_upper):
        tel_positions.append((m.start(), m.end()))

    if not tel_positions:
        return []

    # Cluster nearby hits
    clusters = []
    current_cluster = [tel_positions[0]]
    for pos in tel_positions[1:]:
        if pos[0] - current_cluster[-1][1] < 30:  # allow 30bp gap
            current_cluster.append(pos)
        else:
            clusters.append(current_cluster)
            current_cluster = [pos]
    clusters.append(current_cluster)

    for cluster in clusters:
        start = cluster[0][0]
        end = cluster[-1][1]
        region_len = end - start
        tel_bases = sum(e - s for s, e in cluster)
        tel_frac = tel_bases / region_len if region_len > 0 else 0

        if region_len >= 30 and tel_frac >= min_pct:
            # Determine strand from first hit
            first_seq = seq_upper[cluster[0][0]:cluster[0][1]]
            strand = "G-rich" if first_seq == "TTAGGG" else "C-rich"
            results.append((start, end, strand))

    # Merge overlapping results
    if not results:
        return []
    merged = [results[0]]
    for start, end, strand in results[1:]:
        ps, pe, pstrand = merged[-1]
        if start <= pe + 50:
            merged[-1] = (ps, max(end, pe), pstrand)
        else:
            merged.append((start, end, strand))
    return merged


def analyze_read(read_name: str, seq: str) -> ReadTelomere:
    """Analyze a single read for telomeric content."""
    read_len = len(seq)
    regions = find_telomere_with_gaps(seq)

    prime5_tel = 0
    prime3_tel = 0
    internal_tel = 0
    hits = []

    for start, end, strand in regions:
        tel_len = end - start
        # Classify position
        if start < 500:
            position = "5prime"
            prime5_tel = max(prime5_tel, tel_len)
        elif end > read_len - 500:
            position = "3prime"
            prime3_tel = max(prime3_tel, tel_len)
        else:
            position = "internal"
            internal_tel += tel_len

        hits.append(TelomereHit(
            read_name=read_name,
            read_length=read_len,
            strand=strand,
            start=start,
            end=end,
            tel_length=tel_len,
            position=position,
        ))

    return ReadTelomere(
        read_name=read_name,
        read_length=read_len,
        has_telomere=len(hits) > 0,
        prime5_tel_len=prime5_tel,
        prime3_tel_len=prime3_tel,
        internal_tel_len=internal_tel,
        hits=hits,
    )


def process_bam(bam_path: str, min_read_length: int = 5000) -> list[ReadTelomere]:
    """Process all reads in a BAM file."""
    results = []
    n_total = 0
    n_tel = 0

    print(f"Reading BAM: {bam_path}", file=sys.stderr)
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped and read.query_sequence is None:
                continue
            n_total += 1
            if n_total % 100000 == 0:
                print(f"  Processed {n_total:,} reads, {n_tel:,} with telomere...", file=sys.stderr)

            seq = read.query_sequence
            if seq is None or len(seq) < min_read_length:
                continue

            result = analyze_read(read.query_name, seq)
            if result.has_telomere:
                n_tel += 1
                results.append(result)

    print(f"Done. Total reads: {n_total:,}, with telomere: {n_tel:,}", file=sys.stderr)
    return results


def write_results(results: list[ReadTelomere], sample_name: str, output_dir: Path):
    """Write analysis results to TSV files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Per-read summary
    summary_path = output_dir / f"{sample_name}_telomere_reads.tsv"
    with open(summary_path, "w") as f:
        f.write("read_name\tread_length\tprime5_tel_len\tprime3_tel_len\tinternal_tel_len\tn_hits\n")
        for r in results:
            f.write(f"{r.read_name}\t{r.read_length}\t{r.prime5_tel_len}\t"
                    f"{r.prime3_tel_len}\t{r.internal_tel_len}\t{len(r.hits)}\n")
    print(f"Wrote: {summary_path}", file=sys.stderr)

    # 2. Per-hit detail
    detail_path = output_dir / f"{sample_name}_telomere_hits.tsv"
    with open(detail_path, "w") as f:
        f.write("read_name\tread_length\tstrand\tstart\tend\ttel_length\tposition\n")
        for r in results:
            for h in r.hits:
                f.write(f"{h.read_name}\t{h.read_length}\t{h.strand}\t"
                        f"{h.start}\t{h.end}\t{h.tel_length}\t{h.position}\n")
    print(f"Wrote: {detail_path}", file=sys.stderr)

    # 3. Statistical summary
    stats_path = output_dir / f"{sample_name}_telomere_stats.txt"
    with open(stats_path, "w") as f:
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Total reads with telomere: {len(results)}\n\n")

        # Terminal telomere lengths (most biologically relevant)
        # Filter: >1000bp to exclude subtelomeric noise / ITS fragments
        # Mouse telomeres (129/Ola) are typically 30-150kb
        terminal_5 = [r.prime5_tel_len for r in results if r.prime5_tel_len > 1000]
        terminal_3 = [r.prime3_tel_len for r in results if r.prime3_tel_len > 1000]
        all_terminal = terminal_5 + terminal_3

        if all_terminal:
            arr = np.array(all_terminal)
            f.write("=== Terminal Telomere Length (reads with >1000bp telomere at ends) ===\n")
            f.write(f"  Count:   {len(arr)}\n")
            f.write(f"  Mean:    {np.mean(arr):.0f} bp\n")
            f.write(f"  Median:  {np.median(arr):.0f} bp\n")
            f.write(f"  Std:     {np.std(arr):.0f} bp\n")
            f.write(f"  Min:     {np.min(arr):.0f} bp\n")
            f.write(f"  Max:     {np.max(arr):.0f} bp\n")
            f.write(f"  Q25:     {np.percentile(arr, 25):.0f} bp\n")
            f.write(f"  Q75:     {np.percentile(arr, 75):.0f} bp\n\n")

        if terminal_5:
            arr = np.array(terminal_5)
            f.write("=== 5' Telomere Length ===\n")
            f.write(f"  Count:   {len(arr)}\n")
            f.write(f"  Mean:    {np.mean(arr):.0f} bp\n")
            f.write(f"  Median:  {np.median(arr):.0f} bp\n")
            f.write(f"  Max:     {np.max(arr):.0f} bp\n\n")

        if terminal_3:
            arr = np.array(terminal_3)
            f.write("=== 3' Telomere Length ===\n")
            f.write(f"  Count:   {len(arr)}\n")
            f.write(f"  Mean:    {np.mean(arr):.0f} bp\n")
            f.write(f"  Median:  {np.median(arr):.0f} bp\n")
            f.write(f"  Max:     {np.max(arr):.0f} bp\n\n")

        # Length distribution
        f.write("=== Terminal Telomere Length Distribution ===\n")
        thresholds = [0, 1000, 5000, 10000, 20000, 50000, 100000, 200000]
        for i in range(len(thresholds) - 1):
            lo, hi = thresholds[i], thresholds[i + 1]
            count = sum(1 for x in all_terminal if lo <= x < hi)
            f.write(f"  {lo:>8,} - {hi:>8,} bp:  {count}\n")
        count = sum(1 for x in all_terminal if x >= thresholds[-1])
        f.write(f"  >{thresholds[-1]:>7,} bp:  {count}\n")

    print(f"Wrote: {stats_path}", file=sys.stderr)

    # 4. Quick print to stdout
    with open(stats_path) as f:
        print(f.read())


def main():
    parser = argparse.ArgumentParser(description="Analyze telomere length from PacBio HiFi BAM")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--sample_name", required=True, help="Sample name for output files")
    parser.add_argument("--output_dir", default="results", help="Output directory")
    parser.add_argument("--min_read_length", type=int, default=5000,
                        help="Minimum read length to consider (default: 5000)")
    args = parser.parse_args()

    results = process_bam(args.bam, min_read_length=args.min_read_length)
    write_results(results, args.sample_name, Path(args.output_dir))


if __name__ == "__main__":
    main()
