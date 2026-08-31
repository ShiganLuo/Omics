#!/usr/bin/env python3
"""Find overlaps between narrowPeak and TE annotation GTF.

Uses bedtools intersect to find peak-TE overlaps, then counts
overlaps per TE class (SINE, LINE, LTR, DNA, etc.).

Outputs:
  - overlap BED: full intersect result
  - overlap TSV: sample_id, te_class, overlap_count, total_peaks, overlap_ratio
"""
import argparse
import logging
import os
import subprocess
import sys
from collections import Counter

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("intersect_te")


def parse_gtf_attribute(attr_str, key):
    """Extract a value from GTF attributes string."""
    for part in attr_str.rstrip(";").split(";"):
        part = part.strip()
        if part.startswith(key + " "):
            return part.split('"')[1] if '"' in part else part.split()[-1]
    return None


def count_peaks(peak_file):
    """Count non-comment, non-empty lines in a BED/narrowPeak file."""
    count = 0
    with open(peak_file) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def main():
    p = argparse.ArgumentParser("Intersect peaks with TE annotation")
    p.add_argument("--peak", required=True, help="narrowPeak file")
    p.add_argument("--te-gtf", required=True, help="TE annotation GTF")
    p.add_argument("--bedtools", default="bedtools", help="bedtools binary")
    p.add_argument("--overlap-bed", required=True, help="Output overlap BED")
    p.add_argument("--overlap-tsv", required=True, help="Output overlap count TSV")
    p.add_argument("--sample-id", required=True, help="Sample ID")
    args = p.parse_args()

    total_peaks = count_peaks(args.peak)
    logger.info(f"Total peaks: {total_peaks}")

    # Step 1: bedtools intersect -wa -wb
    cmd = [
        args.bedtools, "intersect",
        "-a", args.peak,
        "-b", args.te_gtf,
        "-wa", "-wb",
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"bedtools intersect failed: {result.stderr}")
        sys.exit(1)

    # Write overlap BED
    with open(args.overlap_bed, "w") as f:
        f.write(result.stdout)

    # Step 2: Count overlaps per TE class
    te_class_counter = Counter()
    peak_te_classes = {}  # peak_id -> set of te_classes (deduplicate)

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 15:
            continue
        # narrowPeak has 10 cols, GTF has 9 cols
        # peak fields: 0-9, GTF fields: 10-18
        peak_id = fields[3]
        gtf_attr = fields[18] if len(fields) > 18 else fields[-1]

        # Try class_id first, then family_id, then gene_id
        te_class = (
            parse_gtf_attribute(gtf_attr, "class_id")
            or parse_gtf_attribute(gtf_attr, "family_id")
            or parse_gtf_attribute(gtf_attr, "gene_id")
            or "Unknown"
        )

        if peak_id not in peak_te_classes:
            peak_te_classes[peak_id] = set()
        peak_te_classes[peak_id].add(te_class)

    # Count unique peaks per TE class
    for peak_id, classes in peak_te_classes.items():
        for te_class in classes:
            te_class_counter[te_class] += 1

    # Write overlap TSV
    with open(args.overlap_tsv, "w") as f:
        f.write("sample_id\tte_class\tpeak_count\ttotal_peaks\toverlap_ratio\n")
        for te_class, count in sorted(te_class_counter.items(), key=lambda x: -x[1]):
            ratio = count / total_peaks if total_peaks > 0 else 0.0
            f.write(f"{args.sample_id}\t{te_class}\t{count}\t{total_peaks}\t{ratio:.6f}\n")

    logger.info(f"TE classes found: {len(te_class_counter)}")
    for te_class, count in te_class_counter.most_common(10):
        ratio = count / total_peaks if total_peaks > 0 else 0.0
        logger.info(f"  {te_class}: {count}/{total_peaks} ({ratio:.2%})")


if __name__ == "__main__":
    main()
