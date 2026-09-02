#!/usr/bin/env python3
"""Peak-centric TE overlap analysis.

For each peak, calculate:
  - te_count: number of TEs overlapping this peak
  - te_classes: TE class types (SINE, LINE, LTR, etc.)
  - te_class_count: number of distinct TE classes
  - te_covered_bases: total bases of the peak covered by TEs (merged)
  - peak_length: total peak length
  - te_coverage_frac: fraction of peak covered by TEs

Usage:
    python peak_centric_overlap.py \\
        --peak peaks.narrowPeak --te-gtf TE.gtf \\
        --output peak_te_summary.tsv
"""
import argparse
import logging
import subprocess
import sys
import tempfile
import os
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("peak_centric_overlap")


def parse_gtf_attribute(attr_str, key):
    """Extract a value from GTF attributes string."""
    for part in attr_str.rstrip(";").split(";"):
        part = part.strip()
        if part.startswith(key + " "):
            return part.split('"')[1] if '"' in part else part.split()[-1]
    return None


def merge_intervals(intervals):
    """Merge overlapping intervals. Input: list of (start, end). Returns merged list."""
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def main():
    p = argparse.ArgumentParser("Peak-centric TE overlap analysis")
    p.add_argument("--peak", required=True, help="narrowPeak file")
    p.add_argument("--te-gtf", required=True, help="TE annotation GTF")
    p.add_argument("--output", required=True, help="Output TSV file")
    p.add_argument("--bedtools", default="bedtools", help="bedtools binary")
    p.add_argument("--sample-id", default="", help="Sample ID to include in output")
    args = p.parse_args()

    # Step 1: bedtools intersect -wa -wb
    cmd = [args.bedtools, "intersect", "-a", args.peak, "-b", args.te_gtf, "-wa", "-wb"]
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"bedtools intersect failed: {result.stderr}")
        sys.exit(1)

    # Step 2: Parse peak info from narrowPeak
    peak_info = {}  # peak_id -> (chrom, start, end)
    with open(args.peak) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 4:
                peak_info[fields[3]] = (fields[0], int(fields[1]), int(fields[2]))

    # Step 3: Parse overlaps — group by peak
    peak_te_intervals = defaultdict(list)  # peak_id -> [(te_start, te_end)]
    peak_te_classes = defaultdict(set)     # peak_id -> {class1, class2, ...}
    peak_te_subfamilies = defaultdict(set) # peak_id -> {subfamily1, ...}

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 15:
            continue
        peak_id = fields[3]
        gtf_start = int(fields[13])
        gtf_end = int(fields[14])
        gtf_attr = fields[18] if len(fields) > 18 else fields[-1]

        te_class = (
            parse_gtf_attribute(gtf_attr, "class_id")
            or parse_gtf_attribute(gtf_attr, "family_id")
            or "Unknown"
        )
        te_subfamily = (
            parse_gtf_attribute(gtf_attr, "gene_id")
            or parse_gtf_attribute(gtf_attr, "family_id")
            or "Unknown"
        )

        peak_te_intervals[peak_id].append((gtf_start, gtf_end))
        peak_te_classes[peak_id].add(te_class)
        peak_te_subfamilies[peak_id].add(te_subfamily)

    # Step 4: Calculate peak-centric metrics
    rows = []
    for peak_id, (pchrom, pstart, pend) in peak_info.items():
        peak_length = pend - pstart
        te_count = len(peak_te_intervals.get(peak_id, []))
        classes = sorted(peak_te_classes.get(peak_id, set()))
        subfamilies = sorted(peak_te_subfamilies.get(peak_id, set()))

        # Merge TE intervals within this peak to compute covered bases
        intervals = peak_te_intervals.get(peak_id, [])
        if intervals:
            # Clip intervals to peak boundaries
            clipped = [(max(s, pstart), min(e, pend)) for s, e in intervals if max(s, pstart) < min(e, pend)]
            merged = merge_intervals(clipped)
            covered = sum(e - s for s, e in merged)
            coverage_frac = covered / peak_length if peak_length > 0 else 0.0
        else:
            covered = 0
            coverage_frac = 0.0

        rows.append({
            "peak_id": peak_id,
            "chrom": pchrom,
            "start": pstart,
            "end": pend,
            "peak_length": peak_length,
            "te_count": te_count,
            "te_class_count": len(classes),
            "te_classes": ",".join(classes) if classes else "none",
            "te_subfamily_count": len(subfamilies),
            "te_subfamilies": ",".join(subfamilies) if subfamilies else "none",
            "te_covered_bases": covered,
            "te_coverage_frac": f"{coverage_frac:.6f}",
        })

    # Step 5: Write output
    header = ["sample_id", "peak_id", "chrom", "start", "end", "peak_length",
              "te_count", "te_class_count", "te_classes",
              "te_subfamily_count", "te_subfamilies",
              "te_covered_bases", "te_coverage_frac"]
    with open(args.output, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            vals = [args.sample_id, r["peak_id"], r["chrom"], r["start"], r["end"],
                    r["peak_length"], r["te_count"], r["te_class_count"], r["te_classes"],
                    r["te_subfamily_count"], r["te_subfamilies"],
                    r["te_covered_bases"], r["te_coverage_frac"]]
            f.write("\t".join(str(v) for v in vals) + "\n")

    total = len(rows)
    with_te = sum(1 for r in rows if r["te_count"] > 0)
    logger.info(f"Total peaks: {total}, peaks with TE overlap: {with_te} ({with_te/total*100:.1f}%)" if total else "No peaks found")


if __name__ == "__main__":
    main()
