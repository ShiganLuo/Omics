#!/usr/bin/env python3
"""Find overlaps between narrowPeak and TE annotation GTF.

Uses bedtools intersect to find peak-TE overlaps, then produces
per-locus overlap data per TE subfamily (gene_id) with TE class info.

Also counts reads from IP and Input BAMs in overlap regions
for enrichment analysis.

Outputs:
  - overlap BED: full intersect result
  - subfamily TSV: per-locus overlap fraction, peak count, TE class, and read counts
"""
import argparse
import logging
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

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


def run_bedtools_coverage(bed_file, bam_file, bedtools="bedtools"):
    """Count reads in BED regions from a BAM file using bedtools coverage -counts.

    Filters out PCR duplicates (flag 0x400) via samtools view -F 1024.
    Returns dict: (chrom, start, end) -> read_count
    """
    import tempfile as _tempfile
    tmp_bam = _tempfile.mktemp(suffix=".bam")
    try:
        # Step 1: filter duplicates with samtools
        samtools_cmd = ["samtools", "view", "-b", "-F", "1024", "-o", tmp_bam, bam_file]
        logger.info(f"Running: {' '.join(samtools_cmd)}")
        r = subprocess.run(samtools_cmd, capture_output=True, text=True)
        if r.returncode != 0:
            logger.error(f"samtools filter failed: {r.stderr}")
            return {}

        # Step 2: count reads with bedtools coverage
        bedtools_cmd = [
            bedtools, "coverage",
            "-a", bed_file,
            "-b", tmp_bam,
            "-counts",
        ]
        logger.info(f"Running: {' '.join(bedtools_cmd)}")
        result = subprocess.run(bedtools_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"bedtools coverage failed: {result.stderr}")
            return {}

        counts = {}
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) >= 4:
                key = (fields[0], int(fields[1]), int(fields[2]))
                counts[key] = int(fields[3])
        return counts
    finally:
        if os.path.exists(tmp_bam):
            os.unlink(tmp_bam)


def main():
    p = argparse.ArgumentParser("Intersect peaks with TE annotation")
    p.add_argument("--peak", required=True, help="narrowPeak file")
    p.add_argument("--te-gtf", required=True, help="TE annotation GTF")
    p.add_argument("--bedtools", default="bedtools", help="bedtools binary")
    p.add_argument("--overlap-bed", required=True, help="Output overlap BED")
    p.add_argument("--subfamily-tsv", required=True, help="Output subfamily locus-level TSV")
    p.add_argument("--sample-id", required=True, help="Sample ID")
    p.add_argument("--ip-bam", default="", help="IP BAM file for read counting")
    p.add_argument("--input-bam", default="", help="Input BAM file for read counting")
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

    # Step 2: Parse overlaps and collect subfamily locus data
    subfamily_loci = []   # (sample_id, subfamily, te_class, te_length, overlap_start, overlap_end, overlap_frac)
    subfamily_peak_ids = defaultdict(set)  # subfamily -> set of unique peak_ids
    overlap_regions = []  # (chrom, start, end) for bedtools coverage

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 15:
            continue
        peak_start = int(fields[1])
        peak_end = int(fields[2])
        peak_id = fields[3]
        gtf_chrom = fields[10]
        gtf_start = int(fields[13])
        gtf_end = int(fields[14])
        gtf_attr = fields[18] if len(fields) > 18 else fields[-1]

        te_class = (
            parse_gtf_attribute(gtf_attr, "class_id")
            or parse_gtf_attribute(gtf_attr, "family_id")
            or parse_gtf_attribute(gtf_attr, "gene_id")
            or "Unknown"
        )
        te_subfamily = (
            parse_gtf_attribute(gtf_attr, "gene_id")
            or parse_gtf_attribute(gtf_attr, "family_id")
            or "Unknown"
        )

        te_length = gtf_end - gtf_start
        overlap_start = max(peak_start, gtf_start)
        overlap_end = min(peak_end, gtf_end)
        overlap_length = max(0, overlap_end - overlap_start)
        overlap_frac = overlap_length / te_length if te_length > 0 else 0.0

        subfamily_loci.append((
            args.sample_id, te_subfamily, te_class, te_length,
            overlap_start, overlap_end, overlap_frac,
        ))
        subfamily_peak_ids[te_subfamily].add(peak_id)
        overlap_regions.append((gtf_chrom, overlap_start, overlap_end))

    # Step 3: Count reads in overlap regions from IP and Input BAMs
    ip_read_counts = {}  # (chrom, start, end) -> count
    input_read_counts = {}

    if args.ip_bam and args.input_bam:
        # Write overlap regions as temporary BED for bedtools coverage
        tmpdir = tempfile.mkdtemp()
        regions_bed = os.path.join(tmpdir, "overlap_regions.bed")
        with open(regions_bed, "w") as f:
            for chrom, start, end in overlap_regions:
                f.write(f"{chrom}\t{start}\t{end}\n")

        logger.info(f"Counting reads in {len(overlap_regions)} overlap regions")
        ip_read_counts = run_bedtools_coverage(regions_bed, args.ip_bam, args.bedtools)
        input_read_counts = run_bedtools_coverage(regions_bed, args.input_bam, args.bedtools)

        # Cleanup
        os.unlink(regions_bed)
        os.rmdir(tmpdir)
    else:
        logger.info("No BAM files provided, skipping read counting")

    # Step 4: Write subfamily locus-level TSV with TE class and read counts
    with open(args.subfamily_tsv, "w") as f:
        header = "sample_id\tte_subfamily\tte_class\tte_length\tinterval_overlap_frac\toverlap_peak_count"
        if ip_read_counts:
            header += "\tip_reads\tinput_reads"
        f.write(header + "\n")

        for sid, sfam, tclass, tlen, ostart, oend, ofrac in subfamily_loci:
            peak_count = len(subfamily_peak_ids.get(sfam, set()))
            line = f"{sid}\t{sfam}\t{tclass}\t{tlen}\t{ofrac:.6f}\t{peak_count}"
            if ip_read_counts:
                ip_reads = 0
                input_reads = 0
                for (chrom, rs, re), count in ip_read_counts.items():
                    if rs == ostart and re == oend:
                        ip_reads = count
                        break
                for (chrom, rs, re), count in input_read_counts.items():
                    if rs == ostart and re == oend:
                        input_reads = count
                        break
                line += f"\t{ip_reads}\t{input_reads}"
            f.write(line + "\n")

    subfamilies = set(sfam for _, sfam, *_ in subfamily_loci)
    logger.info(f"TE subfamilies found: {len(subfamilies)}, total loci: {len(subfamily_loci)}")


if __name__ == "__main__":
    main()
