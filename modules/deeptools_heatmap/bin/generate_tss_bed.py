#!/usr/bin/env python3
"""Generate a BED file of TSS regions from a GTF annotation.

Each line represents a transcription start site ± flank bp window.
Output: chr, start, end, gene_name, strand (tab-separated, BED6).
"""
import argparse
import logging
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("generate_tss_bed")


def parse_gtf_attribute(attr_str, key):
    """Extract a value from GTF attributes string."""
    for part in attr_str.rstrip(";").split(";"):
        part = part.strip()
        if part.startswith(key + " "):
            return part.split('"')[1] if '"' in part else part.split()[-1]
    return None


def main():
    p = argparse.ArgumentParser("Generate TSS BED from GTF")
    p.add_argument("--gtf", required=True, help="Input GTF annotation file")
    p.add_argument("--output", required=True, help="Output BED file path")
    p.add_argument("--flank", type=int, default=1000, help="Flanking region around TSS (bp)")
    p.add_argument("--feature", default="transcript", help="GTF feature type to extract (default: transcript)")
    args = p.parse_args()

    seen = set()
    count = 0

    with open(args.gtf) as gtf, open(args.output, "w") as out:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != args.feature:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2]) if len(fields) > 2 else start
            strand = fields[6]
            attrs = fields[8]

            gene_name = parse_gtf_attribute(attrs, "gene_name") or parse_gtf_attribute(attrs, "gene_id")
            if not gene_name:
                continue

            # TSS is start for + strand, end for - strand
            if strand == "+":
                tss = start
            else:
                tss = end

            bed_start = max(0, tss - args.flank)
            bed_end = tss + args.flank

            key = (chrom, bed_start, bed_end)
            if key in seen:
                continue
            seen.add(key)

            out.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_name}\t0\t{strand}\n")
            count += 1

    logger.info(f"Generated {count} TSS regions (±{args.flank}bp) from {args.gtf}")


if __name__ == "__main__":
    main()
