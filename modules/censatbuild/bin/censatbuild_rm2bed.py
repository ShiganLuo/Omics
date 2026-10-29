#!/usr/bin/env python3
"""Convert RepeatMasker .out satellite hits to clustered array BED (censatbuild).

Parses the RepeatMasker tabular .out, classifies each hit by the consensus
name (major/minor satellite patterns), clusters collinear hits into arrays
(gap <= --gap), and writes:
  - arrays BED: chrom, start, end, class, hit_count, .   (name column = class)
  - hits BED:   raw per-monomer hits with class in name column

Region classes feed modules/fiberseq/censat (name column = region class).

Usage:
    python censatbuild_rm2bed.py --rm-out genome.rm.out \
        --arrays-out arrays.bed --hits-out hits.bed \
        [--gap 10000] [--major-pattern major] [--minor-pattern minor]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict


def classify(name: str, major_pat: str, minor_pat: str) -> str:
    low = name.lower()
    if minor_pat and minor_pat.lower() in low:
        return "minorSat_core"
    if major_pat and major_pat.lower() in low:
        return "majorSat_array"
    return "otherSat"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rm-out", required=True, help="RepeatMasker .out file")
    p.add_argument("--arrays-out", required=True, help="clustered arrays BED")
    p.add_argument("--hits-out", required=True, help="raw hits BED")
    p.add_argument("--gap", type=int, default=10000, help="max gap within an array (bp)")
    p.add_argument("--major-pattern", default="major", help="consensus name pattern: major satellite")
    p.add_argument("--minor-pattern", default="minor", help="consensus name pattern: minor satellite")
    p.add_argument("--class-pattern", default="Satellite",
                   help="keep only hits whose repeat class/family contains this (empty = keep all)")
    args = p.parse_args()

    hits = []  # (chrom, start, end, class)
    with open(args.rm_out) as fh:
        in_table = False
        for line in fh:
            if line.lstrip().startswith("SW") and "perc" in line:
                in_table = True
                continue
            if not in_table or line.startswith("----") or not line.strip():
                continue
            f = line.split()
            if len(f) < 11:
                continue
            try:
                start, end = int(f[5]), int(f[6])
            except ValueError:
                continue
            if args.class_pattern and \
                    args.class_pattern.lower() not in f[10].lower():
                continue
            chrom = f[4]
            cls = classify(f[9], args.major_pattern, args.minor_pattern)
            if start > end:  # minus strand entries are still query coords
                start, end = end, start
            hits.append((chrom, start - 1, end, cls))  # to BED 0-based

    with open(args.hits_out, "w") as out:
        for chrom, start, end, cls in hits:
            out.write(f"{chrom}\t{start}\t{end}\t{cls}\t.\t.\n")

    # cluster per (chrom, class)
    by_key = defaultdict(list)
    for chrom, start, end, cls in hits:
        by_key[(chrom, cls)].append((start, end))

    n_arrays = 0
    with open(args.arrays_out, "w") as out:
        for (chrom, cls), spans in sorted(by_key.items()):
            spans.sort()
            cs, ce, cnt = spans[0][0], spans[0][1], 1
            for s, e in spans[1:]:
                if s - ce <= args.gap:
                    ce = max(ce, e)
                    cnt += 1
                else:
                    out.write(f"{chrom}\t{cs}\t{ce}\t{cls}\t{cnt}\t.\n")
                    n_arrays += 1
                    cs, ce, cnt = s, e, 1
            out.write(f"{chrom}\t{cs}\t{ce}\t{cls}\t{cnt}\t.\n")
            n_arrays += 1

    print(f"[censatbuild_rm2bed] hits={len(hits)} arrays={n_arrays} "
          f"gap={args.gap}", file=sys.stderr)


if __name__ == "__main__":
    main()
