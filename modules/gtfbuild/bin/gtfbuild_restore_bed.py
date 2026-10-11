#!/usr/bin/env python3
"""Restore source annotation columns to a paftools.js liftover BED.

paftools liftover emits BED6 (new coords, old-locus name, 0, strand) and
drops the source's annotation columns. The name column encodes the source
locus as ``chrom_start_end``, so each lifted row can be joined back to its
source row(s) to restore the full annotation.

Output columns: chrom, start, end, name(gene_id), score(.), strand,
family_id, class_id  (name is the BED name column consumed downstream).
A provenance header plus a column-name line ('#' comments) is written first.

Usage:
    python gtfbuild_restore_bed.py --lifted lifted.raw.bed \
        --source source.bed --out lifted.bed [--genome NAME]
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lifted", required=True, help="paftools.js liftover output")
    p.add_argument("--source", required=True, help="source BED (annotation cols)")
    p.add_argument("--out", required=True, help="restored output BED")
    p.add_argument("--genome", default="", help="target genome name (header)")
    args = p.parse_args()

    # key (chrom, start, end) -> list of source rows (rare duplicates kept)
    src = defaultdict(list)
    with open(args.source) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            src[(f[0], f[1], f[2])].append(f)

    n_in = n_out = n_miss = 0
    with open(args.lifted) as fh, open(args.out, "w") as out:
        # provenance header + column names ('#' comments; BED parsers skip them)
        if args.genome:
            out.write(f"#target-genome {args.genome}\n")
        out.write(f"#lifted-from {os.path.basename(args.source)} via WGA "
                  f"liftover (minimap2 asm5 + paftools.js liftover)\n")
        out.write("#note projected coordinates, not an official annotation\n")
        out.write("#chrom\tstart\tend\tname\tscore\tstrand\tfamily_id\tclass_id\n")
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            n_in += 1
            f = line.rstrip("\n").split("\t")
            chrom, start, end, name = f[0], f[1], f[2], f[3]
            key_parts = name.rsplit("_", 2)
            if len(key_parts) != 3:
                n_miss += 1
                continue
            rows = src.get((key_parts[0], key_parts[1], key_parts[2]))
            if not rows:
                n_miss += 1
                continue
            for r in rows:
                # source: chrom start end strand gene_id family_id class_id
                strand = r[3] if len(r) > 3 else "."
                gene_id = r[4] if len(r) > 4 else name
                family = r[5] if len(r) > 5 else "."
                cls = r[6] if len(r) > 6 else "."
                out.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t.\t"
                          f"{strand}\t{family}\t{cls}\n")
                n_out += 1
    print(f"[gtfbuild_restore_bed] lifted={n_in} written={n_out} "
          f"unjoined={n_miss}", file=sys.stderr)


if __name__ == "__main__":
    main()
