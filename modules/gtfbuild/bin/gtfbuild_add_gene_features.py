#!/usr/bin/env python3
"""Add gene-level features to a GTF that only has transcript/exon/CDS lines.

gffread -T converts GFF3 to GTF but drops gene features; downstream tools
(fire_annotate, nucpos) classify peaks against ``gene`` + ``exon`` lines.
This script derives one gene line per gene_id from the transcript spans.

Usage:
    python gtfbuild_add_gene_features.py --gtf in.gtf --out out.gtf
"""
from __future__ import annotations

import argparse
import sys


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", required=True, help="input GTF (no gene features)")
    p.add_argument("--out", required=True, help="output GTF with gene features")
    args = p.parse_args()

    # Pass 1: gene spans + names
    spans = {}  # gene_id -> [chrom, strand, start, end, gene_name, source]
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "mRNA"):
                continue
            attrs = {}
            for kv in f[8].strip().rstrip(";").split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                parts = kv.split(None, 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip('"')
            gid = attrs.get("gene_id")
            if not gid:
                continue
            start, end = int(f[3]), int(f[4])
            if gid not in spans:
                spans[gid] = [f[0], f[6], start, end,
                              attrs.get("gene_name", gid), f[1]]
            else:
                s = spans[gid]
                s[2] = min(s[2], start)
                s[3] = max(s[3], end)

    # Pass 2: emit one gene line before the first transcript of each gene
    seen = set()
    n_gene = 0
    with open(args.gtf) as fh, open(args.out, "w") as out:
        for line in fh:
            if not line.startswith("#"):
                f = line.split("\t")
                if len(f) >= 9 and f[2] in ("transcript", "mRNA"):
                    gid = None
                    for kv in f[8].strip().rstrip(";").split(";"):
                        kv = kv.strip()
                        if kv.startswith("gene_id"):
                            gid = kv.split(None, 1)[1].strip().strip('"')
                    if gid and gid not in seen:
                        seen.add(gid)
                        chrom, strand, start, end, gname, source = spans[gid]
                        out.write(
                            f'{chrom}\t{source}\tgene\t{start}\t{end}\t.\t'
                            f'{strand}\t.\tgene_id "{gid}"; gene_name "{gname}";\n'
                        )
                        n_gene += 1
            out.write(line)
    print(f"[gtfbuild_add_gene_features] wrote {n_gene} gene features to {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
