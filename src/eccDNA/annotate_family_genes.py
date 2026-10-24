#!/usr/bin/env python3
"""Annotate ecDNA breakpoint-family loci with genes (GRCm39 GENCODE).

For each family in the breakpoint_families TSV, lists all genes overlapping
the family's spanned locus, tagging each as ``full`` (gene entirely inside the
locus) or ``partial`` (gene extends beyond it). Writes a gene-level TSV.

Parameters
----------
--families-tsv : TSV from breakpoint_families.py
--ecdna-dir : directory containing <sample>/decoil/*filtered.bed (spans)
--samples : sample names, repeatable
--gtf : GENCODE GTF (GRCm39)
--out : output TSV path
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import subprocess

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
logger = logging.getLogger("annotate_family_genes")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    p.add_argument("--families-tsv", required=True)
    p.add_argument("--ecdna-dir", required=True)
    p.add_argument("--samples", action="append", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args(argv)


def family_spans(families_tsv, ecdna_dir, samples):
    circles = {}
    for s in samples:
        path = os.path.join(ecdna_dir, s, "decoil",
                            "reconstruct.ecDNA.filtered.bed")
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.split("\t")
                circles.setdefault((s, f[3]), []).append(
                    (f[0], int(f[1]), int(f[2])))
    spans = {}
    with open(families_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            members = [(m.split(":circ")[0], m.split(":circ")[1])
                       for m in row["members"].split(";")]
            frags = [f for m in members for f in circles[m]]
            spans[row["family"]] = (frags[0][0],
                                    min(f[1] for f in frags),
                                    max(f[2] for f in frags))
    return spans


def genes_in(chrom, lo, hi, gtf):
    """Genes overlapping [lo, hi]; returns (name, start, end, full/partial)."""
    out = subprocess.run(
        ["bash", "-c",
         f"awk -F'\\t' -v C={chrom} -v L={lo} -v H={hi} "
         f"'$1==C && $3==\"gene\" && $5>=L && $4<=H' {gtf}"],
        capture_output=True, text=True)
    genes = []
    for line in out.stdout.splitlines():
        f = line.split("\t")
        m = re.search(r'gene_name "([^"]+)"', f[8])
        if m:
            name = m.group(1)
        else:
            name = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        gs, ge = int(f[3]), int(f[4])
        tag = "full" if (gs >= lo and ge <= hi) else "partial"
        genes.append((name, gs, ge, tag))
    return genes


def structures_for(chrom, lo, hi, gtf):
    """Exon/UTR blocks of genes overlapping [lo, hi], one row per block.

    Returns list of (gene, block_type, start, end, strand); blocks are the
    union across transcripts (exons merged), plus UTR blocks as-is.
    """
    out = subprocess.run(
        ["bash", "-c",
         f"awk -F'\\t' -v C={chrom} -v L={lo} -v H={hi} "
         f"'$1==C && ($3==\"exon\" || $3==\"UTR\") "
         f"&& $5>=L && $4<=H' {gtf}"],
        capture_output=True, text=True)
    blocks = []
    for line in out.stdout.splitlines():
        f = line.split("	")
        m = re.search(r'gene_name "([^"]+)"', f[8])
        if not m:
            continue
        blocks.append((m.group(1), f[2], int(f[3]), int(f[4]), f[6]))
    return blocks


def merge_exons(blocks):
    """Merge exon rows per gene into canonical intron-exon structure."""
    per = {}
    for gene, btype, start, end, strand in blocks:
        per.setdefault(gene, {"exons": [], "UTRs": [], "strand": strand})
        key = "UTRs" if btype == "UTR" else "exons"
        per[gene][key].append((start, end))
    for gene, d in per.items():
        ex = sorted(d["exons"])
        merged = []
        for s, e in ex:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        d["exons"] = merged
    return per


def main(argv=None):
    args = parse_args(argv)
    spans = family_spans(args.families_tsv, args.ecdna_dir, args.samples)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    struct_path = os.path.splitext(args.out)[0] + "_structures.tsv"
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family", "locus", "gene", "gene_start", "gene_end",
                    "containment"])
        n_blocks = 0
        with open(struct_path, "w", newline="") as fs:
            ws = csv.writer(fs, delimiter="\t")
            ws.writerow(["family", "gene", "block_type", "start", "end",
                         "strand"])
            for fam in sorted(spans, key=int):
                chrom, lo, hi = spans[fam]
                genes = genes_in(chrom, lo, hi, args.gtf)
                logger.info("family %s %s:%s-%s -> %d genes",
                            fam, chrom, f"{lo:,}", f"{hi:,}", len(genes))
                for name, gs, ge, tag in genes:
                    w.writerow([fam, f"{chrom}:{lo:,}-{hi:,}", name, gs, ge,
                                tag])
                blocks = merge_exons(structures_for(chrom, lo, hi, args.gtf))
                for gene, d in blocks.items():
                    for btype in ("exons", "UTRs"):
                        for s, e in d[btype]:
                            ws.writerow([fam, gene, btype, s, e, d["strand"]])
                            n_blocks += 1
    logger.info("Wrote %s (%d structure rows) and %s",
                args.out, n_blocks, struct_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
