#!/usr/bin/env python3
"""Extract small RNA genes from a GENCODE GTF and produce BED + FASTA (±flank bp).

Usage:
    python extract_smallrna.py \
        --gtf gencode.v47.annotation.gtf \
        --fasta hg38.fa \
        --chrom-sizes chrom.sizes \
        --outdir output/genome/smallrna \
        --flank 50 \
        --types miRNA snRNA snoRNA rRNA misc_RNA scRNA scaRNA vaultRNA
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List


def parse_gtf_for_smallrna(gtf: Path, types: set) -> List[dict]:
    """Parse GTF and extract gene entries matching given gene_type."""
    genes = []
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = cols[8]
            gt = _attr(attrs, "gene_type")
            if gt not in types:
                continue
            genes.append({
                "chrom": cols[0],
                "start": int(cols[3]) - 1,  # BED is 0-based
                "end": int(cols[4]),
                "gene_id": _attr(attrs, "gene_id"),
                "score": ".",
                "strand": cols[6],
                "gene_name": _attr(attrs, "gene_name"),
            })
    return genes


def _attr(attrs: str, key: str) -> str:
    """Extract a single value from GTF attributes column."""
    for token in attrs.split(";"):
        token = token.strip()
        if token.startswith(key + " "):
            return token.split('"')[1] if '"' in token else token.split()[-1]
    return ""


def write_bed(genes: List[dict], out: Path):
    """Write BED6+1 format: chrom start end gene_id score strand gene_name."""
    with open(out, "w") as f:
        for g in genes:
            f.write(f"{g['chrom']}\t{g['start']}\t{g['end']}\t"
                    f"{g['gene_id']}\t{g['score']}\t{g['strand']}\t"
                    f"{g['gene_name']}\n")


def bed_to_fasta(bed: Path, fasta: Path, chrom_sizes: Path, flank: int, out: Path,
                 bedtools: str = "bedtools"):
    """BED → slop ±flank bp → getfasta → FASTA with gene names."""
    slop_cmd = [
        bedtools, "slop", "-i", str(bed), "-g", str(chrom_sizes),
        "-b", str(flank), "-s"
    ]
    getfasta_cmd = [
        bedtools, "getfasta", "-fi", str(fasta), "-bed", "stdin",
        "-s", "-name", "-fo", str(out)
    ]
    p1 = subprocess.Popen(slop_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(getfasta_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert p1.stdout is not None
    p1.stdout.close()
    _, err2 = p2.communicate()
    if p1.wait() != 0:
        err1 = p1.stderr.read().decode() if p1.stderr else ""
        print(f"[ERROR] bedtools slop failed: {err1}", file=sys.stderr)
        sys.exit(1)
    if p2.returncode != 0:
        print(f"[ERROR] bedtools getfasta failed: {err2.decode()}", file=sys.stderr)
        sys.exit(1)


def main():
    p = argparse.ArgumentParser("Extract small RNA genes from GENCODE GTF")
    p.add_argument("--gtf", required=True, help="GENCODE GTF file")
    p.add_argument("--fasta", required=True, help="Reference genome FASTA")
    p.add_argument("--chrom-sizes", required=True, help="chrom.sizes file")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--flank", type=int, default=50, help="Flanking bp (default: 50)")
    p.add_argument("--types", nargs="+",
                   default=["miRNA", "snRNA", "snoRNA", "rRNA", "misc_RNA", "scRNA", "scaRNA", "vaultRNA"],
                   help="Gene types to extract")
    p.add_argument("--bedtools", default="bedtools", help="bedtools binary path")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bed_out = outdir / "smallrna_genes.bed"
    fasta_out = outdir / "smallrna_genes_flank.fa"

    # 1. Parse GTF → BED
    types = set(args.types)
    print(f"[extract_smallrna] Filtering gene_type in {types}")
    genes = parse_gtf_for_smallrna(Path(args.gtf), types)
    print(f"[extract_smallrna] Found {len(genes)} small RNA genes")
    if not genes:
        print("[ERROR] No small RNA genes found. Check GTF and --types.", file=sys.stderr)
        sys.exit(1)
    write_bed(genes, bed_out)
    print(f"[extract_smallrna] BED → {bed_out}")

    # 2. BED → FASTA (±flank)
    bed_to_fasta(bed_out, Path(args.fasta), Path(args.chrom_sizes), args.flank, fasta_out, args.bedtools)
    print(f"[extract_smallrna] FASTA (±{args.flank}bp) → {fasta_out}")


if __name__ == "__main__":
    main()
