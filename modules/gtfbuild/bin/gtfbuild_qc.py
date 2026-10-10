#!/usr/bin/env python3
"""QC for LiftOff-based GTF projection (gtfbuild module).

Compares the lifted GFF3/GTF against the source GTF and reports mapping
completeness: gene/transcript mapping rate, partial models, extra copies.

Usage:
    python gtfbuild_qc.py --src-gtf gencode.gtf --lifted-gff3 lifted.gff3 \
        --lifted-gtf lifted.gtf --unmapped unmapped.gff3 \
        --genome T2T_mhaESC_v1.5 --out qc.tsv
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict


def _get_logger() -> logging.Logger:
    logger = logging.getLogger("gtfbuild_qc")
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)s | %(name)s | %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


logger = _get_logger()


def parse_attrs(field: str) -> dict:
    """Parse GTF (key "value";) or GFF3 (key=value;) attributes."""
    attrs = {}
    for kv in field.strip().rstrip(";").split(";"):
        kv = kv.strip()
        if not kv:
            continue
        if "=" in kv and '"' not in kv.split("=", 1)[1]:
            k, v = kv.split("=", 1)
            attrs[k.strip()] = v.strip()
        elif " " in kv:
            k, v = kv.split(None, 1)
            attrs[k.strip()] = v.strip().strip('"')
    return attrs


def count_features(path: str) -> dict:
    """Count gene/transcript/exon features and LiftOff mapping attributes."""
    stats = {
        "gene_ids": set(),
        "transcripts": 0,
        "exons": 0,
        "partial_transcripts": 0,
        "extra_copy_transcripts": 0,
    }
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feature, attrs = f[2], parse_attrs(f[8])
            gid = attrs.get("gene_id") or attrs.get("gene_name")
            if feature == "gene" and gid:
                stats["gene_ids"].add(gid)
            elif feature in ("transcript", "mRNA"):
                stats["transcripts"] += 1
                if attrs.get("partial_mapping") in ("True", "true", "1"):
                    stats["partial_transcripts"] += 1
                try:
                    if int(attrs.get("extra_copy_number", 0)) > 0:
                        stats["extra_copy_transcripts"] += 1
                except ValueError:
                    pass
            elif feature == "exon":
                stats["exons"] += 1
    return stats


def count_unmapped(path: str) -> int:
    n = 0
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                n += 1
    except FileNotFoundError:
        logger.warning(f"unmapped file not found: {path}")
    return n


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--src-gtf", required=True, help="source annotation GTF")
    p.add_argument("--lifted-gff3", required=True, help="LiftOff GFF3 output")
    p.add_argument("--lifted-gtf", required=True, help="converted lifted GTF")
    p.add_argument("--unmapped", default="", help="LiftOff -unmapped output")
    p.add_argument("--genome", default="genome", help="target genome name")
    p.add_argument("--out", required=True, help="output QC TSV")
    args = p.parse_args()

    src = count_features(args.src_gtf)
    lifted = count_features(args.lifted_gff3)
    unmapped = count_unmapped(args.unmapped) if args.unmapped else 0

    src_genes = len(src["gene_ids"])
    lifted_genes = len(lifted["gene_ids"])
    rate = (lifted_genes / src_genes) if src_genes else 0.0

    rows = [
        ("genome", args.genome),
        ("source_gtf", args.src_gtf),
        ("lifted_gff3", args.lifted_gff3),
        ("source_genes", src_genes),
        ("source_transcripts", src["transcripts"]),
        ("source_exons", src["exons"]),
        ("lifted_genes", lifted_genes),
        ("lifted_transcripts", lifted["transcripts"]),
        ("lifted_exons", lifted["exons"]),
        ("unmapped_features", unmapped),
        ("gene_mapping_rate", f"{rate:.4f}"),
        ("partial_transcripts", lifted["partial_transcripts"]),
        ("extra_copy_transcripts", lifted["extra_copy_transcripts"]),
    ]
    with open(args.out, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")
    logger.info(f"QC written to {args.out}: gene_mapping_rate={rate:.4f} "
                f"({lifted_genes}/{src_genes} genes)")


if __name__ == "__main__":
    main()
