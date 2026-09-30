#!/usr/bin/env python3
"""Annotate h5ad var with gene_type, family_id, class_id from TE BED + gene TSV.

Annotates both .var and .raw.var (if raw exists).
Writes to a new file first for verification; --replace overwrites original.

Usage:
    python annotate_gene_type.py \
        --h5ad ovaries_scTE_auto.h5ad \
        --te-bed rheMac10_rmsk_TE.bed \
        --gene-tsv geneIDAnnotation.csv \
        [--replace]
"""
import argparse
import logging
import os
import shutil
import sys
from typing import Dict, List, Tuple

import anndata as ad
import pandas as pd

ad.settings.allow_write_nullable_strings = True

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def parse_te_bed(bed_path: str) -> Dict[str, Tuple[str, str]]:
    """TE BED: chrom, start, end, strand, gene_id, family_id, class_id (with header)."""
    result: Dict[str, Tuple[str, str]] = {}
    with open(bed_path, encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            gene_id = parts[4]
            if gene_id not in result:
                result[gene_id] = (parts[5], parts[6])
    return result


def parse_gene_tsv(tsv_path: str) -> Dict[str, str]:
    """Gene TSV: gene_id, gene_name, gene_type (with header)."""
    result: Dict[str, str] = {}
    with open(tsv_path, encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            gene_name = parts[1].strip()
            gene_type = parts[2].strip()
            if gene_name:
                result[gene_name] = gene_type
    return result


def annotate_var(var: pd.DataFrame,
                 te_map: Dict[str, Tuple[str, str]],
                 gene_map: Dict[str, str]) -> pd.DataFrame:
    """Annotate a var DataFrame with gene_type, family_id, class_id."""
    types: List[str] = []
    families: List[str] = []
    classes: List[str] = []
    n_te = n_gene = n_unknown = 0

    for g in var.index:
        if g in te_map:
            types.append("TE")
            fam, cls = te_map[g]
            families.append(fam)
            classes.append(cls)
            n_te += 1
        elif g in gene_map:
            types.append(gene_map[g] or "unknown")
            families.append("")
            classes.append("")
            n_gene += 1
        else:
            types.append("unknown")
            families.append("")
            classes.append("")
            n_unknown += 1

    var["gene_type"] = pd.Categorical(types)
    var["family_id"] = pd.Categorical(families)
    var["class_id"] = pd.Categorical(classes)
    logger.info("  TE=%d, gene=%d, unknown=%d (total=%d)",
                n_te, n_gene, n_unknown, len(var))
    return var


def main():
    parser = argparse.ArgumentParser(description="Annotate h5ad with TE family/class")
    parser.add_argument("--h5ad", required=True, help="Input h5ad file")
    parser.add_argument("--te-bed", required=True, help="TE BED file")
    parser.add_argument("--gene-tsv", required=True, help="Gene annotation TSV")
    parser.add_argument("--replace", action="store_true",
                        help="Overwrite original file (default: write to .annotated.h5ad)")
    args = parser.parse_args()

    # Parse references
    te_map = parse_te_bed(args.te_bed)
    logger.info("TE BED: %d unique TE subfamilies", len(te_map))
    gene_map = parse_gene_tsv(args.gene_tsv)
    logger.info("Gene TSV: %d entries", len(gene_map))

    # Read h5ad
    logger.info("Reading: %s", args.h5ad)
    adata = ad.read_h5ad(args.h5ad)
    logger.info("Shape: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Annotate main .var
    logger.info("Annotating .var ...")
    adata.var = annotate_var(adata.var, te_map, gene_map)

    # Annotate .raw.var if exists
    if adata.raw is not None:
        logger.info("Annotating .raw.var ...")
        raw_adata = adata.raw.to_adata()
        raw_adata.var = annotate_var(raw_adata.var, te_map, gene_map)
        adata.raw = raw_adata

    # Output
    if args.replace:
        out_path = args.h5ad
        # Backup original
        backup = args.h5ad + ".bak"
        if not os.path.exists(backup):
            shutil.copy2(args.h5ad, backup)
            logger.info("Backup: %s", backup)
    else:
        base, ext = os.path.splitext(args.h5ad)
        out_path = f"{base}.annotated{ext}"

    logger.info("Writing: %s", out_path)
    adata.write_h5ad(out_path)
    logger.info("Done.")

    # Summary
    vc = adata.var["gene_type"].value_counts()
    logger.info("gene_type summary:\n%s", vc.to_string())
    te_count = int((adata.var["gene_type"] == "TE").sum())
    logger.info("TE with family_id: %d",
                (adata.var.loc[adata.var["gene_type"] == "TE", "family_id"] != "").sum())
    logger.info("TE with class_id: %d",
                (adata.var.loc[adata.var["gene_type"] == "TE", "class_id"] != "").sum())


if __name__ == "__main__":
    main()
