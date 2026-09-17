#!/usr/bin/env python3
"""RNA-seq expression normalization from featureCounts output.

Reads a featureCounts count matrix and computes normalized expression values
using one of several methods: CPM, RPKM/FPKM, TPM, or raw counts.

Optionally converts Ensembl gene IDs to gene symbols using a GTF annotation.

Normalization Methods:
    CPM  = (C_i / N) × 10^6          — Counts Per Million
    RPKM = C_i / (L_i/10^3 × N/10^6) — Reads Per Kilobase Million
    FPKM = same as RPKM (alias)
    TPM  = (C_i/L_i) / Σ(C_j/L_j) × 10^6 — Transcripts Per Million

Usage:
    python normalization.py -i counts.tsv -o output.tsv [--method tpm] [--gtf annot.gtf]
"""

import pandas as pd
import re
import sys
import os
import argparse
from pathlib import Path
from typing import Literal, Optional
try:
    from .LogUtil import setup_logger
    from .gene_id2name import convert_featurecounts_gene_ids
except ImportError:
    from LogUtil import setup_logger
    from gene_id2name import convert_featurecounts_gene_ids
logger = setup_logger(__name__)

class RNASeqNormalizer:
    """Normalize featureCounts output to CPM, RPKM, FPKM, or TPM.

    Parameters
    ----------
    gtf_path : str, optional
        Path to GTF annotation for gene_id → gene_name conversion.
        If not provided, gene IDs are used as-is.
    gene_id_col : str
        Column name for gene identifiers in the input file (default: "Geneid").
    """

    def __init__(self, gtf_path: Optional[str] = None, gene_id_col: str = "Geneid"):
        self.gtf_path = gtf_path
        self.gene_id_col = gene_id_col

    @staticmethod
    def compute_cpm(counts_df: pd.DataFrame, length: str = "Length") -> pd.DataFrame:
        """Counts Per Million: CPM_i = (C_i / N) × 10^6, where N = total library size."""
        counts_only = counts_df.drop(columns=[length])
        cpm = counts_only.div(counts_only.sum(axis=0), axis=1) * 1e6
        return cpm

    @staticmethod
    def compute_rpkm(counts_df: pd.DataFrame, length: str = "Length") -> pd.DataFrame:
        """Reads Per Kilobase Million: RPKM_i = C_i / (L_i/10^3 × N/10^6)."""
        counts_only = counts_df.drop(columns=[length])
        gene_length_kb = counts_df[length] / 1000
        rpkm = counts_only.div(gene_length_kb, axis=0)
        rpkm = rpkm.div(rpkm.sum(axis=0) / 1e6, axis=1)
        return rpkm

    @staticmethod
    def compute_tpm(counts_df: pd.DataFrame, length: str = "Length") -> pd.DataFrame:
        """Transcripts Per Million: TPM_i = (C_i/L_i) / Σ(C_j/L_j) × 10^6."""
        counts_only = counts_df.drop(columns=[length])
        gene_length_kb = counts_df[length] / 1000
        rpk = counts_only.div(gene_length_kb, axis=0)
        tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
        return tpm

    @staticmethod
    def extract_sample_name(
        col: str,
        pattern_aligned: str = r'([^/]+?)\.Aligned',
        pattern_fallback: str = r'/([^/]+)\.[^.]+$'
    ) -> str:
        """Extract sample name from featureCounts / BAM column names.

        Rules:
        1. If column contains 'Aligned', extract the part before 'Aligned'
        2. Otherwise, extract filename (without extension) from a file path
        3. If no pattern matches, return the original column name
        """
        m = re.search(pattern_aligned, col)
        if m:
            return m.group(1)
        m = re.search(pattern_fallback, col)
        if m:
            return m.group(1)
        return col

    def run_norm(
        self,
        infile: str,
        gtf: Optional[str] = None,
        method: Literal["cpm", "rpkm", "fpkm", "tpm", "count"] = "tpm",
        convert_to_gene_name: bool = True,
        remove_version: bool = True
    ) -> pd.DataFrame:
        """Normalize featureCounts output.

        Parameters
        ----------
        infile : str
            Path to featureCounts output TSV.
        method : str
            Normalization method: cpm, rpkm, fpkm, tpm, or count (no normalization).
        gtf : str, optional
            Path to GTF annotation (overrides constructor gtf_path).
        convert_to_gene_name : bool
            Convert Ensembl gene IDs to gene symbols using GTF.
        remove_version : bool
            Strip version suffix from gene IDs (e.g. ENSMUSG0000001 → ENSMUSG0000001).

        Returns
        -------
        pd.DataFrame
            Normalized expression matrix.
        """
        gene_id_col = self.gene_id_col
        logger.info(f"Normalization: method={method}, convert_to_gene_name={convert_to_gene_name}, remove_version={remove_version}")
        target_gtf = gtf or self.gtf_path

        df_counts = pd.read_csv(infile, sep="\t", comment='#')
        df_counts.drop(columns=['Chr', 'Start', 'End', 'Strand'], inplace=True, errors='ignore')
        df_counts = df_counts.set_index(gene_id_col)
        df_counts.columns = [self.extract_sample_name(c) for c in df_counts.columns]

        if method == "cpm":
            df = self.compute_cpm(df_counts, length="Length")
        elif method in ["rpkm", "fpkm"]:
            df = self.compute_rpkm(df_counts, length="Length")
        elif method == "tpm":
            df = self.compute_tpm(df_counts, length="Length")
        else:
            df = df_counts.drop(columns=["Length"], errors="ignore")
            logger.info("No normalization applied, only dropped Length column")

        df = df.reset_index()
        if convert_to_gene_name:
            if not target_gtf or not os.path.exists(target_gtf):
                raise ValueError(f"GTF file missing or not found: {target_gtf}")
            df = convert_featurecounts_gene_ids(df, target_gtf)
        else:
            if remove_version:
                df[gene_id_col] = df[gene_id_col].apply(lambda x: x.split('.')[0])

        return df

    def combine_PE_SE(self, PE: str, SE: str) -> pd.DataFrame:
        """Merge paired-end and single-end featureCounts output."""
        gene_id_col = self.gene_id_col
        df_PE = pd.read_csv(PE, sep="\t", comment='#')
        df_PE.columns = [self.extract_sample_name(c) for c in df_PE.columns]
        df_SE = pd.read_csv(SE, sep="\t", comment='#')
        df_SE.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'], inplace=True, errors='ignore')
        df_SE.columns = [self.extract_sample_name(c) for c in df_SE.columns]
        return pd.merge(df_PE, df_SE, on=gene_id_col)


def main() -> None:
    p = argparse.ArgumentParser(
        description="Normalize featureCounts output (CPM/RPKM/FPKM/TPM)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("-i", "--input", required=True,
                   help="featureCounts output file")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV path")
    p.add_argument("--method", default="tpm",
                   choices=["cpm", "rpkm", "fpkm", "tpm", "count"],
                   help="Normalization method (default: tpm)")
    p.add_argument("--gtf", default=None,
                   help="GTF annotation file for gene ID conversion")
    p.add_argument("--gene-id-col", default="Geneid",
                   help="Gene ID column name in input (default: Geneid)")
    p.add_argument("--convert-gene-name", action="store_true", default=False,
                   help="Convert gene IDs to gene names using GTF")
    p.add_argument("--no-remove-version", action="store_true", default=False,
                   help="Keep version suffix in gene IDs (e.g. .1, .2)")
    args = p.parse_args()

    normalizer = RNASeqNormalizer(gtf_path=args.gtf, gene_id_col=args.gene_id_col)
    df = normalizer.run_norm(
        args.input,
        method=args.method,
        convert_to_gene_name=args.convert_gene_name,
        remove_version=not args.no_remove_version,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    logger.info(f"Normalized expression saved to {args.output} ({df.shape[0]} genes × {df.shape[1] - 1} samples)")


if __name__ == "__main__":
    main()
