import pandas as pd
import re
import sys
import os
import logging
from pathlib import Path
from typing import Literal
current_path = Path(__file__).resolve()
scripts_dir = current_path.parents[1]
sys.path.append(str(scripts_dir))

try:
    from annotation.gene_id2name import convert_featurecounts_gene_ids
except ImportError:
    pass

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RNASeqNormalizer:
    def __init__(self, gtf_path=None):
        self.gtf_path = gtf_path

    @staticmethod
    def compute_cpm(counts_df: pd.DataFrame, length: str = "Length"):
        """
        Formula:
        $$CPM = \frac{C_i}{N} \cdot 10^6$$
        where $C_i$ is the count of gene $i$ and $N$ is the total library size.
        """
        counts_only = counts_df.drop(columns=[length])
        cpm = counts_only.div(counts_only.sum(axis=0), axis=1) * 1e6
        return cpm

    @staticmethod
    def compute_rpkm(counts_df: pd.DataFrame, length: str = "Length"):
        """
        Formula:
        $$RPKM = \frac{C_i}{\frac{L_i}{10^3} \cdot \frac{N}{10^6}}$$
        where $L_i$ is gene length in bp.
        """
        counts_only = counts_df.drop(columns=[length])
        gene_length_kb = counts_df[length] / 1000
        rpkm = counts_only.div(gene_length_kb, axis=0)
        rpkm = rpkm.div(rpkm.sum(axis=0) / 1e6, axis=1)
        return rpkm

    @staticmethod
    def compute_tpm(counts_df: pd.DataFrame, length: str = "Length"):
        """
        Formula:
        $$TPM_i = \frac{rpk_i}{\sum rpk} \cdot 10^6$$
        where $rpk_i = \frac{C_i}{L_i/10^3}$.
        """
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
    ):
        """
        Extract sample name from featureCounts / BAM column names.

        Rules:
        1. If the column name contains 'Aligned', extract the part before 'Aligned'
        2. Otherwise, try to extract the filename (without extension) from a file path
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
        gtf: str = None,
        method: Literal["cpm", "rpkm", "fpkm", "tpm", "count"] = "tpm",
        convert_to_gene_name: bool = True,
        remove_version: bool = True
    ) -> pd.DataFrame:
        """
        Perform gene expression normalization on featureCounts output.

        Normalization Methods:
        - CPM: Counts Per Million
        - RPKM: Reads Per Kilobase Million
        - TPM: Transcripts Per Million

        Parameters
        ----------
        infile : str
            Path to featureCounts output.
        method : str
            'cpm', 'rpkm', 'fpkm', 'tpm' or 'count'.
            count: do nothing
        gtf : str
            Path to GTF annotation file.
        convert_to_gene_name : bool
            Whether to convert gene IDs to gene names.
        remove_version : bool
            Whether to remove version numbers from gene IDs.

        Returns
        -------
        pd.DataFrame
            Normalized matrix or count matrix
        """
        logger.info(f"Starting normalization using method: {method},convert_to_gene_name: {convert_to_gene_name}, remove_version: {remove_version}")
        target_gtf = gtf or self.gtf_path
        df_counts = pd.read_csv(infile, sep="\t", comment='#')
        df_counts.drop(columns=['Chr', 'Start', 'End', 'Strand'], inplace=True, errors='ignore')
        df_counts = df_counts.set_index('Geneid')
        df_counts.columns = [self.extract_sample_name(c) for c in df_counts.columns]
        
        if method == "cpm":
            df = self.compute_cpm(df_counts, length="Length")
        elif method in ["rpkm", "fpkm"]:
            df = self.compute_rpkm(df_counts, length="Length")
        elif method == "tpm":
            df = self.compute_tpm(df_counts, length="Length")
        else:
            df_counts.drop(columns = ["Length"])
            logger.info("don't apply any normalization methods,only drop Length column")
            
        df = df.reset_index()
        if convert_to_gene_name:
            if not target_gtf or os.path.getsize(target_gtf) == 0:
                raise ValueError("GTF file is missing or empty.")
            df = convert_featurecounts_gene_ids(df, target_gtf)
        else:
            if remove_version:
                df['Geneid'] = df['Geneid'].apply(lambda x: x.split('.')[0])
        
        return df

    def combine_PE_SE(self, PE: str, SE: str):
        df_PE = pd.read_csv(PE, sep="\t", comment='#')
        df_PE.columns = [self.extract_sample_name(c) for c in df_PE.columns]
        df_SE = pd.read_csv(SE, sep="\t", comment='#')
        df_SE.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'], inplace=True)
        df_SE.columns = [self.extract_sample_name(c) for c in df_SE.columns]
        return pd.merge(df_PE, df_SE, on="Geneid")

if __name__ == "__main__":
    human_gtf = "/data/pub/zhousha/Reference/human/GENCODE/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf"
    normalizer = RNASeqNormalizer(gtf_path=human_gtf)
    
    tpm = normalizer.run_norm(
        "/data/pub/zhousha/Totipotent20251031/RNAseqML/matrix/huam_all_count.tsv",
        method="tpm",
        convert_to_gene_name=False,
        remove_version=True
    )
    tpm.to_csv("/data/pub/zhousha/Totipotent20251031/RNAseqML/matrix/human_all_tpm.tsv", sep="\t", index=False)
    df = pd.read_csv("/data/pub/zhousha/Totipotent20251031/RNAseqML/matrix/human_all_tpm.tsv", sep="\t", index_col=0)
