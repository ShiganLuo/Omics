import pandas as pd
import re
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
scripts_dir = current_path.parents[1]
sys.path.append(str(scripts_dir))
from annotation.gene_id2name import convert_featurecounts_gene_ids

def compute_cpm(counts_df:pd.DataFrame,length:str):
    counts_only = counts_df.drop(columns=[length])
    cpm = counts_only.div(counts_only.sum(axis=0), axis=1) * 1e6
    return cpm

def compute_rpkm(counts_df:pd.DataFrame,length:str):
    """
    rpkm or fpkm, SE or PE
    """
    counts_only = counts_df.drop(columns=[length])
    gene_length_kb = counts_df[length] / 1000  # bp -> kb
    rpkm = counts_only.div(gene_length_kb, axis=0)
    rpkm = rpkm.div(rpkm.sum(axis=0) / 1e6, axis=1)
    return rpkm

def compute_tpm(counts_df:pd.DataFrame,length:str):
    counts_only = counts_df.drop(columns=[length])
    gene_length_kb = counts_df[length] / 1000
    rpk = counts_only.div(gene_length_kb, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    return tpm

def extract_sample_name(
    col: str,
    pattern_aligned: str = r'([^/]+?)Aligned',
    pattern_fallback: str = r'/([^/]+)\.[^.]+$'
):
    """
    Extract sample name from featureCounts / BAM column names.

    Rules:
    1. If the column name contains 'Aligned', extract the part before 'Aligned'
       e.g. GSM5837959Aligned.sortedByCoord.out -> GSM5837959
    2. Otherwise, try to extract the filename (without extension) from a file path
       e.g. /path/to/SRR123456.fastq.gz -> SRR123456
    3. If no pattern matches, return the original column name
    """

    # 1. Prefer extracting the part before 'Aligned'
    m = re.search(pattern_aligned, col)
    if m:
        return m.group(1)

    # 2. Fallback: extract filename without extension from path
    m = re.search(pattern_fallback, col)
    if m:
        return m.group(1)

    # 3. No match: return original value
    return col


def run_norm(
    infile: str,
    gtf: str,
    method: str = "tpm"
) -> pd.DataFrame:
    """
    Perform gene expression normalization on featureCounts output.

    This function reads a featureCounts count matrix, preprocesses it by
    removing genomic coordinate columns, normalizes raw counts using the
    specified method (CPM / RPKM / FPKM / TPM), and converts gene IDs to
    gene names based on the provided GTF annotation.

    Parameters
    ----------
    infile : str
        Path to the featureCounts output file (tab-separated).
        The file must contain at least the following columns:
        'Geneid', 'Length', and sample count columns.

    method : str
        Normalization method to apply.
        Supported methods:
        - 'cpm'   : Counts Per Million
        - 'rpkm'  : Reads Per Kilobase Million
        - 'fpkm'  : Fragments Per Kilobase Million (treated the same as RPKM)
        - 'tpm'   : Transcripts Per Million

    gtf : str
        Path to the GTF annotation file used to map gene IDs
        (e.g. Ensembl IDs) to gene names.

    Returns
    -------
    pd.DataFrame
        A normalized gene expression matrix with gene names as rows
        and samples as columns.

    Raises
    ------
    ValueError
        If an unsupported normalization method is provided.
    """
    df_counts = pd.read_csv(infile,
                            sep="\t",
                            comment='#')
    df_counts.drop(columns=['Chr', 'Start', 'End', 'Strand'],inplace=True)
    df_counts = df_counts.set_index('Geneid')
    df_counts.columns = [extract_sample_name(c) for c in df_counts.columns]
    df = pd.DataFrame()
    if method == "cpm":
        df = compute_cpm(df_counts,length="Length")
    elif method == "rpkm" or method == "fpkm":
        df = compute_rpkm(df_counts,length="Length")
    elif method == "tpm":
        df = compute_tpm(df_counts,length="Length")
    else:
        raise ValueError("please input correct normalization method")
    df = df.reset_index()
    df_new = convert_featurecounts_gene_ids(df,gtf)
    return df_new

def combine_PE_SE(PE:str,SE:str):
    df_PE = pd.read_csv(PE,sep="\t",comment='#')
    df_PE.columns = [extract_sample_name(c) for c in df_PE.columns]
    df_SE = pd.read_csv(SE,sep="\t",comment='#')
    df_SE.drop(columns=['Chr', 'Start', 'End', 'Strand','Length'],inplace=True)
    df_SE.columns = [extract_sample_name(c) for c in df_SE.columns]
    df_new = pd.merge(df_PE,df_SE,on="Geneid")
    return df_new




if __name__ == "__main__":
    human_gtf = "/disk5/luosg/Reference/GENCODE/human/GRCh38/gencode.v49.primary_assembly.basic.annotation.gtf"
    mouse_gtf = "/disk5/luosg/Reference/GENCODE/mouse/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf"
    # tpm1 = run_norm("/home/luosg/Data/genomeStability/output/count/featureCounts/human_paired_count.tsv","tpm")
    # tpm1.to_csv("/home/luosg/Data/genomeStability/output/count/featureCounts/human_paired_tpm.tsv",sep="\t")
    # tpm2 = run_norm("/home/luosg/Data/genomeStability/output/count/featureCounts/human_single_count.tsv","tpm")
    # tpm2.to_csv("/home/luosg/Data/genomeStability/output/count/featureCounts/human_single_tpm.tsv",sep="\t")  
    # df = combine_PE_SE("/home/luosg/Data/genomeStability/output/counts/featureCounts/human/human_paired_count.tsv",
    #                    "/home/luosg/Data/genomeStability/output/counts/featureCounts/human/human_single_count.tsv")
    # df.to_csv("/home/luosg/Data/genomeStability/output/counts/featureCounts/human/huam_all_count.tsv",sep="\t",index=False)
    df_tpm1 = run_norm("/disk5/luosg/GCN2_20251224/output/counts/featureCounts/mouse/GCN2pub_paired_count.tsv",mouse_gtf,"tpm")
    df_tpm2 = run_norm("/disk5/luosg/GCN2_20251224/output/counts/featureCounts/mouse/GCN2seq_paired_count.tsv",mouse_gtf,"tpm")
    print(df_tpm1.shape)
    print(df_tpm2.shape)
    df_new = pd.merge(df_tpm1,df_tpm2,on="gene_name",how="inner")
    print(df_new.shape)
    df_new.to_csv("/disk5/luosg/GCN2_20251224/output/counts/featureCounts/mouse/GCN2_paired_count.tpm",sep="\t",index=False)




    