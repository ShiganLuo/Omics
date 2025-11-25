from pdb import run
import pandas as pd
import re
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

def extract_sample_name(col):
    m = re.search(r'output/2pass/([^/]+)/', col)
    if m:
        return m.group(1)
    else:
        return col 



def run_norm(infile:str,method:str):
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
    return df

def combine_PE_SE(PE:str,SE:str):
    df_PE = pd.read_csv(PE,sep="\t",comment='#')
    df_PE.columns = [extract_sample_name(c) for c in df_PE.columns]
    df_SE = pd.read_csv(SE,sep="\t",comment='#')
    df_SE.drop(columns=['Chr', 'Start', 'End', 'Strand','Length'],inplace=True)
    df_SE.columns = [extract_sample_name(c) for c in df_SE.columns]
    df_new = pd.merge(df_PE,df_SE,on="Geneid")
    return df_new




if __name__ == "__main__":
    # tpm1 = run_norm("/home/luosg/Data/genomeStability/output/count/featureCounts/human_paired_count.tsv","tpm")
    # tpm1.to_csv("/home/luosg/Data/genomeStability/output/count/featureCounts/human_paired_tpm.tsv",sep="\t")
    # tpm2 = run_norm("/home/luosg/Data/genomeStability/output/count/featureCounts/human_single_count.tsv","tpm")
    # tpm2.to_csv("/home/luosg/Data/genomeStability/output/count/featureCounts/human_single_tpm.tsv",sep="\t")  
    # df = combine_PE_SE("/home/luosg/Data/genomeStability/output/counts/featureCounts/human/human_paired_count.tsv",
    #                    "/home/luosg/Data/genomeStability/output/counts/featureCounts/human/human_single_count.tsv")
    # df.to_csv("/home/luosg/Data/genomeStability/output/counts/featureCounts/human/huam_all_count.tsv",sep="\t",index=False)
    df_tpm = run_norm("/home/luosg/Data/genomeStability/output/counts/featureCounts/human/huam_all_count.tsv","tpm")
    df_tpm.to_csv("/home/luosg/Data/genomeStability/output/counts/featureCounts/human/human_all_tpm.tsv",sep="\t")



    