import pandas as pd
def compute_cpm(counts_df):
    counts_only = counts_df.drop(columns=['length'])
    cpm = counts_only.div(counts_only.sum(axis=0), axis=1) * 1e6
    return cpm

def compute_rpkm(counts_df):
    counts_only = counts_df.drop(columns=['length'])
    gene_length_kb = counts_df['length'] / 1000  # bp -> kb
    rpkm = counts_only.div(gene_length_kb, axis=0)
    rpkm = rpkm.div(rpkm.sum(axis=0) / 1e6, axis=1)
    return rpkm

def compute_tpm(counts_df):
    counts_only = counts_df.drop(columns=['length'])
    gene_length_kb = counts_df['length'] / 1000
    rpk = counts_only.div(gene_length_kb, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    return tpm

if __name__ == "__main__":
    df_counts = pd.DataFrame({
    'gene_id': ['gene1','gene2','gene3'],
    'length': [1000, 2000, 500],  # 基因长度，单位 bp
    'sample1': [100, 150, 80],
    'sample2': [120, 130, 90]
    })

    df_counts = df_counts.set_index('gene_id')