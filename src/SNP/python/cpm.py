import pandas as pd
import argparse
import pandas as pd
import numpy as np
from scipy import stats

def calc_tmm_normalization_factors(df):
    epsilon = 1e-10
    m_values = np.log2(np.clip(df.values[:, :, None] / (df.values[:, None, :] + epsilon), epsilon, None))

    trimmed_m_values = np.percentile(m_values, 95, axis=2) - np.percentile(m_values, 5, axis=2)
    
    normalization_factors = np.median(trimmed_m_values, axis=0)
    
    return normalization_factors

def compute_cpm(df: pd.DataFrame, outfile: str, normalization_factors=None, prior_count=2):
    library_sizes = df.sum(axis=0)

    if normalization_factors is not None:
        library_sizes /= normalization_factors

    cpm = df.div(library_sizes, axis=1) * 1e6
    
    cpm.to_csv(outfile, index=True, header=True)
    
    return cpm




if __name__ == '__main__':
    infile = "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/counts/mouseTEcount.cntTable"
    outfile = "a.csv"
    df = pd.read_csv(infile,sep="\t",index_col=0)
    normalization_factors = calc_tmm_normalization_factors(df)

    cpm = compute_cpm(df, outfile)

    print(cpm)
