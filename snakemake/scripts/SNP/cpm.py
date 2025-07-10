import pandas as pd
import argparse
import pandas as pd
import numpy as np
from scipy import stats

def calc_tmm_normalization_factors(df):
    # 1. 计算基因之间的 M 值 (log2(count_i / count_j))，避免除以零
    epsilon = 1e-10  # 设置一个小的常数来防止除零错误
    m_values = np.log2(np.clip(df.values[:, :, None] / (df.values[:, None, :] + epsilon), epsilon, None))

    # 2. 修剪 M 值，去掉前后 5% 的极端值
    # 通过对每一对基因之间的比较进行修剪，避免外部极端值影响
    trimmed_m_values = np.percentile(m_values, 95, axis=2) - np.percentile(m_values, 5, axis=2)
    
    # 3. 计算每个样本的标准化因子，使用中位数作为标准化因子
    normalization_factors = np.median(trimmed_m_values, axis=0)
    
    return normalization_factors

def compute_cpm(df: pd.DataFrame, outfile: str, normalization_factors=None, prior_count=2):
    # 计算总计数（库大小），如果没有标准化因子，直接使用列的总和作为库大小
    library_sizes = df.sum(axis=0) # 加入 prior_count 来避免零值

    # 如果提供了标准化因子，将其应用于库大小
    if normalization_factors is not None:
        library_sizes /= normalization_factors

    # 计算 CPM：将每个基因的计数值除以库大小并乘以 1e6
    cpm = df.div(library_sizes, axis=1) * 1e6
    
    # 将计算结果保存为 CSV 文件
    cpm.to_csv(outfile, index=True, header=True)
    
    return cpm




if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description="caculate cpm for TEtranscipts,support multi sample")
    # parser.add_argument('--input', type=str, required=True, help='the gene counts matrix')
    # parser.add_argument('--output', type=str, required=True, help='the cpm outfile')
    # args = parser.parse_args()
    # infile = args.input
    # outfile = args.output
    infile = "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/counts/mouseTEcount.cntTable"
    outfile = "a.csv"
    df = pd.read_csv(infile,sep="\t",index_col=0)
    # 计算TMM标准化因子
    normalization_factors = calc_tmm_normalization_factors(df)

    # 计算CPM
    cpm = compute_cpm(df, outfile)

    # 输出标准化后的CPM
    print(cpm)
