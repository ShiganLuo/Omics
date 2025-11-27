import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import zscore

def read_gmt(gmt_path):
    """简单读取 GMT，返回 dict: name -> list(genes)"""
    g = {}
    with open(gmt_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                name = parts[0]
                genes = [x.upper() for x in parts[2:] if x]
                g[name] = genes
    return g

def filter_genes_for_matrix(df: pd.DataFrame, genes: list, min_std: float = 1e-6):
    """返回在 df 中存在且方差大于 min_std 的基因列表"""
    present = [g for g in genes if g in df.index]
    good = [g for g in present if df.loc[g].std(ddof=0) > min_std]
    return good

# 方法 A：Mean z-score（基因先做 z-score，再求均值）
def compute_gsi_mean_zscore(df_resid: pd.DataFrame, gsi_genes: list, min_std: float = 1e-6):
    """
    df_resid: genes x samples (DataFrame)
    gsi_genes: list of gene names (matching df index)
    returns: pd.Series of GSI per sample
    """
    genes = filter_genes_for_matrix(df_resid, gsi_genes, min_std=min_std)
    if len(genes) == 0:
        raise ValueError("No GSI genes found with sufficient variance.")
    # compute z per gene across samples (axis=1 since rows=genes)
    # using ddof=0 behavior via scipy/zscore default: here use pandas:
    z = df_resid.loc[genes].apply(lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-9), axis=1)
    gsi = z.mean(axis=0)  # mean over genes -> series indexed by samples
    return gsi

# 方法 B：PCA PC1（对基因做可选标准化）
def compute_gsi_pca_pc1(df_resid: pd.DataFrame, gsi_genes: list, min_std: float = 1e-6, scale_genes: bool = True):
    """
    returns: pd.Series of PC1 per sample (signed so it correlates positively with mean zscore)
    """
    genes = filter_genes_for_matrix(df_resid, gsi_genes, min_std=min_std)
    if len(genes) == 0:
        raise ValueError("No GSI genes found with sufficient variance.")
    mat = df_resid.loc[genes].T  # samples x genes for sklearn
    if scale_genes:
        # Standardize columns (genes)
        mat = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-9)
    pca = PCA(n_components=1, random_state=0)
    pc1 = pca.fit_transform(mat)[:, 0]
    pc1_ser = pd.Series(pc1, index=mat.index, name="PC1")
    return pc1_ser, pca.explained_variance_ratio_[0]

# 可选方法 C：ssGSEA via gseapy（如果安装了 gseapy）
def compute_gsi_ssgsea(df_log2: pd.DataFrame, gene_set: list, sample_axis='columns'):
    """
    使用 gseapy.ssgsea（需要 gseapy 安装）
    df_log2: genes x samples (index genes)
    gene_set: list of genes (names matching index)
    returns: pd.Series of ssgsea scores (samples)
    """
    try:
        import gseapy as gp
    except Exception as e:
        raise ImportError("gseapy not installed. pip install gseapy") from e

    # prepare gmt-like dict for gseapy
    # gseapy expects samples x genes in DataFrame, and gene names as index
    # gseapy.ssgsea requires a dict of gene sets or gmt file.
    gs_dict = {"GSI": gene_set}
    # gseapy.ssgsea wants DataFrame: samples x genes
    data = df_log2.T  # samples x genes
    # run
    res = gp.ssgsea(data=data, gene_sets=gs_dict, sample_norm_method='rank', outdir=None, verbose=False)
    # result is DataFrame: gene_set x sample
    scores = res.res2d.loc["GSI"]
    return pd.Series(scores, index=data.index)


# ------------- utility: build combined GSI and compare -------------
def build_and_compare_gsi(df_resid: pd.DataFrame, gsi_genes: list, min_std: float = 1e-6, scale_genes=True):
    """
    compute mean_zscore and pca_pc1, align sign and return DataFrame with both scores and correlation
    """
    gsi_mean = compute_gsi_mean_zscore(df_resid, gsi_genes, min_std=min_std)
    pc1, var_explained = compute_gsi_pca_pc1(df_resid, gsi_genes, min_std=min_std, scale_genes=scale_genes)
    # align sign so that PC1 correlates positively with mean_zscore
    corr = pc1.corr(gsi_mean)
    if corr < 0:
        pc1 = -pc1
        corr = -corr
    res = pd.DataFrame({"GSI_mean_z": gsi_mean, "GSI_PC1": pc1})
    return res, corr, var_explained

if __name__ == "__main__":
    # 读取残差矩阵
    df_resid = pd.read_csv("residual.tsv", index_col=0, sep="\t")
    # 读取 GSI gene list（示例：从 gmt）
    gsi_dict = read_gmt("/home/luosg/Data/genomeStability/data/final.gmt")
    # 假设 final.gmt 中有一个条目 'GSI' 或你手动合并感兴趣的基因
    gsi_genes = list({g for genes in gsi_dict.values() for g in genes})  # 所有不重复基因，或换成某个 pathway

    # 计算并比较
    res_df, corr_val, var_exp = build_and_compare_gsi(df_resid, gsi_genes, min_std=1e-6, scale_genes=True)
    print("PC1 vs mean_zscore corr:", corr_val)
    print("PC1 variance explained:", var_exp)
    res_df.to_csv("gsi_scores_comparison.tsv", sep="\t")