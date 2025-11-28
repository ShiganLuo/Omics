import pandas as pd
import numpy as np

def project_new_samples(
    df_new: pd.DataFrame,
    df_train: pd.DataFrame,
    loadings: pd.DataFrame,
    lv1_ref: pd.Series = None,
    preproc_scale: bool = False
):
    """
    将训练集得到的潜在因子载荷矩阵应用到新样本，得到 LV 得分。

    参数
    ----------
    df_new : pd.DataFrame
        新样本表达矩阵，行=基因，列=样本
    df_train : pd.DataFrame
        训练集表达矩阵（行=基因，列=样本），用于计算均值/标准差或 log2(TPM+1)
    loadings : pd.DataFrame
        训练得到的载荷矩阵，行=基因，列=LV，如 PCA.components_.T
    lv1_ref : pd.Series, optional
        训练集 LV1，用于方向对齐
    preproc_scale : bool
        是否对基因做 z-score（列标准化）

    返回
    -------
    projected_factors : pd.DataFrame
        新样本投影得到的 LV 得分，行=样本，列=LV
    """

    # 对齐基因
    common_genes = df_new.index.intersection(loadings.index)
    if len(common_genes) == 0:
        raise ValueError("新样本和载荷矩阵没有公共基因")
    
    X_new = df_new.loc[common_genes].T.copy()      # 样本 x 基因
    L = loadings.loc[common_genes].values          # 基因 x LV

    # 可选 z-score
    if preproc_scale:
        mean_train = df_train.loc[common_genes].mean(axis=1)
        std_train = df_train.loc[common_genes].std(axis=1) + 1e-8
        X_new = (X_new - mean_train.values) / std_train.values

    # 投影到 LV
    projected = np.dot(X_new.values, L)            # 样本 x LV
    projected_factors = pd.DataFrame(
        projected,
        index=X_new.index,
        columns=loadings.columns
    )

    # LV1 方向对齐
    if lv1_ref is not None:
        common_samples = projected_factors.index.intersection(lv1_ref.index)
        corr = np.corrcoef(projected_factors.loc[common_samples, "LV1"],
                           lv1_ref.loc[common_samples])[0, 1]
        if np.isfinite(corr) and corr < 0:
            projected_factors *= -1

    return projected_factors
