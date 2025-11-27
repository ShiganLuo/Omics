import pandas as pd
from pycombat import Combat

def batch_correct_pycombat(df_log2: pd.DataFrame,
                           df_sample: pd.DataFrame,
                           batch_col: str,
                           covariate_cols: list = None,
                           sample_id_col: str = "Sample_id"):
    """
    使用 pyComBat 对 log2(TPM+1) 表达矩阵做批次 (batch) 校正。

    Parameters
    ----------
    df_log2 : pd.DataFrame
        行 = 基因，列 = 样本 (样本 ID 对应 df_sample[sample_id_col])
    df_sample : pd.DataFrame
        包含样本批次信息和协变量信息
    batch_col : str
        主批次列名 (df_sample 中)
    covariate_cols : list of str, optional
        要作为协变量 (covariates) 的列名列表 (例如：["Cell_line", "condition"])
    sample_id_col : str
        样本 ID 列名，在 df_sample 中，对应 df_log2 的列名

    Returns
    -------
    df_corrected : pd.DataFrame
        批次校正后的表达矩阵 (same shape as df_log2)
    """

    # —— 1. 对齐样本 —— #
    df_sample_unique = df_sample.drop_duplicates(subset=[sample_id_col])
    df_sample_sel = df_sample_unique[df_sample_unique[sample_id_col].isin(df_log2.columns)]

    # reorder df_log2 columns to match df_sample order
    sel_ids = list(df_sample_sel[sample_id_col])
    df_log2_sel = df_log2.loc[:, sel_ids]

    # —— 2. 构建 batch labels 和 covariate dataframe —— #
    batch_labels = df_sample_sel[batch_col].values

    covariates = None
    if covariate_cols:
        cov_dict = {}
        for col in covariate_cols:
            cov_dict[col] = df_sample_sel[col].astype(str).values
        covariates = pd.DataFrame(cov_dict, index=sel_ids)

    # —— 3. 调用 pyComBat —— #
    # pyComBat expects data as samples × features (i.e. transpose)
    Y = df_log2_sel.T  # now rows = samples, cols = genes

    combat = Combat()
    if covariates is not None:
        Y_corr = combat.fit_transform(Y=Y, b=batch_labels, C=covariates)
    else:
        Y_corr = combat.fit_transform(Y=Y, b=batch_labels)

    # —— 4. 转回 genes × samples 格式 —— #
    df_corr = pd.DataFrame(Y_corr.T, index=df_log2_sel.index, columns=df_log2_sel.columns)
    return df_corr




if __name__ == "__main__":
    df_sample = pd.read_csv("/home/luosg/Data/genomeStability/data/target_fq.tsv",sep="\t")
    df_log2 = pd.read_csv("/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229.tsv",index_col=0,sep="\t")
    df_corrected = batch_correct_pycombat(
        df_log2 = df_log2,
        df_sample = df_sample,
        batch_col = "Study_name",
        covariate_cols = ["Cell_line", "condition"],
        sample_id_col = "Sample_id"
    )
    df_corrected.to_csv("/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229_debatch.tsv",sep="\t")

