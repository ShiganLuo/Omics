import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, FactorAnalysis, NMF

# ------------------------ 工具函数 ------------------------
def minmax_scale_df(df):
    """对每一列做 0-1 缩放；若常数列则保持 0"""
    return df.apply(lambda x: (x - x.min()) / (x.max() - x.min() + 1e-8), axis=0)


def latent_factor_modeling_with_alignment(
        df_resid: pd.DataFrame,
        n_components: int = 5,
        scale_gene: bool = False,
        scale_output: str = None,
        nmf_random_state: int = 42):
    """
    同时运行 NMF / PCA / FA，并将 PCA/FA 方向对齐到 NMF。
    返回结果包含 factors / loadings / extra
    """
    X = df_resid.T.copy()  # 行=样本, 列=基因

    if scale_gene:
        X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-8)

    results = {}

    # ------------------------ 1) NMF ------------------------
    X_pos = X - X.min().min() + 1e-8
    nmf_model = NMF(n_components=n_components, init="nndsvda",
                    random_state=nmf_random_state, max_iter=2000)
    nmf_factors_arr = nmf_model.fit_transform(X_pos)
    nmf_loadings_arr = nmf_model.components_.T

    nmf_factors = pd.DataFrame(nmf_factors_arr, index=X.index,
                               columns=[f"LV{i+1}" for i in range(n_components)])
    nmf_loadings = pd.DataFrame(nmf_loadings_arr, index=X.columns,
                                columns=[f"LV{i+1}" for i in range(n_components)])

    results["nmf"] = {"factors": nmf_factors.copy(),
                      "loadings": nmf_loadings.copy(),
                      "extra": {"model": nmf_model}}

    baseline_lv1 = nmf_factors["LV1"]

    # ------------------------ 2) PCA ------------------------
    pca_model = PCA(n_components=n_components, random_state=nmf_random_state)
    pca_factors_arr = pca_model.fit_transform(X)
    pca_loadings_arr = pca_model.components_.T

    pca_factors = pd.DataFrame(pca_factors_arr, index=X.index,
                               columns=[f"LV{i+1}" for i in range(n_components)])
    pca_loadings = pd.DataFrame(pca_loadings_arr, index=X.columns,
                                columns=[f"LV{i+1}" for i in range(n_components)])

    # 对齐方向
    common = pca_factors.index.intersection(baseline_lv1.index)
    corr = np.corrcoef(pca_factors.loc[common, "LV1"],
                       baseline_lv1.loc[common])[0, 1]
    if np.isfinite(corr) and corr < 0:
        pca_factors *= -1
        pca_loadings *= -1

    results["pca"] = {"factors": pca_factors.copy(),
                      "loadings": pca_loadings.copy(),
                      "extra": {"explained_variance_ratio": pca_model.explained_variance_ratio_}}

    # ------------------------ 3) Factor Analysis ------------------------
    fa_model = FactorAnalysis(n_components=n_components)
    fa_factors_arr = fa_model.fit_transform(X)
    fa_loadings_arr = fa_model.components_.T

    fa_factors = pd.DataFrame(fa_factors_arr, index=X.index,
                              columns=[f"LV{i+1}" for i in range(n_components)])
    fa_loadings = pd.DataFrame(fa_loadings_arr, index=X.columns,
                               columns=[f"LV{i+1}" for i in range(n_components)])

    corr = np.corrcoef(fa_factors.loc[common, "LV1"],
                       baseline_lv1.loc[common])[0, 1]
    if np.isfinite(corr) and corr < 0:
        fa_factors *= -1
        fa_loadings *= -1

    results["fa"] = {"factors": fa_factors.copy(),
                     "loadings": fa_loadings.copy(),
                     "extra": {}}

    # ------------------------ 4) 可选 0-1 缩放 ------------------------
    if scale_output is not None:
        for method in ["nmf", "pca", "fa"]:
            if scale_output in ["factors", "both"]:
                results[method]["factors"] = minmax_scale_df(results[method]["factors"])
            if scale_output in ["loadings", "both"]:
                results[method]["loadings"] = minmax_scale_df(results[method]["loadings"])

    return results



# ------------------------ 5) 保存结果 ------------------------
def save_latent_results(results, outdir):
    os.makedirs(outdir, exist_ok=True)
    for method, res in results.items():
        res["factors"].to_csv(os.path.join(outdir, f"{method}_factors.csv"))
        res["loadings"].to_csv(os.path.join(outdir, f"{method}_loadings.csv"))
        extra = res["extra"]
        if extra:
            # 以单行保存
            pd.DataFrame([extra]).to_csv(os.path.join(outdir, f"{method}_extra.csv"), index=False)


if __name__ == "__main__":
    # ------------------------ 6) 调用示例 ------------------------
    df_resid = pd.read_csv("/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229_514_resid.tsv",sep="\t",index_col=0)
    results = latent_factor_modeling_with_alignment(df_resid, n_components=5,
                                                    scale_gene=True, scale_output="both")
    save_latent_results(results, "/home/luosg/Data/genomeStability/output/result/latent")
