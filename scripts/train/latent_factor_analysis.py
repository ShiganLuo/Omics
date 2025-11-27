#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
端到端潜在变量建模（PCA / FA / NMF）

输入：
    df_resid: 已校正后的残差矩阵（行=基因，列=样本）

输出：
    sample_factors.csv
    gene_loadings.csv
    pca_scree_plot.png（仅 PCA 时）
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, FactorAnalysis, NMF


def latent_factor_modeling(df_resid,
                           method="pca",
                           n_components=5,
                           scale=False):
    """
    对残差矩阵执行潜在变量建模 (PCA / Factor Analysis / NMF)
    """

    X = df_resid.T.copy()

    if scale:
        X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-8)

    extra = {}

    if method.lower() == "pca":
        model = PCA(n_components=n_components)
        sample_factors = model.fit_transform(X)
        gene_loadings = model.components_.T
        extra["explained_variance_ratio"] = model.explained_variance_ratio_

    elif method.lower() == "fa":
        model = FactorAnalysis(n_components=n_components)
        sample_factors = model.fit_transform(X)
        gene_loadings = model.components_.T

    elif method.lower() == "nmf":
        X_pos = X - X.min().min() + 1e-8
        model = NMF(n_components=n_components, init="nndsvda", random_state=42)
        sample_factors = model.fit_transform(X_pos)
        gene_loadings = model.components_.T

    else:
        raise ValueError("method 必须为 pca / fa / nmf")

    sample_factors = pd.DataFrame(
        sample_factors,
        index=X.index,
        columns=[f"LV{i+1}" for i in range(n_components)]
    )

    gene_loadings = pd.DataFrame(
        gene_loadings,
        index=df_resid.index,
        columns=[f"LV{i+1}" for i in range(n_components)]
    )

    return sample_factors, gene_loadings, extra


def plot_pca_scree(explained_ratio, out="pca_scree_plot.png"):
    plt.figure(figsize=(6, 4))
    plt.plot(range(1, len(explained_ratio) + 1), explained_ratio, marker="o")
    plt.xlabel("PC")
    plt.ylabel("Explained variance ratio")
    plt.title("PCA Scree Plot")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resid", required=True,
                        help="输入残差矩阵 TSV（行=基因，列=样本）")
    parser.add_argument("--method", default="pca",
                        help="pca / fa / nmf")
    parser.add_argument("--k", type=int, default=5,
                        help="提取潜在因子数")
    parser.add_argument("--scale", action="store_true",
                        help="是否对基因 Z-score")
    parser.add_argument("--prefix", default="latent",
                        help="输出前缀")
    args = parser.parse_args()

    # ---------- 读取矩阵 ----------
    df_resid = pd.read_csv(args.resid, sep="\t", index_col=0)

    # ---------- 模型 ----------
    sample_factors, gene_loadings, extra = latent_factor_modeling(
        df_resid=df_resid,
        method=args.method,
        n_components=args.k,
        scale=args.scale
    )

    # ---------- 输出 ----------
    sample_factors.to_csv(f"{args.prefix}_sample_factors.tsv", sep="\t")
    gene_loadings.to_csv(f"{args.prefix}_gene_loadings.tsv", sep="\t")

    if args.method == "pca":
        plot_pca_scree(extra["explained_variance_ratio"],
                       out=f"{args.prefix}_pca_scree.png")

    print("✔ 潜在因子分析完成")
    print(f"  Sample factors → {args.prefix}_sample_factors.tsv")
    print(f"  Gene loadings → {args.prefix}_gene_loadings.tsv")
    if args.method == "pca":
        print(f"  Scree plot → {args.prefix}_pca_scree.png")


if __name__ == "__main__":
    main()
