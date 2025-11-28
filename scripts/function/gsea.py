import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Optional
def run_gsva(
        count_matrix:str,
        gmt_pathway:str,
        outdir:str
):
    es = gp.gsva(
    data = count_matrix,
    gene_sets = gmt_pathway,
    outdir = outdir
    )
    return es

def run_gsea(rnk_series:pd.Series,gmt_pathway:str,outdir:str):
    gsea_results = gp.prerank(
        rnk=rnk_series,               # 排名 Series
        gene_sets=gmt_pathway,          # 基因集或GMT文件路径
        outdir=outdir,    # 结果输出目录
        min_size=5,
        max_size=500,
        permutation_num=1000,         # 排列次数
        seed=42,                      # 随机种子
        format='png',                 # 输出图表的格式
        no_plot=False,                # 生成富集图
        threads=4                   # 使用多核加速
    )
    return gsea_results



def plot_gsea_from_csv(csv_path: str, ranked_genes: List[str], out_png: str, top_n: Optional[int] = None):
    """
    从 GSEA CSV 文件绘制累积 ES 曲线，并保存图片。
    默认绘制所有通路，可选绘制前 top_n 条显著通路。

    Parameters:
    -----------
    csv_path : str
        GSEA prerank 结果 CSV 文件路径
    ranked_genes : List[str]
        所有基因按排名顺序的列表
    out_png : str
        保存输出图片路径
    top_n : Optional[int]
        绘制前多少条显著通路，默认 None 表示绘制所有
    """
    # 读取 CSV
    df = pd.read_csv(csv_path)

    # 找到 FDR 列
    fdr_col_candidates = [c for c in df.columns if 'fdr' in c.lower()]
    if not fdr_col_candidates:
        raise ValueError("CSV 中未找到 FDR 列")
    fdr_col = fdr_col_candidates[0]

    # 如果指定 top_n，则按 FDR 排序取前 top_n
    if top_n is not None:
        df = df.sort_values(fdr_col).head(top_n)

    plt.figure(figsize=(14, 6))

    N = len(ranked_genes)
    for _, row in df.iterrows():
        term = row['Term']
        nes = row['NES']

        # 命中基因
        lead_genes = str(row['Lead_genes']).split(";")

        # 构建 hit 向量
        hits = np.array([1 if g in lead_genes else 0 for g in ranked_genes])
        nh = hits.sum()
        no = N - nh
        if nh == 0:
            continue

        # 计算累积 ES
        running_es = np.cumsum(hits / nh - (1 - hits) / no)

        # 绘制
        plt.plot(running_es, label=f"{term} (NES={nes:.2f})")

    plt.xlabel("Ranked Genes")
    plt.ylabel("Enrichment Score (ES)")
    plt.title("GSEA Pathways")
    plt.legend(fontsize=6, ncol=2)  # legend 调小字体，分两列
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"图已保存到 {out_png}")

if __name__ == "__main__":
    # run_gsva(
    #     count_matrix= "/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229.tsv",
    #     gmt_pathway= "/home/luosg/Data/genomeStability/data/final.gmt",
    #     outdir= "/home/luosg/Data/genomeStability/output/result/gsva230"
    # )
    ################ gsea preprank
    # df = pd.read_csv("/home/luosg/Data/genomeStability/output/result/latent/fa_loadings.csv",index_col=0)
    # geneRank = df["LV1"]
    # gsea_results = run_gsea(geneRank,"/home/luosg/Data/genomeStability/data/final.gmt","/home/luosg/Data/genomeStability/output/result/gsea")
    ########## multiple curve plot
    # df_gsea = pd.read_csv("/home/luosg/Data/genomeStability/output/result/gsea/gseapy.gene_set.prerank.report.csv")
    df_rnk = pd.read_csv("/home/luosg/Data/genomeStability/output/result/gsea/prerank_data.rnk",header=None,sep="\t")
    geneRnk = df_rnk[0].to_list()
    plot_gsea_from_csv("/home/luosg/Data/genomeStability/output/result/gsea/gseapy.gene_set.prerank.report.csv",
                       geneRnk,"/home/luosg/Data/genomeStability/output/result/gsea/all.png",5)

    
