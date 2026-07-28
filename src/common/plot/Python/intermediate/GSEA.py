import anndata as ad
import scanpy as sc
import numpy as np
import logging
import sys
import pandas as pd
import decoupler
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap

logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # 指定输出到 stdout 而不是 stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)
def GSEA_rankGene(
        adata:ad.AnnData,
        control:str,
        gsea_outfile:str,
        reactome_file:str = "/disk5/luosg/scRNAseq/data/reactiome.csv",
        reference:str = "rest",
        groupby:str="celltype",
        method:str = "t-test",
        use_raw:bool = False,
        decoupler_times = 10000000

):
    logging.info(f"parements: {control} vs {reference}, groupby: {groupby},method: {method}, use_raw: {use_raw}")
    sc.tl.rank_genes_groups(adata, 
                            groupby = groupby, 
                            groups = control, 
                            method=method,
                            reference = reference,
                            key_added=method,
                            use_raw=use_raw)
    celltype_condition = control
    # extract scores
    t_stats = (
        # Get dataframe of DE results for condition vs. rest
        sc.get.rank_genes_groups_df(adata, celltype_condition, key="t-test")
        # Subset to highly variable genes
        .set_index("names")
        # .loc[adata.var["highly_variable"]]
        # Sort by absolute score
        .sort_values("scores", key=np.abs, ascending=False)[
            # Format for decoupler
            ["scores"]
        ]
        .rename_axis([celltype_condition], axis=1)
    )
    # sample:stim cell population:FCGR3A+ Monocytes vs rest all cell
    logging.info("Convert gene names from lowercase to uppercase for gsea")
    t_stats.index = t_stats.index.str.upper()

    logging.info(f"the deg shape:{t_stats.shape}")
    reactome = pd.read_csv(reactome_file,sep=",")
    geneset_size = reactome.groupby("geneset").size()
    logging.info(f"filter pathway from {reactome_file}, standard: gene number is in (15,500)")
    gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 500)]
    logging.info(f"rename reactome dataframe: genesymbol:target,geneset:source for decoupler")
    reactome.rename(columns={"genesymbol":"target","geneset":"source"},inplace=True)

    scores, pvals = decoupler.mt.gsea(t_stats.T,reactome[reactome["source"].isin(gsea_genesets)],times=decoupler_times)
    gsea_result = pd.concat({"score": scores.T, "pval": pvals.T}, axis=1).droplevel(level=1,axis=1).sort_values("pval")
    gsea_result.reset_index(names = "source",inplace=True)
    gsea_result.to_csv(gsea_outfile,sep="\t",index=False)
    return gsea_result


def GSEAPlot(
        gsea_result: pd.DataFrame,
        outplot: str,
        gsea_obvious_outfile: str | None = None,
        pval_cutoff: float = 0.05,
        top_flag:bool = False,
        top_number:int = 5
):
    # 过滤
    df = gsea_result[gsea_result["pval"] < pval_cutoff].copy()
    if gsea_obvious_outfile is not None:
        df.to_csv(gsea_obvious_outfile, sep="\t", index=False)

    df["pval_plot"] = -np.log10(df["pval"].replace(0, 1e-30))

    df = df.rename(columns={"score": "ES"}).sort_values("ES", ascending=True)
    if top_flag:
        df = pd.concat([df.head(top_number), df.tail(top_number)])
    # ✅ 自动 wrap 超长 pathway 名
    def wrap_label(x, width=40):
        return "\n".join(textwrap.wrap(x, width=width))

    df["source_wrapped"] = df["source"].apply(lambda x: wrap_label(x, width=35))

    # ✅ 根据 pathway 数量和名称长度动态调整画布大小
    n_pathways = df.shape[0]
    max_label_len = df["source"].str.len().max()

    height = max(4, n_pathways * 0.5)             # 行数影响高度
    width = 10 + max(0, (max_label_len - 35) / 3) # 过长名称自动加宽

    plt.figure(figsize=(width, height))

    # 绘图
    sns.scatterplot(
        data=df,
        x="ES",
        y="source_wrapped",
        size="pval_plot",
        hue="ES",
        sizes=(40, 300),
        palette="coolwarm",
        legend="brief"
    )

    plt.axvline(0, linestyle='--', color='gray')

    plt.legend(
        title=None,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=False,
        borderaxespad=0
    )

    plt.xlabel("Enrichment Score (ES)", fontsize=12)
    plt.ylabel("Pathway", fontsize=12)

    plt.title("GSEA Enrichment Analysis", fontsize=14)

    plt.tight_layout()
    plt.savefig(outplot, dpi=300, bbox_inches="tight")
    plt.close()


def GSEA_BarPlot_Enhanced(
    gsea_result: pd.DataFrame,
    outplot: str,
    pval_cutoff: float = 0.05,
    sort_by: str = "ES",    # "ES" or "pval"
    top_flag: bool = False,
    top_number: int = 5,
    label_wrap: int = 40,
    cmap_name: str = "viridis"
):
    """
    增强版 GSEA barplot：
    - 中心对称：负 ES 左、正 ES 右
    - 按 ES 或 pval 自动排序
    - 颜色按显著性渐变（-log10(pval)）
    - ggplot 风格
    """

    # ------------------------
    # 1. 过滤和预处理
    # ------------------------
    df = gsea_result[gsea_result["pval"] < pval_cutoff].copy()
    df = df.rename(columns={"score": "ES"}).sort_values("ES")
    if top_flag:
        df = pd.concat([df.head(top_number), df.tail(top_number)])
    df["logP"] = -np.log10(df["pval"].replace(0, 1e-7))

    # 名称换行
    df["pathway_wrapped"] = df["source"].apply(
        lambda s: "\n".join(textwrap.wrap(s, label_wrap))
    )


    # ------------------------
    # 4. 图尺寸
    # ------------------------
    height = max(6, len(df) * 0.45)
    fig, ax = plt.subplots(figsize=(12, height))

    # ------------------------
    # 5. 颜色映射
    # ------------------------
    cmap = sns.color_palette(cmap_name, as_cmap=True)
    norm = plt.Normalize(df["logP"].min(), df["logP"].max())
    df["color"] = df["logP"].apply(lambda x: cmap(norm(x)))

    # ------------------------
    # 6. 中心对称 barplot
    # ------------------------
    ax.barh(
        y=df["pathway_wrapped"],
        width=df["ES"],
        color=df["color"],
        edgecolor="black",
        linewidth=0.5
    )

    # ------------------------
    # 7. 中心线（x=0）
    # ------------------------
    ax.axvline(0, color="black", linewidth=1.2)

    # ------------------------
    # 8. 坐标轴与标签
    # ------------------------
    ax.set_xlabel("Enrichment Score (ES)")
    ax.set_ylabel("Pathway")

    ax.set_title(
        "GSEA Enrichment Barplot\n"
        "Left: Negative ES | Right: Positive ES",
        fontsize=14,
        pad=15
    )

    # ------------------------
    # 9. Colorbar (p-value 显著性)
    # ------------------------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("-log10(P-value)")

    # ------------------------
    # 10. 美化
    # ------------------------
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    ax.grid(False, axis="y")

    plt.tight_layout()
    plt.savefig(outplot, dpi=300, bbox_inches="tight")
    plt.close()




def runGSEA(
        adata:ad.AnnData,
        control:str,
        gsea_outfile:str,
        gsea_obviousfile:str,
        gsea_plot:str
):
    gsea_result = GSEA_rankGene(adata,control=control,gsea_outfile=gsea_outfile)
    GSEAPlot(gsea_result,gsea_obvious_outfile=gsea_obviousfile,outplot=gsea_plot)


if __name__ == "__main__":

    # runGSEA()
    gsea_result = pd.read_csv("/disk5/luosg/scRNAseq/output/combine/Intestine/GSEA/table/combined_groupE2_chang_10XSC3_Intestinal_stem_cell-combined_groupCKO_chang_10XSC3_Intestinal_stem_cell_DEG_gene_gsea.tsv",sep="\t")
    # print(gsea_result.head())
    # GSEAPlot(gsea_result,outplot="/disk5/luosg/scRNAseq/output/combine/Intestine/GSEA/plotTop/B_cell_CKO_WT_gsea.png",top_flag=True,top_number=5)


    GSEA_BarPlot_Enhanced(
        gsea_result=gsea_result,
        outplot="/disk5/luosg/scRNAseq/output/combine/Intestine/GSEA/plotTop/Plasmacytoid_dendritic_cells_E2_CKO_gsea.png",
        sort_by="pval",
        top_flag=True,
        top_number=5
    )


