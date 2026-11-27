import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from typing import List, Dict, Literal, Optional

PlotFormat = str

VALID_FORMATS = ("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg",
                 "jpeg", "pgf", "raw", "rgba")


def enrich_go_kegg_analysis(gene_list: List[str], organism: Literal['human', 'mouse'] = 'human', outdir='enrichment_results') -> Dict[str, pd.DataFrame]:
    """
    只做GO和KEGG富集分析，返回结果字典，不绘图。
    """
    os.makedirs(outdir, exist_ok=True)
    results_dict = {}
    go_categories = ['GO_Biological_Process_2021',
                     'GO_Cellular_Component_2021',
                     'GO_Molecular_Function_2021']
    for cat in go_categories:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=cat,
            organism=organism,
            outdir=outdir,
            cutoff=0.05
        )
        df = enr.results
        results_dict[cat] = df
    kegg_cat = 'KEGG_2021_Human' if organism.lower()=='human' else 'KEGG_2019_Mouse'
    enr_kegg = gp.enrichr(
        gene_list=gene_list,
        gene_sets=kegg_cat,
        organism=organism,
        outdir=outdir,
        cutoff=0.05
    )
    df_kegg = enr_kegg.results
    results_dict['KEGG'] = df_kegg
    for key, df in results_dict.items():
        df.to_csv(os.path.join(outdir, f'{key}_enrichment.csv'), index=False)
    print(f"富集分析完成，结果保存在文件夹: {outdir}")
    return results_dict

def plot_enrichment_bar(df: pd.DataFrame, out_base: str, top_n=10, kind='GO',
                        title=None, image_formats: Optional[List[PlotFormat]] = None):
    """
    只负责绘制GO/KEGG富集柱状图,每种请求格式输出一份 "{out_base}.{fmt}"。

    Parameters
    ----------
    df : pd.DataFrame
        单个富集库的结果表(enrichr results)。
    out_base : str
        输出路径前缀(不含扩展名),每种格式写成 ``{out_base}.{fmt}``。
    top_n : int
        按 Adjusted P-value 取前 N 条绘图。
    kind : str
        ``'GO'``/``'KEGG'``,决定配色与轴标签。
    title : str, optional
        图标题;None 时用 ``{kind} Top {top_n} Enrichment``。
    image_formats : List[str], optional
        Save formats. Default ``["png"]``.
    """
    if image_formats is None:
        image_formats = ["png"]
    bad = [f for f in image_formats if f not in VALID_FORMATS]
    if bad:
        raise ValueError(
            f"Unsupported image format(s): {', '.join(bad)}; "
            f"valid: {', '.join(VALID_FORMATS)}")
    if df.empty:
        print(f"{kind} enrichment result is empty, skip plot.")
        return
    top_df = df.sort_values('Adjusted P-value').head(top_n).copy()
    top_df['-log10(FDR)'] = -np.log10(top_df['Adjusted P-value'])
    if kind == 'KEGG':
        colors = top_df['Adjusted P-value'].apply(lambda x: '#762a83' if x < 0.01 else '#c2a5cf')
        ylabel = 'Pathway'
    else:
        colors = top_df['Adjusted P-value'].apply(lambda x: '#1b7837' if x < 0.01 else '#a6dba0')
        ylabel = 'Term'
    plt.figure(figsize=(9, 0.6*len(top_df)+2))
    bar = sns.barplot(
        x='-log10(FDR)', y=ylabel, data=top_df,
        palette=colors, edgecolor='black'
    )
    for i, (v, n) in enumerate(zip(top_df['-log10(FDR)'], top_df['Overlap'])):
        bar.text(v + 0.05, i, f'n={n}', va='center', fontsize=10, color='black')
    plt.xlabel('-log10(FDR)', fontsize=13, fontweight='bold')
    plt.ylabel(ylabel, fontsize=13, fontweight='bold')
    if title is None:
        title = f'{kind} Top {top_n} Enrichment'
    plt.title(title, fontsize=15, fontweight='bold')
    plt.tight_layout()
    for fmt in image_formats:
        out_path = f"{out_base}.{fmt}"
        plt.savefig(out_path, dpi=150)
        print(f"Saved {out_path}")
    plt.close()

if __name__ == "__main__":
    df = pd.read_csv("/data/pub/zhousha/20260411_MERIPseq/output/exomePeak/con_sig_diff_peak_name.xls", sep="\t")
    gene_list = df["gene_name"].dropna().unique().tolist()
    outdir = "/data/pub/zhousha/20260411_MERIPseq/output/exomePeak/plot/py"
    results = enrich_go_kegg_analysis(gene_list, outdir=outdir, organism="mouse")
    for key, df in results.items():
        kind = 'KEGG' if key == 'KEGG' else 'GO'
        plot_enrichment_bar(df, os.path.join(outdir, f'{key}_top10'), top_n=10,
                            kind=kind, title=None, image_formats=["png"])
