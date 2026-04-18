import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from typing import List, Dict, Literal


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
    # 保存结果表
    for key, df in results_dict.items():
        df.to_csv(os.path.join(outdir, f'{key}_enrichment.csv'), index=False)
    print(f"富集分析完成，结果保存在文件夹: {outdir}")
    return results_dict

def plot_enrichment_bar(df: pd.DataFrame, outpath: str, top_n=10, kind='GO', title=None):
    """
    只负责绘制GO/KEGG富集柱状图。
    """
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
    plt.savefig(outpath, dpi=150)
    plt.close()

if __name__ == "__main__":
    # --- 使用示例 ---
    df = pd.read_csv("/data/pub/zhousha/20260411_MERIPseq/output/exomePeak/con_sig_diff_peak_name.xls", sep="\t")
    gene_list = df["gene_name"].dropna().unique().tolist()
    outdir = "/data/pub/zhousha/20260411_MERIPseq/output/exomePeak/plot/py"
    results = enrich_go_kegg_analysis(gene_list, outdir=outdir, organism="mouse")
    # 分别绘图
    for key, df in results.items():
        kind = 'KEGG' if key == 'KEGG' else 'GO'
        plot_enrichment_bar(df, os.path.join(outdir, f'{key}_top10.png'), top_n=10, kind=kind, title=None)
