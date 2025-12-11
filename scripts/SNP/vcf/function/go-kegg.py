import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

def enrich_go_kegg(gene_list, organism='Human', outdir='enrichment_results', top_n=10):
    """
    输入基因列表，做 GO 和 KEGG 富集分析，并输出结果和柱状图。
    
    Args:
        gene_list (list): 基因名列表，例如 ['TP53','BRCA1',...]
        organism (str): 'Human' 或 'Mouse'
        outdir (str): 输出结果的文件夹
        top_n (int): 绘制柱状图显示的前 top N 项
        
    Returns:
        dict: {'GO_BP': df, 'GO_CC': df, 'GO_MF': df, 'KEGG': df}
    """
    
    os.makedirs(outdir, exist_ok=True)
    results_dict = {}
    
    # --- GO 富集分析 ---
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
        
        # 可视化 top_n
        if not df.empty:
            top_df = df.sort_values('Adjusted P-value').head(top_n)
            plt.figure(figsize=(8,6))
            sns.barplot(x=-np.log10(top_df['Adjusted P-value']), y=top_df['Term'], color='skyblue')
            plt.xlabel('-log10(FDR)')
            plt.title(f'{cat} Top {top_n} Enrichment')
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f'{cat}_top{top_n}.png'))
            plt.close()
    
    # --- KEGG 富集分析 ---
    kegg_cat = 'KEGG_2021_Human' if organism.lower()=='human' else 'KEGG_2021_Mouse'
    enr_kegg = gp.enrichr(
        gene_list=gene_list,
        gene_sets=kegg_cat,
        organism=organism,
        description='KEGG',
        outdir=outdir,
        cutoff=0.05
    )
    df_kegg = enr_kegg.results
    results_dict['KEGG'] = df_kegg
    
    if not df_kegg.empty:
        top_df = df_kegg.sort_values('Adjusted P-value').head(top_n)
        plt.figure(figsize=(8,6))
        sns.barplot(x=-np.log10(top_df['Adjusted P-value']), y=top_df['Term'], color='lightgreen')
        plt.xlabel('-log10(FDR)')
        plt.title(f'KEGG Top {top_n} Enrichment')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f'KEGG_top{top_n}.png'))
        plt.close()
    
    # 保存结果表
    for key, df in results_dict.items():
        df.to_csv(os.path.join(outdir, f'{key}_enrichment.csv'), index=False)
    
    print(f"富集分析完成，结果保存在文件夹: {outdir}")
    return results_dict

if __name__ == "__main__":
    # --- 使用示例 ---
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/ci8CLC_gene_counts.csv",index_col=0)
    gene_list = df.index.to_list()
    results = enrich_go_kegg(gene_list,outdir="/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/go")
