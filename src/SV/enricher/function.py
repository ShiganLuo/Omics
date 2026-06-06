import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import logging
from typing import List, Optional
logger = logging.getLogger(__name__)

def enrich_go(
        gene_list: List[str],
        organism: str = 'Human',
        outdir: str = 'enrichment_results',
        top_n: int = 10,
        cutoff: float = 0.05,
        image_formats: Optional[List[str]] = None
    ) -> dict:
    """
    Perform GO enrichment analysis and visualize top enriched terms.

    Runs Gene Ontology (Biological Process, Cellular Component, Molecular
    Function) enrichment via gseapy's ``enrichr`` and saves bar plots of
    the top-ranked terms.

    Parameters
    ----------
    gene_list : list of str
        Gene symbols to analyze for enrichment.
    organism : str, default="Human"
        Organism name (``"Human"`` or ``"Mouse"``).
    outdir : str, default="enrichment_results"
        Directory to save enrichment results and plots.
    top_n : int, default=10
        Number of top enriched terms to visualize in bar plots.
    cutoff : float, default=0.05
        Adjusted p-value cutoff for significance.
    image_formats : list of str, optional
        Output image formats (e.g. ``["png", "pdf"]``). Defaults to
        ``["png"]``.

    Returns
    -------
    dict
        Mapping of GO category name to its enrichment results DataFrame.
        Keys: ``"GO_Biological_Process_2021"``, ``"GO_Cellular_Component_2021"``,
        ``"GO_Molecular_Function_2021"``.
    """
    if image_formats is None:
        image_formats = ['png']
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
            cutoff=cutoff
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
            for fmt in image_formats:
                plt.savefig(os.path.join(outdir, f'{cat}_top{top_n}.{fmt}'))
            plt.close()
    
def enrich_kegg(
        gene_list: List[str],
        organism: str = 'Human',
        outdir: str = 'enrichment_results',
        top_n: int = 10,
        cutoff: float = 0.05,
        image_formats: Optional[List[str]] = None
    ) -> dict:
    """
    Perform KEGG pathway enrichment analysis and visualize top terms.

    Runs KEGG pathway enrichment via gseapy's ``enrichr`` and saves bar
    plots of the top-ranked pathways.

    Parameters
    ----------
    gene_list : list of str
        Gene symbols to analyze for enrichment.
    organism : str, default="Human"
        Organism name (``"Human"`` or ``"Mouse"``).
    outdir : str, default="enrichment_results"
        Directory to save enrichment results and plots.
    top_n : int, default=10
        Number of top enriched pathways to visualize in bar plots.
    cutoff : float, default=0.05
        Adjusted p-value cutoff for significance.
    image_formats : list of str, optional
        Output image formats (e.g. ``["png", "pdf"]``). Defaults to
        ``["png"]``.

    Returns
    -------
    dict
        Mapping of ``"KEGG"`` to its enrichment results DataFrame.
    """
    if image_formats is None:
        image_formats = ['png']
    results_dict = {}
    kegg_cat = 'KEGG_2021_Human' if organism.lower()=='human' else 'KEGG_2021_Mouse'
    enr_kegg = gp.enrichr(
        gene_list=gene_list,
        gene_sets=kegg_cat,
        organism=organism,
        description='KEGG',
        outdir=outdir,
        cutoff=cutoff
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
        for fmt in image_formats:
            plt.savefig(os.path.join(outdir, f'KEGG_top{top_n}.{fmt}'))
        plt.close()
    
    for key, df in results_dict.items():
        df.to_csv(os.path.join(outdir, f'{key}_enrichment.csv'), index=False)
    
    logger.info(f"富集分析完成，结果保存在文件夹: {outdir}")
    return results_dict

if __name__ == "__main__":
    # --- 使用示例 ---
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/ci8CLC_gene_counts.csv",index_col=0)
    gene_list = df.index.to_list()
    results = enrich_go_kegg(gene_list,outdir="/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/go")
