import pandas as pd
import numpy as np
from scipy.stats import poisson
import matplotlib.pyplot as plt
import seaborn as sns
import os
from typing import Dict,Set
def gsi_mutation_enrichment(gene_count_file, GSI_gene, outdir='GSI_enrichment'):
    """
    利用全能干细胞特有变异，比较 GSI 基因突变数与背景基因突变数的差异。
    使用 Poisson 检验。
    
    Args:
        gene_count_file (str): 基因突变计数 CSV 文件，行索引为基因名，列为突变数
        GSI_gene (list or set): GSI 基因列表
        outdir (str): 输出结果文件夹
        
    Returns:
        dict: {'gsi_total': , 'expected_total': , 'p_value': , 'gsi_mean': , 'bg_mean': }
    """
    
    os.makedirs(outdir, exist_ok=True)
    
    # 读取基因突变计数
    gene_count = pd.read_csv(gene_count_file, index_col=0)
    gene_count.columns = ['mut_count']
    
    # 提取 GSI 和背景基因突变数
    GSI_gene = set(GSI_gene)
    gsi_mut_counts = gene_count.loc[gene_count.index.isin(GSI_gene), 'mut_count']
    bg_mut_counts = gene_count.loc[~gene_count.index.isin(GSI_gene), 'mut_count']
    
    # 背景平均突变数
    bg_mean = bg_mut_counts.mean()
    
    # GSI 总突变数及预期
    gsi_total = gsi_mut_counts.sum()
    gsi_num_genes = len(gsi_mut_counts)
    expected_total = bg_mean * gsi_num_genes
    
    # Poisson 右尾检验
    p_value = 1 - poisson.cdf(gsi_total - 1, expected_total)
    
    # 可视化
    df_plot = pd.DataFrame({
        'Gene_Set': ['GSI', 'Background'],
        'Mean_Mut_Count': [gsi_mut_counts.mean(), bg_mean]
    })
    
    plt.figure(figsize=(6,4))
    sns.barplot(x='Gene_Set', y='Mean_Mut_Count', data=df_plot, palette=['skyblue','lightgray'])
    plt.ylabel('Mean mutation count per gene')
    plt.title('GSI genes vs Background genes mutation count')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'GSI_vs_Background.png'))
    plt.close()
    
    # 输出统计结果
    result_dict = {
        'gsi_total': gsi_total,
        'expected_total': expected_total,
        'p_value': p_value,
        'gsi_mean': gsi_mut_counts.mean(),
        'bg_mean': bg_mean
    }
    
    pd.DataFrame([result_dict]).to_csv(os.path.join(outdir, 'GSI_enrichment_results.csv'), index=False)
    
    print(f"GSI 突变富集分析完成，结果保存在 {outdir}")
    return result_dict

def read_gmt_file(gmt_file_path: str) -> Dict[str, Set[str]]:
    """
    读取 GMT 文件并将其内容解析为基因集字典。

    Args:
        gmt_file_path: GMT 文件的完整路径。

    Returns:
        一个字典，键为基因集名称 (Gene Set Name)，值为该基因集包含的基因集合 (Set of Genes)。
    """
    gene_sets: Dict[str, Set[str]] = {}
    
    try:
        with open(gmt_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    gene_set_name = parts[0]
                    genes = set(parts[2:])
                    
                    genes.discard('') 
                    
                    if genes:
                        gene_sets[gene_set_name] = genes
                        
        print(f"成功读取 {len(gene_sets)} 个基因集。")
        return gene_sets
        
    except FileNotFoundError:
        print(f"错误：文件未找到，请检查路径: {gmt_file_path}")
        return {}
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return {}

def GSI_gene(
        gmt_file:str
):
    all_gene_sets = read_gmt_file(gmt_file)
    total_gmt_genes = set()
    for gene_set in all_gene_sets.values():
        total_gmt_genes.update(gene_set)

    return total_gmt_genes

if __name__ == "__main__":
    human_files = ["/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/ci8CLC_gene_counts.csv",
                   "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/hTBLC/hTBLC_gene_counts.csv"]
    mouse_files = ["/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/TLSC/TLSC_gene_counts.csv",
                   "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ciTotiSC/ciTotiSC_gene_counts.csv"]
    human_gmt = "/disk5/luosg/Totipotent20251031/data/geneset/GSI_human.gmt"
    mouse_gmt = "/disk5/luosg/Totipotent20251031/data/geneset/GSI_mouse.gmt"
    genes_gsi = GSI_gene(human_gmt)
    gsi_mutation_enrichment("/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/ci8CLC/ci8CLC_gene_counts.csv",genes_gsi,"./")

    # for infile in human_files:
    #     genes_gsi = GSI_gene(human_gmt)
    #     gsi_mutation_enrichment(infile,genes_gsi,"./")
