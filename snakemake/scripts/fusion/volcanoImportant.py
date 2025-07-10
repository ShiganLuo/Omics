import pandas as pd
import re
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def extract_with_pandas(gtf_path: str):
    genes = dict()
    
    # 使用 pandas 读取 GTF 文件，跳过注释行
    gtf_df = pd.read_csv(gtf_path, sep='\t', comment='#', header=None,dtype={5: str})
    
    # 设定 GTF 文件的列名
    gtf_df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    
    # 筛选出 'gene' 类型的行
    genes_df = gtf_df[gtf_df['feature'] == 'transcript']
    
    # 提取 gene_id 和 gene_name
    for _, row in genes_df.iterrows():
        attributes = row['attribute']
        # 提取 gene_id 和 gene_name
        gene_id = next((x.split('"')[1] for x in attributes.split(';') if 'gene_id' in x), None)
        gene_name = next((x.split('"')[1] for x in attributes.split(';') if 'ref_gene_id' in x), None)
        if gene_id and gene_name:
            genes[gene_id] = gene_name
    
    return genes


def filterTE(infile:str,TE:str):
    TEfusionGeneSet = set()
    with open(infile, 'r', encoding='utf-8') as file:
        for line in file:
            if re.search(TE, line):
                gene = line.split("\t")[3]
                TEfusionGeneSet.add(gene)
    return TEfusionGeneSet

def volcano(infile:str,annFile:str,outImage:str,Title:str,TEfusionGeneSet:set=None,padjThreshold:int = 0.05,l2fcThreshold:int = 0.58):
    df1 = pd.read_table(infile,sep = "\t",index_col=0)
    df1  = df1.sort_values('padj',ascending=True)
    df2 = pd.read_table(annFile,sep = "\t")
    df = pd.merge(
        df1, 
        df2,
        left_on=df1.index, 
        right_on=df2.columns[0],
        how='inner'  # 连接方式：inner, left, right, outer
    )
    # print(df.head())
    # 计算-log10(padj)
    df['-log10(padj)'] = -np.log10(df['padj'])
    # 设置颜色
    ### 方案一：常规
    # df['color'] = ['red' if (p < padjThreshold and fc > l2fcThreshold) else 
    #             'blue' if (p < padjThreshold and fc < -l2fcThreshold) else
    #             'green' if (abs(fc) > l2fcThreshold) else 'grey'
    #             for p, fc in zip(df['padj'], df['log2FoldChange'])]
    ###方案二：融合基因
    df['color'] = ['red' if (gene in TEfusionGeneSet) else
                   'grey'
                   for gene in df['gene_id']]
    # 绘制火山图
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=df, x='log2FoldChange', y='-log10(padj)', 
                    hue='color',palette={'red':'red', 'blue':'blue', 'green':'green', 'grey':'grey'},
                    s=8, alpha=0.7)
    # 添加阈值线
    plt.axhline(y=-np.log10(padjThreshold), color='black', linestyle='--', linewidth=0.8)
    plt.axvline(x=l2fcThreshold, color='black', linestyle='--', linewidth=0.8)
    plt.axvline(x=-l2fcThreshold, color='black', linestyle='--', linewidth=0.8)

    # # 标记部分显著基因
    # significant = df[(df['padj'] < padjThreshold) & (abs(df['log2FoldChange']) > l2fcThreshold)].head(5)
    # for i, row in significant.iterrows():
    #     plt.text(row['log2FoldChange'], row['-log10(padj)'], row['gene_name'], 
    #             fontsize=8, ha='center', va='bottom')


    plt.xlabel('log2(Fold Change)', fontsize=12)
    plt.ylabel('-log10(adjusted p-value)', fontsize=12)
    # plt.title('Volcano Plot of Differential Expression', fontsize=14)
    plt.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', label=f'{Title}-fusion gene',
                markerfacecolor='red', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='other genes',
                markerfacecolor='grey', markersize=8),
    ], loc='upper right')

    plt.tight_layout()
    plt.savefig(outImage, dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="volcano plot for TE fusion gene")
    parser.add_argument('--input', type=str, required=True, help='Path to StgTEOverlap.bed file')
    parser.add_argument('--output', type=str, required=True, help='Path to output Image file')
    parser.add_argument('--gtf', type=str, required=True, help='Path to gtf file of stringTie')
    parser.add_argument('--volcano', type=str, required=True, help='Path to the results of DESeq2,not filter')
    parser.add_argument('--annFile', type=str, required=True, help='Path to the annotation for gene_id')
    parser.add_argument('--TE', type=str, required=True,nargs='+',help='Path to the annotation for gene_id')
    parser.add_argument('--graphTitle', type=str, required=True,help='Path to the annotation for gene_id')
    args = parser.parse_args()
    fusionFile = args.input
    outImage = args.output
    gtf = args.gtf
    volcanoFile = args.volcano
    annFile = args.annFile
    targetTE = args.TE
    graphTitle = args.graphTitle
    # gene_dict = extract_with_pandas('/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/2pass/mouse.gtf')
    # # print(gene_dict)
    # fusionFile = "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/2pass/mouse_StgTEOverlap.bed"
    # MERVL = set.union(filterTE(fusionFile,"MERVL-int"),filterTE(fusionFile,"MT2_Mm"))
    # # print(MERVL)
    # fusionGene = {gene_dict[gene_id] for gene_id in  MERVL if gene_id in gene_dict}
    # # print(fusionGene)

    # volcanoFile = "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/DESeq2/TEcount_Gene.csv"
    # annFile = "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv"
    # outImage = 'volcano_plot.png'
    # volcano(volcanoFile,annFile,outImage,TEfusionGeneSet=fusionGene)
    gene_dict = extract_with_pandas(gtf)
    # print(gene_dict)
    TEBed = set() 
    for i in targetTE:
        TEBed = set.union(TEBed,filterTE(fusionFile,i))

    # print(MERVL)
    fusionGene = {gene_dict[gene_id] for gene_id in  TEBed if gene_id in gene_dict}
    # print(fusionGene)
    volcano(volcanoFile,annFile,outImage,Title=graphTitle,TEfusionGeneSet=fusionGene)

