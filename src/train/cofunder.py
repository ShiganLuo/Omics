import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import logging
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,  # 指定输出到 stdout 而不是 stderr
	datefmt='%Y-%m-%d %H:%M:%S'
)


def marker_gene_cofunder_overlap(
        df_log2: pd.DataFrame, 
        gsi_genes: list,
        confounder_gene_sets: dict,
        corr_threshold: float = 0.5,
        allow_high_corr: bool = False,
        min_std: float = 1e-6,
        filter_method: str = "mean"  # "max" 或 "mean"
    ):
    """
    用 GSI gene 与 confounder gene 的 overlap 构建 confounder score，并去除高相关 confounder

    Parameters
    ----------
    df_log2 : pd.DataFrame
        log2(TPM+1) 矩阵，行=基因，列=样本
    gsi_genes : list
        高相关 GSI 基因列表
    confounder_gene_sets : dict
        confounder gene set 字典，key=信号名, value=gene list
    corr_threshold : float
        高相关 confounder 丢弃阈值
    allow_high_corr : bool
        是否允许高相关 confounder 参与回归
    min_std : float
        z-score 计算中最小标准差阈值
    filter_method : str
        高相关 confounder 筛选方式："max" 或 "mean"

    Returns
    -------
    df_resid : pd.DataFrame
        去除 confounder 影响后的 GSI 基因表达 residual
    """

    # -------------------------------
    # 1. 提取 GSI 基因，剔除方差极小基因
    # -------------------------------
    gsi_present = [g for g in gsi_genes if g in df_log2.index]
    gsi_var = [g for g in gsi_present if df_log2.loc[g].std(ddof=0) > min_std]
    if len(gsi_var) == 0:
        raise ValueError("没有满足方差阈值的 GSI 基因")
    df_genes = df_log2.loc[gsi_var]

    # -------------------------------
    # 2. 计算 confounder score（只用 overlap 基因）
    # -------------------------------
    confounder_score = pd.DataFrame(index=df_log2.columns)
    confounder_overlap_dict = {}

    for name, geneset in confounder_gene_sets.items():
        genes_present = [g for g in geneset if g in df_log2.index]
        genes_var = [g for g in genes_present if df_log2.loc[g].std(ddof=0) > min_std]

        # 只保留 overlap 基因
        overlap = list(set(genes_var) & set(gsi_var))
        if len(overlap) == 0:
            continue
        z = df_log2.loc[overlap].apply(lambda x: (x - x.mean()) / (x.std(ddof=0)+min_std), axis=1)
        confounder_score[name] = z.mean(axis=0)
        confounder_overlap_dict[name] = overlap

    if confounder_score.shape[1] == 0:
        logging.info("Warning: 没有有效 confounder 基因，返回原始 GSI 矩阵")
        return df_genes.copy()

    # -------------------------------
    # 3. 丢弃高相关 confounder（可选）
    # -------------------------------
    if not allow_high_corr:
        to_keep = []
        for name in confounder_score.columns:
            corrs = df_genes.apply(lambda x: confounder_score[name].corr(x), axis=1)
            if filter_method == "max":
                corr_val = corrs.abs().max()
            elif filter_method == "mean":
                corr_val = corrs.abs().mean()
            else:
                raise ValueError("filter_method 必须为 'max' 或 'mean'")
            
            if corr_val < corr_threshold:
                to_keep.append(name)
            else:
                logging.info(f"Confounder '{name}' dropped due to high correlation (corr={corr_val:.2f})")
        confounder_score = confounder_score[to_keep]
        if confounder_score.shape[1] == 0:
            logging.info("Warning: All confounders are discarded; the original GSI matrix is returned.")
            return df_genes.copy()

    # -------------------------------
    # 4. 回归 residual
    # -------------------------------
    df_resid = pd.DataFrame(index=df_genes.index, columns=df_genes.columns, dtype=float)
    for gene in df_genes.index:
        y = df_genes.loc[gene].values
        X = confounder_score.values
        X = sm.add_constant(X)
        model = sm.OLS(y, X).fit()
        df_resid.loc[gene] = model.resid

    return df_resid, confounder_overlap_dict


def read_gmt(gmt:str) -> dict:
    gene_sets = {}
    try:
        with open(gmt,'r',encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >=3:
                    set_name = parts[0]
                    genes = parts[2:]
                    gene_sets[set_name] = [g for g in genes if g]
        return gene_sets
    except FileNotFoundError:
        logging.warning(f"file not found: {gmt}")
        return None
    except Exception as e:
        logging.warning(f"some unexpected wrong with reading file: {e} ")
    return None
    
def check_confounder_overlap_corr(df_log2: pd.DataFrame, 
                                       gsi_genes: list, 
                                       confounder_gene_sets: dict,
                                       min_std: float = 1e-6):
    """
    检查 confounder gene 与 GSI gene 的相关性，并记录每个重叠基因与每个 GSI gene 的相关性

    Returns
    -------
    result : pd.DataFrame
        每个 confounder 的:
        - 重叠基因数量
        - 重叠基因列表
        - overlap_gene_corr_dict: {confounder_gene: {gsi_gene: corr, ...}, ...}
        - mean_corr_with_GSI: confounder score 与每个 GSI gene 平均相关性
    """
    # 过滤方差过小的 GSI 基因
    gsi_present = [g for g in gsi_genes if g in df_log2.index]
    gsi_var = [g for g in gsi_present if df_log2.loc[g].std(ddof=0) > min_std]
    if len(gsi_var) == 0:
        raise ValueError("没有满足方差阈值的 GSI 基因")
    df_gsi = df_log2.loc[gsi_var]
    
    result = []

    for name, genes in confounder_gene_sets.items():
        # 只保留存在且方差足够的基因
        genes_present = [g for g in genes if g in df_log2.index]
        genes_var = [g for g in genes_present if df_log2.loc[g].std(ddof=0) > min_std]
        
        overlap = list(set(genes_var) & set(gsi_var))
        n_overlap = len(overlap)
        
        overlap_gene_corr_dict = {}
        if len(genes_var) == 0:
            mean_corr = None
        else:
            # 标准化 confounder genes
            z = df_log2.loc[genes_var].apply(lambda x: (x - x.mean()) / (x.std(ddof=0)+min_std), axis=1)
            score = z.mean(axis=0)
            # confounder score 与每个 GSI gene 相关性
            corrs = df_gsi.apply(lambda x: score.corr(x), axis=1)
            mean_corr = corrs.mean()
            
            # 每个 overlap 基因与每个 GSI gene 的相关性
            for g in overlap:
                g_expr = df_log2.loc[g]
                g_corrs = df_gsi.apply(lambda x: g_expr.corr(x), axis=1).to_dict()
                overlap_gene_corr_dict[g] = g_corrs

        result.append({
            'confounder': name,
            'n_overlap': n_overlap,
            'overlap_genes': overlap,
            'overlap_gene_corr': overlap_gene_corr_dict,
            'mean_corr_with_GSI': mean_corr
        })
        
    return pd.DataFrame(result)



def plot_confounder_gsi_corr_heatmap(df_corr_result: pd.DataFrame,outplot:str, figsize=(12,8), cmap='RdBu_r'):
    """
    将 confounder gene 与 GSI gene 的相关性绘制热图

    Parameters
    ----------
    df_corr_result : pd.DataFrame
        check_confounder_overlap_corr_full 的结果
    figsize : tuple
        图像大小
    cmap : str
        colormap
    """
    # 构建 heatmap 矩阵
    heatmap_data = {}
    confounder_labels = []

    for idx, row in df_corr_result.iterrows():
        conf_name = row['confounder']
        overlap_dict = row['overlap_gene_corr']  # {conf_gene: {gsi_gene: corr}}
        for conf_gene, gsi_corrs in overlap_dict.items():
            heatmap_data[conf_gene] = gsi_corrs
            confounder_labels.append(conf_name)
    
    # 转为 DataFrame，行=confounder gene, 列=GSI gene
    heatmap_df = pd.DataFrame.from_dict(heatmap_data, orient='index')
    
    # 按 confounder 排序
    heatmap_df['confounder'] = confounder_labels
    heatmap_df.sort_values('confounder', inplace=True)
    confounder_labels_sorted = heatmap_df['confounder'].values
    heatmap_df.drop(columns=['confounder'], inplace=True)
    
    # 绘制热图
    plt.figure(figsize=figsize)
    sns.heatmap(heatmap_df, cmap=cmap, center=0, linewidths=0.5)
    plt.xlabel('GSI genes')
    plt.ylabel('Confounder genes')
    
    # 在 y 轴标记 confounder 来源
    y_ticks = np.arange(0.5, heatmap_df.shape[0]+0.5)
    plt.yticks(y_ticks, heatmap_df.index)
    
    # 可以加 colorbar title
    plt.title('Confounder gene vs GSI gene correlation')
    plt.tight_layout()
    plt.savefig(outplot,dpi=300)



if __name__ == "__main__":
    confounder_gene_sets = {
    'prolif': ['MKI67', 'PCNA', 'TOP2A',
                'MCM2', 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MCM7',
                'CCNB1', 'CCNB2', 'CDK1',
                'CDC20', 'BUB1', 'AURKA', 'TTK', 'RRM2', 'GMNN'
                ],
    'senescence': ['CDKN2A', 'CDKN1A', 'IL6', 'IL8',
                    'TNFRSF10B', 'SERPINE1', 'IGFBP3',
                    'LMNB1', 'MMP3'
                    ],
    'quiescence': ['CCND1', 'RB1', 'FOXO3', 'CDKN1B', 'TSC2'
                    ]
    }
    df_log2 = pd.read_csv("/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229.tsv",sep="\t",index_col=0)
    gsi_gene_dict = read_gmt("/home/luosg/Data/genomeStability/data/final.gmt")
    gsi_genes = list({gene for genes in gsi_gene_dict.values() for gene in genes})
    logging.info((f"总共有 {len(gsi_genes)} 个非重复gsi基因"))
    df_resid,confounder_overlap_dict = marker_gene_cofunder_overlap(df_log2=df_log2,
                                     gsi_genes=gsi_genes,
                                     confounder_gene_sets=confounder_gene_sets,
                                     corr_threshold=0.5)
    df_resid.to_csv("/home/luosg/Data/genomeStability/output/result/matrix/log2tpm229_514_resid.tsv",sep="\t")
    # df_check = check_confounder_overlap_corr(df_log2, gsi_genes, confounder_gene_sets)
    # df_check.to_csv("a.tsv",sep="\t")
    # plot_confounder_gsi_corr_heatmap(df_check,"/home/luosg/Data/genomeStability/output/result/cofounder/cofounder_gsi_heatmap.png")
