from itertools import combinations
from statannotations.Annotator import Annotator
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
import itertools
import numpy as np
def violin(
        df:pd.DataFrame,
        outplot:str
):
    unique_cell_types = df['Cell_type'].unique().tolist()
    box_pairs = list(combinations(unique_cell_types, 2))
    # combinations(unique_cell_types, 2) 会生成所有不重复的 (组A, 组B) 组合

    plt.figure(figsize=(10, 8))
    ax = sns.violinplot(
        x='Cell_type',
        y='LV1',
        data=df,
        palette='Set2',
        inner='box',
        linewidth=1.5
    )
    annotator = Annotator(
        ax,
        box_pairs,
        data=df,
        x='Cell_type',
        y='LV1',
        order=unique_cell_types # 确保 order 与图中的顺序一致
    )
    # 注意：进行多重比较时，通常需要进行多重检验校正 (Multiple Testing Correction)。
    annotator.configure(
        test='Mann-Whitney',
        text_format='star',
        loc='inside',
        # 多重检验校正：例如使用 Benjamini/Hochberg FDR
        comparisons_correction='fdr_bh', # <--- 关键参数！
        pvalue_thresholds=[
            [0.05, "*"],
            [0.01, "**"],
            [0.001, "***"],
            [0.0001, "****"]
        ]
    )

    annotator.apply_and_annotate()
    ax.set_xlabel('') 
    plt.ylabel('GSI')
    plt.tight_layout()
    plt.savefig(outplot,dpi=300)



def violin_by_method_celltype_step(df_dict: dict, outplot: str):
    all_dfs = []
    for method, df in df_dict.items():
        tmp = df.copy()
        tmp['Method'] = method
        all_dfs.append(tmp)
    combined_df = pd.concat(all_dfs, axis=0)

    methods = combined_df['Method'].unique().tolist()
    cell_types = combined_df['Cell_type'].unique().tolist()
    palette = sns.color_palette("Set2", n_colors=len(cell_types))
    cell_color_map = dict(zip(cell_types, palette))

    plt.figure(figsize=(12, 8))
    ax = sns.violinplot(
        x='Method',
        y='LV1',
        hue='Cell_type',
        data=combined_df,
        palette=cell_color_map,
        inner='box',
        linewidth=1.5,
        dodge=True
    )
    # 计算每个方法+cell_type的小提琴中心点
    violin_positions = {}
    for method_idx, method in enumerate(methods):
        cells_in_method = combined_df[combined_df['Method']==method]['Cell_type'].unique()
        n_cell_types = len(cells_in_method)
        offsets = np.linspace(-0.2, 0.2, n_cell_types)
        for i, ct in enumerate(cells_in_method):
            violin_positions[(method, ct)] = method_idx + offsets[i]

    # 遍历每个方法做两两检验
    for method in methods:
        subset = combined_df[combined_df['Method'] == method]
        cells_in_method = subset['Cell_type'].unique().tolist()
        if len(cells_in_method) < 2:
            continue
        pairs = list(itertools.combinations(cells_in_method, 2))
        step = subset['LV1'].max() * 0.05

        for i, (c1, c2) in enumerate(pairs):
            data1 = subset[subset['Cell_type']==c1]['LV1']
            data2 = subset[subset['Cell_type']==c2]['LV1']
            stat, p = mannwhitneyu(data1, data2, alternative='two-sided')

            x1 = violin_positions[(method, c1)]
            x2 = violin_positions[(method, c2)]
            y_base = max(data1.max(), data2.max())
            y = y_base + (i+1)*step
            y_top = y + step*0.3  # 梯形顶部高度

            # 画梯形：竖直-水平-竖直
            ax.plot([x1, x1], [y, y_top], color='black', lw=2)
            ax.plot([x1, x2], [y_top, y_top], color='black', lw=2)
            ax.plot([x2, x2], [y_top, y], color='black', lw=2)

            # 显示星号
            ax.text((x1+x2)/2, y_top + step*0.05,
                    "****" if p<0.0001 else
                    "***" if p<0.001 else
                    "**" if p<0.01 else
                    "*" if p<0.05 else "ns",
                    ha='center', va='bottom', fontsize=14)

    ax.set_xlabel('Method',fontsize=16)
    ax.set_ylabel('LV1',fontsize=16)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.legend(fontsize=13)
    plt.tight_layout()
    plt.savefig(outplot, dpi=300)
    plt.close()
    print(f"图已保存到 {outplot}")

def combine_run():
    Latent = {
        "pca": "/disk5/luosg/Totipotent20251031/RNAseqML/latent/fa_factors.csv",
        "nmf": "/disk5/luosg/Totipotent20251031/RNAseqML/latent/nmf_factors.csv",
        "factor": "/disk5/luosg/Totipotent20251031/RNAseqML/latent/pca_factors.csv"
    }
    outdir="/disk5/luosg/Totipotent20251031/RNAseqML/latent"
    df_dict = {}
    for method, infile in Latent.items():
        df = pd.read_csv(infile, index_col=0)
        df_annotation = pd.read_csv("/disk5/luosg/Totipotent20251031/data/RNAseqML/target_fq.tsv", sep="\t")
        df_annotation.drop_duplicates(subset=["Sample_id"], inplace=True)
        Sample_to_Status_map = df_annotation.set_index('Sample_id')['Status']
        new_index = df.index.map(Sample_to_Status_map)
        df.index = new_index
        df_scaled = (df - df.min()) / (df.max() - df.min())
        df_scaled = df_scaled.reset_index(names="Cell_type")
        df_dict[method] = df_scaled

    outplot = f"{outdir}/all_methods_violin.png"
    violin_by_method_celltype_step(df_dict, outplot)
if __name__ == "__main__":
    # single_run()
    combine_run()
