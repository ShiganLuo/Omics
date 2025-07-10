import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import re
import argparse
from functools import reduce

def plotIntervalDistribution(df:pd.DataFrame,outfile:str):
    # 要绘制的组别
    # 颜色映射
    cmap = plt.colormaps.get_cmap("tab20") 
    unique_funcs = df["Func"].unique()
    color_map = {func: cmap(i / len(unique_funcs)) for i, func in enumerate(unique_funcs)}

    groups = df["group"].unique()
    num_groups = len(groups)
    # 设置图像布局
    fig, axes = plt.subplots(5, 2, figsize=(9 * 2, 6 * 5),dpi=300)
    axes_flat = axes.flatten()
    if num_groups == 1:
        axes = [axes]  # 统一处理 axes 变量

    def autopct_format(pct):
        return f'{pct:.1f}%' if pct >= 5 else ''  
    # 遍历不同的 group
    for ax, group in zip(axes_flat, groups):
        subset = df[df["group"] == group]
        func_counts = subset.groupby("Func")["Number"].sum()  # 按 Func 统计 Number 总和

        # 画饼图并存储返回值
        colors = [color_map[func] for func in func_counts.index]
        wedges, texts, autotexts = ax.pie(
            func_counts, 
            autopct=autopct_format,
            startangle=140, 
            colors=colors, 
            wedgeprops={'edgecolor': 'black'},
            labeldistance=1.1,
            textprops={'fontsize': 17},
            radius = 1.2
        )
        if group == 'mESC1' or group == 'mESC2':
            group = 'mESC'
        ax.set_title(f"Functional Distribution\n{group}",fontsize = 20)
        ax.set_ylabel("")  # 去除 y 轴标签
        # 添加图例
    # fig.legend(
    #     wedges, func_counts.index, 
    #     title="Func", 
    #     loc="center",  # 图例居中
    #     bbox_to_anchor=(0.5, 0.5),  # 调整位置
    #     ncol=1,  # 控制列数
    #     fontsize=20
    # )
    # plt.subplots_adjust(
    #     left=0.05,    # 左侧边距
    #     right=0.95,   # 右侧边距
    #     bottom=0.05,  # 底部边距
    #     top=0.95,     # 顶部边距
    #     wspace=0.1,   # 水平间距
    #     hspace=0.3    # 垂直间距
    # )
    
    # 调整布局
    plt.tight_layout()
    fig.savefig(outfile)

if __name__ == '__main__':
    infiles = ["/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/SNP/IntervalDistribution.csv",
                "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE183522/SNP/IntervalDistribution.csv",
                "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE185005/SNP/IntervalDistribution.csv",
                "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE204801/SNP/IntervalDistribution.csv",
                "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE224794/SNP/IntervalDistribution.csv"]
    dfAll = []
    for infile in infiles:
        df = pd.read_csv(infile,sep=",")
        print(df['group'].unique())
        dfAll.append(df)
    dfTotal = reduce(lambda left, right: pd.concat([left,right],axis=0),dfAll)
    # print(dfTotal.head())
    plotIntervalDistribution(dfTotal,"a.png")
    # dfTotal = reduce(lambda left, right: pd.merge(left, right, on = merge_on, how='outer'), df_list)