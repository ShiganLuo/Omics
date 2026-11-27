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
    cmap = plt.colormaps.get_cmap("tab20") 
    unique_funcs = df["Func"].unique()
    color_map = {func: cmap(i / len(unique_funcs)) for i, func in enumerate(unique_funcs)}

    groups = df["group"].unique()
    num_groups = len(groups)
    fig, axes = plt.subplots(5, 2, figsize=(9 * 2, 6 * 5),dpi=300)
    axes_flat = axes.flatten()
    if num_groups == 1:
        axes = [axes]

    def autopct_format(pct):
        return f'{pct:.1f}%' if pct >= 5 else ''  
    for ax, group in zip(axes_flat, groups):
        subset = df[df["group"] == group]
        func_counts = subset.groupby("Func")["Number"].sum()

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
        ax.set_ylabel("")
    
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
    plotIntervalDistribution(dfTotal,"a.png")
