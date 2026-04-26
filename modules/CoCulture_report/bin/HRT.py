import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging
import os
from typing import List, Literal
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)
def Double_bar(df: pd.DataFrame, outImage: str, fig_width: float = 12, fig_height: float = 7):
    """
    Plot gene and TE counts as double bar plot
    df: columns = feature, count, group
    fig_width, fig_height: 图像宽高，可控制
    """
    import matplotlib.pyplot as plt
    import numpy as np

    groups = df["group"].unique()
    colors_map = {"gene": plt.cm.Dark2, "TE": plt.cm.Pastel1}

    total_bars = sum(df[df["group"] == g].shape[0] for g in groups)
    # 最小宽度保证少量数据也不太窄
    fig_width = max(fig_width, total_bars * 0.5)
    fig_height = max(fig_height, 6)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    start = 0
    xticks = []
    xlabels = []

    y_max = df["count"].max() * 1.2

    for group in groups:
        group_df = df[df["group"] == group].reset_index(drop=True)
        n = len(group_df)
        if n == 0:
            continue

        colors = colors_map[group](np.linspace(0, 1, n))
        for i, (feat, cnt) in enumerate(zip(group_df["feature"], group_df["count"])):
            ax.bar(start + i, cnt, color=colors[i])
            # 相对 y_max 偏移，自动适应
            ax.text(start + i, cnt + 0.02 * y_max, f"{cnt:.0f}", ha="center", va="bottom", fontsize=9,rotation=60)

        xticks.extend(range(start, start + n))
        xlabels.extend(group_df["feature"])

        # group 横线，文字用 transform 固定在轴坐标系
        ax.hlines(y=y_max * 0.95, xmin=start, xmax=start + n - 1, color="black", linewidth=1.5)
        ax.text(start + (n - 1)/2, 0.96, group, ha="center", va="bottom",
                fontsize=12, transform=ax.get_xaxis_transform())

        start += n + 2  # 每组间隔

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=10)

    ax.set_ylim(0, y_max)
    ax.set_ylabel("Average Count", fontsize=12)

    # 动态底部留白，防止 xlabels 被遮挡
    bottom_margin = max(0.2, min(0.3, 0.05 * total_bars))
    plt.subplots_adjust(bottom=bottom_margin)

    plt.tight_layout()
    plt.savefig(outImage, dpi=300)
    plt.close()

def process_group(
        df:pd.DataFrame, 
        samples:List, 
        min_count:int=1,
        top_n:int=10,
        extract_te_index:Literal[0,1,2,3]=0
    ) -> pd.DataFrame:
    """
    Function: process data for one group, filter by min_count, select top_n features
    df should have columns: feature, sample1, sample2, ...
     - samples: list of sample columns to consider
     - min_count: minimum count threshold for all samples
     - top_n: number of top features to select based on average count
     - extract_te_index: if >0, extract the specified index from TE feature name (split by ":") as feature name
    Return: DataFrame with columns: feature, count, group
     - feature: gene name or TE family
     - count: average count across samples
     - group: "gene" or "TE"
    """
    logger.info(f"Processing group with samples: {samples}, min_count: {min_count}, top_n: {top_n}")
    df = df.copy()
    df["is_TE"] = df.iloc[:,0].str.contains(r":")  # 含冒号视为TE
    df.to_csv("debug_input.csv", index=False)  # 输出中间文件，便于调试
    # 先处理基因
    gene_df = df[~df["is_TE"]].copy()
    gene_df["feature"] = gene_df.iloc[:,0]
    gene_df = gene_df[samples + ["feature"]]
    gene_df = gene_df[gene_df[samples].apply(lambda x: all(x > min_count), axis=1)]
    gene_df["count"] = gene_df[samples].mean(axis=1)
    gene_df["group"] = "gene"
    gene_df = gene_df[["feature","count","group"]]
    gene_df = gene_df.sort_values("count", ascending=False).head(top_n)


    # 处理TE
    te_df = df[df["is_TE"]].copy()
    te_df["feature"] = te_df.iloc[:,0].str.split(":").str[extract_te_index]
    te_df = te_df[samples + ["feature"]]
    te_df = te_df.groupby("feature")[samples].sum().reset_index()  # 同亚家族合并
    te_df = te_df[te_df[samples].apply(lambda x: all(x > min_count), axis=1)]
    te_df["count"] = te_df[samples].mean(axis=1)
    te_df["group"] = "TE"
    te_df = te_df[["feature","count","group"]]
    te_df = te_df.sort_values("count", ascending=False).head(top_n)

    return pd.concat([gene_df, te_df], ignore_index=True)

def prepare_data(
        infile1:str, 
        infile2:str, 
        outdir:str, 
        genome1_samples:List[str]=None, 
        genome2_samples:List[str]=None, 
        min_count:int=1, 
        top_n:int=10,
        mode:Literal["TEcount","TElocal"]="TEcount",
        level:Literal["loci","subfamily","family","class"]="subfamily"
    ):
    """
    Function: prepare data for plotting, save intermediate files, call plotting function
     - infile1: TEtranscripts output file for genome1
     - infile2: TEtranscripts output file for genome2
     - outdir: output directory for intermediate files and plots
     - genome1_samples: list of sample names for genome1
     - genome2_samples: list of sample names for genome2
     - min_count: minimum count threshold for filtering features
     - top_n: number of top features to select for plotting
     - mode: type of TEtranscripts output file, affects how TE feature names are processed
    Return: None (saves intermediate files and plots to outdir)
     - intermediate files: genome1_data.csv, genome2_data.csv
     - plots: genome1.png, genome2.png
    Note:
     - mode determines how TE feature names are processed. TEcount files typically have format family:class:superfamily, while TElocal files have format family:class:superfamily:location. The extract_te_index parameter in process_group function is used to specify which part of the TE feature name to extract as the TE family name for plotting.
    """
    df1 = pd.read_csv(infile1, sep="\t")
    df2 = pd.read_csv(infile2, sep="\t")
    genome1 = os.path.basename(os.path.dirname(infile1))
    genome2 = os.path.basename(os.path.dirname(infile2))
    logger.info(f"Loaded data for {genome1} count: {infile1}, rows: {len(df1)}")
    logger.info(f"Loaded data for {genome2} count: {infile2}, rows: {len(df2)}")
    if mode == "TEcount":
        if level == "subfamily":
            extract_te_index = 0
        elif level == "family":
            extract_te_index = 1
        elif level == "class":
            extract_te_index = 2
        else:
            logger.error(f"Invalid level choice: {level}. Please choose from loci, subfamily, family, or class.")
            return
    elif mode == "TElocal":
        if level == "loci":
            extract_te_index = 3
        elif level == "subfamily":
            extract_te_index = 0
        elif level == "family":
            extract_te_index = 1
        elif level == "class":
            extract_te_index = 2
        else:
            logger.error(f"Invalid level choice: {level}. Please choose from loci, subfamily, family, or class.")
            return
    logger.info("Data loaded successfully.")
    if genome1_samples:
        df_plot1 = process_group(df2, genome1_samples, min_count, top_n, extract_te_index)
        outfile1 = os.path.join(outdir, f"{genome1}_2_{genome2}.data.csv")
        df_plot1.to_csv(outfile1, index=False)
        outImage1 = os.path.join(outdir, f"{genome1}_2_{genome2}.png")
        Double_bar(df_plot1, outImage1)

    if genome2_samples:
        df_plot2 = process_group(df1, genome2_samples, min_count, top_n, extract_te_index)
        outfile2 = os.path.join(outdir, f"{genome2}_2_{genome1}.data.csv")
        df_plot2.to_csv(outfile2, index=False)
        outImage2 = os.path.join(outdir, f"{genome2}_2_{genome1}.png")
        Double_bar(df_plot2, outImage2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot gene and TE counts")
    parser.add_argument("-i1","--infile1", required=True,help="TEtranscripts output file, genome1")
    parser.add_argument("-i2","--infile2", required=True,help="TEtranscripts output file, genome2")
    parser.add_argument("-g1","--genome1_samples", nargs="+", help="genome1 sample names")
    parser.add_argument("-g2","--genome2_samples", nargs="+", help="genome2 sample names")
    parser.add_argument("-o","--outdir", required=True,help="Output directory")
    parser.add_argument("-mc","--min_count", type=float, default=5)
    parser.add_argument("--top", type=int, default=10, help="Number of top features to display")
    parser.add_argument("-m","--mode", choices=["TEcount","TElocal"], default="TEcount", help="Type of TEtranscripts output file")
    parser.add_argument("-l", "--level", choices=["loci","subfamily","family","class"], default="subfamily", help="TE classification level to extract for plotting")
    args = parser.parse_args()

    prepare_data(args.infile1, args.infile2, args.outdir, args.genome1_samples, args.genome2_samples, args.min_count, args.top)