import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
def read_gmt(file):
    gene_sets = {}
    with open(file) as f:
        for line in f:
            parts = line.strip().split("\t")
            gene_set_name = parts[0]
            genes = parts[2:]  # GMT第3列开始是基因
            gene_sets[gene_set_name] = genes
    return gene_sets


def plot_heatmap_with_colgroup(
        data: pd.DataFrame,
        col_group_map: dict,
        row_labels: dict = None,
        outplot: str = "heatmap.png",
        figsize:tuple = (6,6)
    ):
    """
    不聚类热图，x轴列名染色显示分组，列分组图例在颜色条下方，每个分组独占一行
    """
    if row_labels:
        data = data.rename(index=row_labels)

    cmap = sns.color_palette("RdBu_r", as_cmap=True)
    vmin, vmax = np.nanmin(data.values), np.nanmax(data.values)

    fig, ax = plt.subplots(figsize=figsize)

    # 右侧颜色条
    cbar_ax = fig.add_axes([0.88, 0.25, 0.03, 0.7])
    sns.heatmap(
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        ax=ax,
        cbar_ax=cbar_ax,
        linewidths=0.5,
        linecolor='gray'
    )

    # 列分组颜色
    col_groups = pd.Series(col_group_map)
    unique_groups = col_groups.unique()
    palette = sns.color_palette("Set2", n_colors=len(unique_groups))
    group_colors = dict(zip(unique_groups, palette))

    # x轴标签染色
    for xtick_label in ax.get_xticklabels():
        label = xtick_label.get_text()
        if label in col_groups:
            xtick_label.set_color(group_colors[col_groups[label]])
            xtick_label.set_rotation(45)
            xtick_label.set_horizontalalignment('right')

    # 调整子图区域，留出下边距和右边距
    fig.subplots_adjust(bottom=0.2, right=0.78, top=0.95)

    # 绘制分组图例，位于右侧颜色条下方，每个分组占一行
    legend_handles = [plt.Line2D([0], [0], color=color, lw=6) for color in group_colors.values()]
    fig.legend(
        handles=legend_handles,
        labels=list(group_colors.keys()),
        loc='center right',
        bbox_to_anchor=(1.0, 0.15),  # 下方偏移
        ncol=1,
        frameon=False,
        title=""
    )

    plt.savefig(outplot, dpi=300)
    plt.close()

    
if __name__ == "__main__":
    gmt_file = "/disk5/luosg/Totipotent20251031/data/geneset/Totipotency.gmt"
    gene_sets = read_gmt(gmt_file)
    tpm_human = pd.read_csv("/disk5/luosg/Totipotent20251031/output/counts/featureCounts/Homo_sapiens/Homo_sapiens_paired_name_tpm.tsv", sep = "\t",index_col=0)
    pre_ZGA = gene_sets["pre-ZGA(Homo)"]
    sample_list = ["GSM7032352","GSM7719025","GSM7719026","GSM7032363","GSM7032364","GSM7032365","GSM7032366","GSM7032367"]
    tpm_hTBLC = tpm_human.loc[pre_ZGA,sample_list]
    tpm_hTBLC_z = tpm_hTBLC.apply(
        lambda x: (x - x.mean()) / x.std(ddof=0), axis=1
    )
    print(tpm_hTBLC)
    col_rename_map = {
        "GSM7032352":"hESC",
        "GSM7719025":"hESC",
        "GSM7719026":"hESC",
        "GSM7032363":"hTBLC",
        "GSM7032364":"hTBLC",
        "GSM7032365":"hTBLC",
        "GSM7032366":"hTBLC",
        "GSM7032367":"hTBLC"
    }
    plot_heatmap_with_colgroup(tpm_hTBLC_z , 
                               col_rename_map, 
                               outplot="/disk5/luosg/Totipotent20251031/output/result/hTBLC/totipotency/totipotency_heatmap.png",
                               figsize=(6.14,6))