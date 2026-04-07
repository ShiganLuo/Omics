import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from random import sample
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.colors import to_rgb
import numpy as np

def plot_heatmap_with_group_colors(
        heatmap_data,
        outplot,
        row_rename_map: dict,
        column_rename_map: dict,
        show_group_labels: bool = True,
        cmap: str = "coolwarm",

        # unified GridSpec weights
        y_colorbar_width: float = 1,
        heatmap_width: float = 50,
        right_cbar_width: float = 2,

        heatmap_height: float = 50,
        x_colorbar_height: float = 0.6,

        # ---- NEW ----
        min_group_fraction_for_label: float = 0.03   # 小组不显示文字
):
    """
    Only show text for large groups.
    Small groups appear as legend on the far right.
    """

    # ---- 1. rename ----
    heatmap_data = heatmap_data.copy()
    heatmap_data.rename(index=row_rename_map, inplace=True)
    heatmap_data.rename(columns=column_rename_map, inplace=True)

    col_order = [c for c in column_rename_map.values() if c in heatmap_data.columns]
    heatmap_data = heatmap_data.loc[:, col_order]

    # ---- 2. group list ----
    row_groups = list(heatmap_data.index)
    n_rows = len(row_groups)

    # ---- 3. 统计每个组行数 ----
    group_sizes = {}
    for g in row_groups:
        group_sizes[g] = group_sizes.get(g, 0) + 1

    total_n = float(n_rows)

    # ---- 4. 按 group 大小降序排序 ----
    sorted_groups = sorted(group_sizes.keys(), key=lambda g: group_sizes[g], reverse=True)

    # reorder rows
    ordered_index = []
    for g in sorted_groups:
        ordered_index.extend([idx for idx in heatmap_data.index if idx == g])

    heatmap_data = heatmap_data.loc[ordered_index, :]

    # update row groups
    row_groups = list(heatmap_data.index)
    n_rows, n_cols = heatmap_data.shape

    # ---- 5. define colors ----
    color_pool = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
    ]

    unique_row_groups = sorted_groups
    row_color_map = {g: color_pool[i % len(color_pool)] for i, g in enumerate(unique_row_groups)}

    # col group colors
    uniq_col = list(dict.fromkeys(list(heatmap_data.columns)))
    palette = sns.color_palette("hls", len(uniq_col))
    col_color_map = {g: palette[i] for i, g in enumerate(uniq_col)}

    # ---- 6. row group blocks ----
    blocks = []
    start = 0
    for i in range(1, n_rows):
        if row_groups[i] != row_groups[i - 1]:
            blocks.append((row_groups[i - 1], start, i))
            start = i
    blocks.append((row_groups[-1], start, n_rows))

    # ---- 6.1 识别大组、小组 ----
    big_groups = [g for g in sorted_groups if group_sizes[g] / total_n >= min_group_fraction_for_label]
    small_groups = [g for g in sorted_groups if g not in big_groups]

    # ---- 7. GridSpec ----
    fig = plt.figure(figsize=(8, 16))
    gs = gridspec.GridSpec(
        2, 3,
        width_ratios=[y_colorbar_width, heatmap_width, right_cbar_width],
        height_ratios=[heatmap_height, x_colorbar_height],
        wspace=0,
        hspace=0
    )

    ax_ybar = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_cbar = fig.add_subplot(gs[0, 2])
    ax_blank = fig.add_subplot(gs[1, 0])
    ax_xbar = fig.add_subplot(gs[1, 1])
    ax_blank_r = fig.add_subplot(gs[1, 2])

    ax_blank.axis("off")
    ax_blank_r.axis("off")

    vmin = np.percentile(heatmap_data.values, 5)
    vmax = np.percentile(heatmap_data.values, 95)
    # ---- 8. heatmap ----
    im = ax_heat.imshow(
        heatmap_data.values, 
        aspect='auto', 
        origin='upper', 
        cmap=cmap,
        vmin = vmin,
        vmax=vmax)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])

    # ---- 9. colorbar ----
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.ax.text(0.5, 1.02, "log2(cpm+1)", ha='center', va='bottom', transform=cbar.ax.transAxes)

    # ---- 10. y-axis color strip ----
    ax_ybar.set_xlim(0, 1)
    ax_ybar.set_ylim(0, n_rows)
    ax_ybar.invert_yaxis()

    for g, s, e in blocks:
        ax_ybar.add_patch(
            patches.Rectangle((0, s), 1, e - s, facecolor=row_color_map[g], edgecolor='none')
        )

    # ---- 大 group：左侧显示文字 ----
    if show_group_labels:
        for g, s, e in blocks:
            if g in big_groups:     # only large groups
                ax_ybar.text(
                    -0.1, (s + e) / 2, g,
                    ha="right", va="center", fontsize=8
                )

    ax_ybar.axis("off")

    # ---- 11. x-axis color strip ----
    col_groups = list(heatmap_data.columns)
    ax_xbar.set_xlim(0, n_cols)
    ax_xbar.set_ylim(0, 1)

    for i, g in enumerate(col_groups):
        ax_xbar.add_patch(
            patches.Rectangle((i, 0), 1, 1, facecolor=col_color_map[g], edgecolor='none')
        )

    if show_group_labels:
        for g in uniq_col:
            pos = [i for i, x in enumerate(col_groups) if x == g]
            mid = (pos[0] + pos[-1] + 1) / 2
            ax_xbar.text(mid, 0, g, ha="center", va="top", rotation=90, fontsize=8)

    ax_xbar.axis("off")

    # ---- 12. 小 group legend（最右边） ----
    if len(small_groups) > 0:
        handles = []
        labels = []
        for g in small_groups:
            handles.append(plt.Rectangle((0, 0), 1, 1, color=row_color_map[g]))
            labels.append(g)

        ax_cbar.legend(
            handles, labels,
            title="Groups",
            loc="upper left",
            bbox_to_anchor=(5, 1),
            fontsize=8,
            title_fontsize=9
        )

    # ---- SAVE ----
    plt.savefig(outplot, dpi=300, bbox_inches='tight')
    plt.close(fig)


def TE_heatmap(
        TEcount:str,
        sample_group:dict,
        outpng:str
):
    """
    sample_group: dict, {sample:group,……}
    """
    df = pd.read_csv(TEcount,sep="\t",index_col=0)
    pattern = r"^([^:]+):([^:]+):[^:]+$"
    df = df.loc[df.index.str.match(pattern),:]
    subfamilies = df.index.str.extract(pattern)
    subfamilies.columns = ["subfamily", "family"]
    subfamily_to_family = dict(zip(subfamilies["subfamily"], subfamilies["family"]))
    df.index = subfamilies["subfamily"].values
    print(df)
    cpm = df.div(df.sum(axis=0), axis=1) * 1e6
    log2_cpm = np.log2(cpm + 1)
    plot_heatmap_with_group_colors(log2_cpm,outpng,subfamily_to_family,sample_group)




if __name__ == "__main__":
    sample_group = {
        "GSM4777760": "2iL",
        "GSM4777761": "2iL",
        "GSM4777766": "2iL-F",
        "GSM4777767": "2iL-F",
        "GSM4777762": "A2iL",
        "GSM4777763": "A2iL",
        "GSM4777764": "LCDM",
        "GSM4777765": "LCDM",
        "GSM4777758": "SL",
        "GSM4777759": "SL"
    }
    # TE_heatmap("/disk5/luosg/Totipotent20251031/PRJNA663159/TEcount/mouse/allTEcount.cntTable",
    #            sample_group)
    sample_group = {
        "GSM5924234":"prEpiSC",
        "GSM5924235":"prEpiSC",
        "GSM5924236":"ci8CLC",
        "GSM5924237":"ci8CLC",
        "GSM5924238":"ci8CLC",
        "GSM5924239":"ci8CLC",
        "GSM7032352":"hESC",
        "GSM7719025":"hESC",
        "GSM7719026":"hESC",
        "GSM7032363":"hTBLC",
        "GSM7032364":"hTBLC",
        "GSM7032365":"hTBLC",
        "GSM7032366":"hTBLC",
        "GSM7032367":"hTBLC"
    }
    TE_heatmap("/disk5/luosg/Totipotent20251031/output/counts/TEcount/human/all_TEcount.cntTable",
               sample_group,
               "heatmap.png")

