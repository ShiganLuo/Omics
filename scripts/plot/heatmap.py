from os import path
from random import sample
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from random import sample
def plot_heatmap_with_group_colors(
        heatmap_data:pd.DataFrame,
        outplot:str,
        row_rename_map: dict,
        column_rename_map: dict,
        show_group_labels:bool = True,   # 是否显示分组文字
        title:str = "Heatmap of Positions",
        xlabel:str = "Positions",
        colorbar_width:float = 0.02,
        x_annotation_height:float = 0.2
):
    ### data preprocessing
    heatmap_data.rename(index=row_rename_map, inplace=True)
    heatmap_data.rename(columns=column_rename_map, inplace=True)
    heatmap_data = heatmap_data.sort_index()
    new_col_order = [new_name for new_name in column_rename_map.values()]
    heatmap_data = heatmap_data[new_col_order]
    group_colors = {}
    for group in heatmap_data.index:
        if group not in group_colors:
            group_colors[group] = sample(
                [
                    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                    '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
                ],
                1
            )[0]
    
    n_rows = heatmap_data.shape[0]
    fig = plt.figure(figsize=(7, 11))
    gs = gridspec.GridSpec(2, 2, width_ratios=[colorbar_width, 1],height_ratios = [1,x_annotation_height],wspace=0,hspace=0.3)
    # ---------------------
    # 颜色条
    ax_bar = fig.add_subplot(gs[0, 0])
    group_positions = {}
    for idx, group in enumerate(heatmap_data.index):
        ax_bar.add_patch(
            plt.Rectangle((0.6, idx), 0.3, 1, facecolor=group_colors[group])
        )
        if group not in group_positions:
            group_positions[group] = []
        group_positions[group].append(idx)
    # 添加分组标签
    if show_group_labels:
        for group, positions in group_positions.items():
            mid_idx = positions[len(positions)//2]
            ax_bar.text(
                0, mid_idx + 0.5, group,  # 文字在颜色条左边
                ha='right', va='center',
                fontsize=8
            )
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, n_rows)
    ax_bar.invert_yaxis()
    ax_bar.axis('off')
    # ---------------------

    # ---------------------
    # 热图
    ax_heat = fig.add_subplot(gs[0, 1])
    sns.heatmap(
        heatmap_data,
        ax=ax_heat,
        cmap='coolwarm',
        linewidths=0.01,
    )
    for label in ax_heat.get_xticklabels():
        label.set_rotation(60)
        label.set_fontsize(8)   # ★ 设置大小
        label.set_ha('center')
    ax_heat.set_yticks([])  # 去掉 y 轴
    ax_heat.set_ylabel('')
    ax_heat.set_xlabel(xlabel)
    ax_heat.set_title(title)
    cbar = ax_heat.collections[0].colorbar
    cbar.ax.text(
        0.5, 1.02,  # x=0.5 colorbar 中心, y=1.05 颜色条上方
        'ES',
        ha='center', va='bottom',
        rotation=0, fontsize=8,
        transform=cbar.ax.transAxes
    )
    # ---------------------

    # ---------------------
    # 在热图下方添加 x 轴分组图例
    ax_legend = fig.add_subplot(gs[1, :]) # 下方图例，与热图同列
    ax_legend.axis('off')
    short_names = list(heatmap_data.columns)
    pathway_rename_map = {v:k for k,v in column_rename_map.items()}
    legend_texts = [
        f"{abbr}: {pathway_rename_map[abbr]}"
        for abbr in short_names
        if abbr in pathway_rename_map
    ]

    max_chars_per_row = 60  # 每行最大字符数
    rows = []
    current_row = []
    current_len = 0

    for text in legend_texts:
        text_len = len(text)
        if current_len + text_len > max_chars_per_row and current_row:
            # 当前行满了，换行
            rows.append(current_row)
            current_row = [text]
            current_len = text_len
        else:
            current_row.append(text)
            current_len += text_len + 2  # 预留间隔
    if current_row:
        rows.append(current_row)
    n_rows = len(rows)
    # 绘制文字，左对齐
    for row_idx, row in enumerate(rows):
        x_pos = 0  # 每行从左边开始
        y_pos = n_rows - row_idx
        for text in row:
            ax_legend.text(
                x_pos, y_pos,
                text,
                ha='left', va='center',  # 左对齐
                rotation=0,
                fontsize=8
            )
            x_pos += len(text) + 2  # 累加下一个文字起始位置

    # 设置坐标范围，让文字显示完整
    ax_legend.set_xlim(0, max_chars_per_row)
    ax_legend.set_ylim(0, len(rows))
    # ---------------------
    plt.savefig(outplot, dpi=300)
    plt.close()

def plot_heatmap_with_group_colors_auto_orient(
        heatmap_data: pd.DataFrame,
        outplot: str,
        row_rename_map: dict,
        column_rename_map: dict,
        title: str = "Heatmap",
        xlabel: str = "",
        force_horizontal: bool = False,
        force_vertical: bool = False,
):
    """
    绘制带有“分组颜色条 + 列注释”的热图，并根据数据形状自动选择方向。

    设计原则
    --------
    1. 分组语义始终锚定在“原始行”（未转置前）
    2. rotate 只改变可视方向，不改变分组含义
    3. 长轴用于连续分组颜色条，短轴用于文字注释
    4. heatmap 的数值 colorbar 永远是垂直的
    数据形状：(样本*含义)
    """

    if force_horizontal and force_vertical:
        raise ValueError("force_horizontal and force_vertical cannot both be True")

    # =========================
    # 数据准备（语义锚点）
    # =========================
    data = heatmap_data.copy()
    data.rename(index=row_rename_map, inplace=True)
    data.rename(columns=column_rename_map, inplace=True)
    data = data.sort_index()
    data = data[list(column_rename_map.values())]

    original_rows = list(data.index)
    n_rows, n_cols = data.shape

    # 行分组（始终基于原始行）
    group_blocks = {}
    for i, g in enumerate(original_rows):
        group_blocks.setdefault(g, []).append(i)

    palette = [
        "#e6194b", "#3cb44b", "#ffe119", "#4363d8",
        "#f58231", "#911eb4", "#46f0f0", "#f032e6"
    ]
    group_colors = {
        g: palette[i % len(palette)]
        for i, g in enumerate(group_blocks)
    }

    # =========================
    # 布局方向
    # =========================
    if force_horizontal:
        rotate = False
    elif force_vertical:
        rotate = True
    else:
        rotate = n_rows > n_cols

    plot_data = data.T if rotate else data

    # =========================
    # Figure / GridSpec
    # =========================

    if not rotate:
        cell = 0.35
        fig_w = max(6, plot_data.shape[1] * cell * 3)
        fig_h = max(4, plot_data.shape[0] * cell * 0.3)

        fig = plt.figure(figsize=(fig_w + 0.6, fig_h))
        gs = gridspec.GridSpec(
            2, 3,
            height_ratios=[1, 0.3],
            width_ratios=[0.05, 1, 0.08],
            hspace=0.0,
            wspace=0.0
        )
        ax_group_long = fig.add_subplot(gs[0, 0])
        ax_group_short = fig.add_subplot(gs[1,1])
        ax_heat  = fig.add_subplot(gs[0, 1])
        ax_cbar  = fig.add_subplot(gs[0, 2])
    else:
        cell = 0.35
        fig_w = max(6, plot_data.shape[1] * cell * 0.2)
        fig_h = max(4, plot_data.shape[0] * cell * 2)

        fig = plt.figure(figsize=(fig_w + 0.6, fig_h))
        gs = gridspec.GridSpec(
            3, 2,
            height_ratios=[0.05, 1, 0.2],
            width_ratios=[1, 0.02],
            hspace=0.0,
            wspace=0.02
        )
        ax_group_long = fig.add_subplot(gs[0, 0])
        ax_group_short = fig.add_subplot(gs[2, 0])
        ax_heat  = fig.add_subplot(gs[1, 0])
        ax_cbar  = fig.add_subplot(gs[1, 1])

    # =========================
    # Heatmap
    # =========================
    sns.heatmap(
        plot_data,
        ax=ax_heat,
        cmap="coolwarm",
        linewidths=0.01,
        cbar=False
    )

    ax_heat.set_title(title)
    ax_heat.set_ylabel(None)

    if not rotate:
        # x轴刻度保留
        ax_heat.set_xlabel(xlabel)
        ax_heat.set_yticks([])

    else:
        # y轴刻度保留
        ax_heat.set_xlabel("")
        ax_heat.set_xticks([])
        ax_heat.tick_params(axis="y", labelsize=12,labelrotation=0)

    # =========================
    # 分组颜色条或注释（刻度注释层）
    # =========================
    if not rotate:
        # 长轴颜色条
        draw_group_color_bar(
            ax=ax_group_long,
            group_blocks=group_blocks,
            group_colors=group_colors,
            n_units=n_rows,
            orientation="y")
        # 短轴注释
        draw_xaxis_text_legend(
            ax=ax_group_short,
            data=data,
            column_rename_map=column_rename_map,
            max_chars_per_row=30
            )
    else:
        # 长轴颜色条
        draw_group_color_bar(
            ax=ax_group_long,
            group_blocks=group_blocks,
            group_colors=group_colors,
            n_units=n_rows,
            orientation="x")

        # 短轴注释
        draw_xaxis_text_legend(
            ax=ax_group_short,
            data=data,
            column_rename_map=column_rename_map
            )
    

    # =========================
    # heamtp colorbar
    # =========================
    plt.colorbar(
        ax_heat.collections[0],
        cax=ax_cbar,
        orientation="vertical"
    )
    ax_cbar.tick_params(labelsize=12)

    plt.savefig(outplot, dpi=300, bbox_inches="tight")
    plt.close()

def draw_group_color_bar(
    ax,
    group_blocks: dict,
    group_colors: dict,
    n_units: int,
    orientation: str = "x",
    show_group_labels: bool = True,
    fontsize: int = 12,
    label_offset: float = 1.05,
):
    """
    绘制分组颜色条（支持 x / y 轴方向）

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        用于绘制分组颜色条的轴
    group_blocks : dict
        {group_name: [indices]}，每个 group 对应连续区间
    group_colors : dict
        {group_name: color}
    n_units : int
        在主轴方向上的单元数量（x 方向用列数，y 方向用行数）
    orientation : {"x", "y"}
        颜色条沿哪个轴方向绘制
    show_group_labels : bool
        是否绘制分组名称
    fontsize : int
        分组标签字体大小
    label_offset : float
        标签偏移量（x 模式是 y 偏移，y 模式是 x 偏移）
    """

    if orientation not in {"x", "y"}:
        raise ValueError("orientation must be 'x' or 'y'")

    ax.axis("off")

    # =========================
    # x 轴方向（原始行为）
    # =========================
    if orientation == "x":
        ax.set_xlim(0, n_units)
        ax.set_ylim(0, 1)

        for g, idxs in group_blocks.items():
            start, end = min(idxs), max(idxs) + 1

            ax.add_patch(
                plt.Rectangle(
                    (start, 0),
                    end - start,
                    1,
                    facecolor=group_colors[g],
                )
            )

            if show_group_labels:
                ax.text(
                    (start + end) / 2,
                    label_offset,
                    g,
                    ha="center",
                    va="bottom",
                    fontsize=fontsize,
                )

    # =========================
    # y 轴方向（新支持）
    # =========================
    else:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, n_units)
        ax.invert_yaxis()  # 与 heatmap 的 y 轴方向保持一致

        for g, idxs in group_blocks.items():
            start, end = min(idxs), max(idxs) + 1

            ax.add_patch(
                plt.Rectangle(
                    (0, start),
                    1,
                    end - start,
                    facecolor=group_colors[g],
                )
            )

            if show_group_labels:
                ax.text(
                    -label_offset,
                    (start + end) / 2,
                    g,
                    ha="right",
                    va="center",
                    fontsize=fontsize,
                )


def draw_xaxis_text_legend(
    ax,
    data,
    column_rename_map: dict,
    max_chars_per_row: int = 130,
    fontsize: int = 12,
):
    """
    在给定 axis 上绘制 x 轴缩写 → 全名 的多行文字注释（左对齐）

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        用于绘制文字注释的轴（通常位于 heatmap 下方）
    data : pd.DataFrame
        热图使用的数据（用于获取 columns 顺序）
    column_rename_map : dict
        原始名 -> 缩写名 的映射
    max_chars_per_row : int
        每一行允许的最大字符数
    fontsize : int
        字体大小
    """

    ax.axis("off")

    # 缩写列表（按 heatmap 当前列顺序）
    short_names = list(data.columns)

    # 反向映射：缩写 -> 全名
    pathway_rename_map = {v: k for k, v in column_rename_map.items()}

    # 构造 legend 文本
    legend_texts = [
        f"{abbr}: {pathway_rename_map[abbr]}"
        for abbr in short_names
        if abbr in pathway_rename_map
    ]

    # =========================
    # 自动分行（按字符数）
    # =========================
    rows = []
    current_row = []
    current_len = 0

    for text in legend_texts:
        text_len = len(text)
        if current_len + text_len > max_chars_per_row and current_row:
            rows.append(current_row)
            current_row = [text]
            current_len = text_len
        else:
            current_row.append(text)
            current_len += text_len + 2  # 预留间隔

    if current_row:
        rows.append(current_row)

    # =========================
    # 绘制文字（左对齐）
    # =========================
    n_rows = len(rows) - 1

    for row_idx, row in enumerate(rows):
        x_pos = 0
        y_pos = n_rows - row_idx

        for text in row:
            ax.text(
                x_pos,
                y_pos,
                text,
                ha="left",
                va="center",
                rotation=0,
                fontsize=fontsize,
            )
            x_pos += len(text) + 2

    # 坐标范围，保证完整显示
    ax.set_xlim(0, max_chars_per_row)
    ax.set_ylim(0, len(rows))

def imshow_heatmap(
        df:pd.DataFrame
):
    fig, ax = plt.subplots()

if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/RNAseqML/gsva230/gseapy.gene_set.gsva.report.csv")
    # print(df.head())
    df_group = pd.read_csv("/disk5/luosg/Totipotent20251031/data/RNAseqML/target_fq.tsv",sep="\t")
    sample_rename_map = dict(zip(df_group['Sample_id'],df_group['Status']))
    pathway_rename_map = {
        "DNA Double Strand Break Response": "DSBR",
        "HDR through Homologous Recombination (HRR)": "HRR",
        "Nonhomologous End-Joining (NHEJ)": "NHEJ",
        "Mismatch Repair": "MMR",
        "Nucleotide Excision Repair": "NER",
        "Base Excision Repair": "BER",
        "Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template":"TLS",
        "Fanconi Anemia Pathway": "FA",
        "Activation of ATR in response to replication stress": "RS",
        "Telomere Maintenance": "TM",
        "Mitotic Spindle Checkpoint": "MSC",
        "Cell Cycle Checkpoints": "CCC"
    }
    df_wide = df.pivot(index='Name', columns='Term', values='ES')
    # df_wide.to_csv("a.csv")
    row_colors = {'Tumor':'red','Normal_cell':'blue',"Aged_cell":'green',"Stem_cell":"black"}
    # plot_heatmap_with_group_colors(df_wide,
    #              outplot = "/home/luosg/Data/genomeStability/output/result/gsva230/heatmap230.png",
    #              row_rename_map = sample_rename_map,
    #              column_rename_map = pathway_rename_map,
    #              title="",
    #              xlabel="Pathways")
    print(df_wide.head())
    plot_heatmap_with_group_colors_auto_orient(df_wide,
                outplot = "/disk5/luosg/Totipotent20251031/RNAseqML/gsva230/heatmap230.png",
                row_rename_map = sample_rename_map,
                column_rename_map = pathway_rename_map,
                title="",
                xlabel="Pathways",
                force_horizontal=False)
