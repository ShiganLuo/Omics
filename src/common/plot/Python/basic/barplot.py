import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import textwrap
from pathlib import Path
from typing import List, Union, Dict
from scipy.stats import chi2_contingency

def plot_group_type_comparison(
    df: pd.DataFrame,
    out_png: str,
    group_col: str = "group",
    type_col: str = "svtype",
    count_col: str = "count",
    group_order: tuple = ("Control", "Experiment"),
    type_order: tuple = ("BND", "DEL", "DUP", "INS", "INV"),
    legend_map: Optional[Dict[str, str]] = None,
    figsize: tuple = (9, 6),
    ylabel: str = "SV count",
    dpi: int = 300,
    use_broken_axis: bool = True,
    do_test: bool = True,
    test_method: str = "chi2",
    colors: Optional[Dict[str, str]] = None,
    bracket_gap: float = 20.0,
    tick_depth: float = 6.0,
    bracket_lw: float = 0.8,
) -> None:
    """
    Plot a grouped bar chart comparing SV type counts across multiple groups.

    Performs chi-squared tests per SV type across all groups and optionally
    annotates pairwise significance on a broken y-axis.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format DataFrame with group, SV type, and count columns.
    out_png : str
        Output file path for the saved image.
    group_col : str, optional
        Column identifying sample groups. Default is ``"group"``.
    type_col : str, optional
        Column identifying SV types. Default is ``"svtype"``.
    count_col : str, optional
        Column with counts. Default is ``"count"``.
    group_order : tuple of str, optional
        Ordered group names to display. Any number of groups supported.
    type_order : tuple of str, optional
        Ordered SV types to display.
    legend_map : dict, optional
        Mapping from group name to legend label. Defaults to identity.
    figsize : tuple, optional
        Figure size ``(width, height)``. Default is ``(9, 5)``.
    ylabel : str, optional
        Y-axis label. Default is ``"SV count"``.
    dpi : int, optional
        Resolution in dots per inch. Default is ``300``.
    use_broken_axis : bool, optional
        If ``True``, use a split y-axis to accommodate large count
        differences. Default is ``True``.
    do_test : bool, optional
        If ``True``, perform chi-squared tests and annotate significance.
        Default is ``True``.
    test_method : str, optional
        Statistical test method: ``"chi2"`` or ``"fisher"``.
        Default is ``"chi2"``.
    colors : dict, optional
        Mapping from group name to color. Auto-generated from ``tab10``
        if ``None``.
    bracket_gap : float, optional
        Vertical gap in display points between stacked bracket units
        within the same SV type. Default is ``10.0``.
    tick_depth : float, optional
        Depth of the outward diagonal ticks in display points. Larger
        values produce steeper ticks. Default is ``6.0``.
    bracket_lw : float, optional
        Line width for bracket lines and ticks. Default is ``0.8``.

    Returns
    -------
    None
        Saves the figure to *out_png* at the specified DPI.
    """
    outdir = os.path.dirname(out_png)
    os.makedirs(outdir, exist_ok=True)

    n_groups = len(group_order)

    if legend_map is None:
        legend_map = {g: g for g in group_order}

    pivot = (
        df.pivot(index=type_col, columns=group_col, values=count_col)
        .reindex(type_order)
        .fillna(0)
    )

    def p_to_star(p: float) -> str:
        if p < 1e-4:
            return "****"
        elif p < 1e-3:
            return "***"
        elif p < 1e-2:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    # ---------- significance: pairwise chi2 ----------
    all_pairs: Dict[str, List[tuple]] = {}
    if do_test:
        for sv in pivot.index:
            pairs = []
            for i, g1 in enumerate(group_order):
                for g2 in group_order[i + 1:]:
                    idx1 = list(group_order).index(g1)
                    idx2 = list(group_order).index(g2)
                    table = np.array([
                        [pivot.loc[sv, g1], pivot[g1].sum() - pivot.loc[sv, g1]],
                        [pivot.loc[sv, g2], pivot[g2].sum() - pivot.loc[sv, g2]],
                    ])
                    # Skip degenerate tables
                    if table.min() < 0 or table.sum() == 0 or (table.sum(axis=0) == 0).any() or (table.sum(axis=1) == 0).any():
                        continue
                    try:
                        if test_method == "fisher":
                            _, p = fisher_exact(table)
                        else:
                            _, p, _, _ = chi2_contingency(table)
                    except ValueError:
                        continue
                    pairs.append((idx1, idx2, p_to_star(p)))
            all_pairs[sv] = pairs

    # ---------- colors ----------
    if colors is None:
        cmap = plt.cm.tab10
        colors = {g: cmap(i / max(n_groups - 1, 1)) for i, g in enumerate(group_order)}

    # ---------- axis setup ----------
    x = np.arange(len(pivot.index))
    total_width = 0.8
    bar_width = total_width / n_groups
    offsets = [bar_width * (i - (n_groups - 1) / 2) for i in range(n_groups)]

    # Auto-detect whether broken axis is needed
    global_max = pivot.values.max()
    sig_sv_list = [sv for sv, pairs in all_pairs.items()
                   if any(s != "ns" for _, _, s in pairs)]
    if sig_sv_list:
        sig_max = pivot.loc[sig_sv_list].values.max()
    else:
        sig_max = np.median(pivot.values)

    # If the max bar is not significantly taller than the rest, skip broken axis
    need_broken = use_broken_axis and (global_max > sig_max * 2.0)

    if need_broken:
        low_max = sig_max * 1.15
        high_min = low_max * 1.1
        high_max = global_max * 1.15

        fig, (ax_top, ax_bottom) = plt.subplots(
            2, 1, sharex=True,
            figsize=figsize,
            gridspec_kw={"height_ratios": [1, 3]},
        )

        axes = (ax_top, ax_bottom)

        for ax in axes:
            for i, g in enumerate(group_order):
                ax.bar(x + offsets[i], pivot[g], bar_width,
                       color=colors[g], label=legend_map[g])

        ax_bottom.set_ylim(0, low_max)
        ax_top.set_ylim(high_min, high_max)

        # broken marks
        d = 0.008
        ax_top.plot((-d, +d), (-d, +d), transform=ax_top.transAxes,
                     color="black", clip_on=False)
        ax_bottom.plot((-d, +d), (1 - d, 1 + d), transform=ax_bottom.transAxes,
                        color="black", clip_on=False)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        for i, g in enumerate(group_order):
            ax.bar(x + offsets[i], pivot[g], bar_width,
                   color=colors[g], label=legend_map[g])

        axes = (ax,)
        ax_top = ax_bottom = ax

    # ---------- significance brackets (宝盖头 style) ----------
    if do_test:
        LEG_PT = 8
        TEXT_PT = 3

        fig.canvas.draw()

        for i_sv, sv in enumerate(pivot.index):
            pairs = all_pairs.get(sv, [])
            if not pairs:
                continue

            y_base = max(pivot.loc[sv, g] for g in group_order)
            bracket_offset = 0.0

            for idx1, idx2, star in pairs:
                if need_broken:
                    ax = ax_bottom if y_base <= low_max else ax_top
                else:
                    ax = ax_top

                trans = ax.transData
                inv = ax.transData.inverted()

                _, y_disp = trans.transform((0, y_base))
                y_hat_disp = y_disp + LEG_PT + bracket_offset
                y_text_disp = y_hat_disp + TEXT_PT
                _, y_hat = inv.transform((0, y_hat_disp))
                _, y_text = inv.transform((0, y_text_disp))

                x1 = x[i_sv] + offsets[idx1]
                x2 = x[i_sv] + offsets[idx2]

                is_sig = star != "ns"
                fs = 12 if is_sig else 9
                fw = "bold" if is_sig else "normal"

                # 宝盖头 bracket as single polyline to avoid line-cap overlap
                tick_x = bar_width * 0.25
                _, y_tick = inv.transform((0, y_hat_disp - tick_depth))

                ax.plot(
                    [x1 - tick_x, x1, x2, x2 + tick_x],
                    [y_tick, y_hat, y_hat, y_tick],
                    lw=bracket_lw, c="black",
                )

                ax.text((x1 + x2) / 2, y_text, star,
                        ha="center", va="bottom",
                        fontsize=fs, fontweight=fw, color="black")

                bracket_offset += tick_depth + TEXT_PT + bracket_gap

    # ---------- formatting ----------
    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(pivot.index)
    ax_bottom.set_ylabel(ylabel)

    if need_broken:
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
        ax_top.spines["bottom"].set_visible(False)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].legend(frameon=False)

    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()


def TEfamily(
    DEG_file: str,
    rmsk_file: str,
    out: str
):
    df_TE = pd.read_csv(DEG_file, sep="\t")
    df_TE.reset_index(names="subfamily", inplace=True)
    df_rmsk = pd.read_csv(rmsk_file, sep="\t", header=None)
    
    df_annotation = df_rmsk[[10, 11]].drop_duplicates(keep="first")
    df_annotation.columns = ["subfamily", "family"]
    
    new_df = pd.merge(df_TE, df_annotation, on="subfamily", how="left")
    
    # 筛选显著
    significant_df = new_df[
        (abs(new_df["logFC"]) > 0.58) & 
        (-np.log10(new_df["PValue"]) > 1.3)
    ].copy()
    
    if significant_df.empty:
        raise ValueError("The filtered DataFrame of significant TEs is empty.")
    
    # 直接分组计数：同时保留方向信息
    significant_df["regulation"] = np.where(significant_df["logFC"] > 0, "Up", "Down")
    result = (
        significant_df.groupby(["family", "regulation"])
        .size()
        .reset_index(name="count")
    )
    
    # 保存 TSV
    result_pivot = result.pivot(index="family", columns="regulation", values="count").fillna(0).astype(int)
    result_pivot.to_csv(f"{out}_TEfamily.tsv", sep="\t", header=True)

    # --- 绘图 ---
    sns.set_style("white")  # 去掉网格线
    fig, ax = plt.subplots(figsize=(12, 6))

    # 分开取调色板：深色给Up，浅色给Down
    up_families = result[result["regulation"] == "Up"]["family"].unique()
    down_families = result[result["regulation"] == "Down"]["family"].unique()
    up_colors = sns.color_palette("dark", len(up_families))
    down_colors = sns.color_palette("pastel", len(down_families))

    # 上调在左侧
    for i, (fam, cnt) in enumerate(result[result["regulation"] == "Up"][["family", "count"]].values):
        ax.bar(i, cnt, color=up_colors[i])
        ax.text(i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    # 下调在右侧
    offset = len(up_families) + 2  # 空2列作为分隔
    for i, (fam, cnt) in enumerate(result[result["regulation"] == "Down"][["family", "count"]].values):
        ax.bar(offset + i, cnt, color=down_colors[i])
        ax.text(offset + i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    # 设置 x 轴刻度
    xticks = list(range(len(up_families))) + list(range(offset, offset + len(down_families)))
    xlabels = list(up_families) + list(down_families)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right")

    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Top TE Families: Upregulated vs Downregulated", fontsize=16)

    # 重新计算 Y 轴上限，横线位置稍高于最大值

    y_max = max(result["count"])
    # 在上方标注 Up 和 Down
    ax.hlines(y=y_max*1.1, xmin=0, xmax=len(up_families)-1, color="black", linewidth=1.5)
    ax.text((len(up_families)-1)/2, y_max*1.13, "Upregulated", ha="center", va="bottom", fontsize=12)

    ax.hlines(y=y_max*1.1, xmin=offset, xmax=offset+len(down_families)-1, color="black", linewidth=1.5)
    ax.text(offset + (len(down_families)-1)/2, y_max*1.13, "Downregulated", ha="center", va="bottom", fontsize=12)

    y_max = max(result["count"]) * 1.25
    ax.set_ylim(0, y_max)
    plt.tight_layout()
    plt.savefig(f"{out}_TEfamily.png", dpi=300)
    plt.close()

def GSEA_BarPlot_Enhanced(
    gsea_result: pd.DataFrame,
    outplot: str,
    pval_cutoff: float = 0.05,
    sort_by: str = "ES",    # "ES" or "pval"
    top_flag: bool = False,
    top_number: int = 5,
    label_wrap: int = 40,
    cmap_name: str = "viridis"
):
    """
    增强版 GSEA barplot：
    - 中心对称：负 ES 左、正 ES 右
    - 按 ES 或 pval 自动排序
    - 颜色按显著性渐变（-log10(pval)）
    - ggplot 风格
    """

    # ------------------------
    # 1. 过滤和预处理
    # ------------------------
    df = gsea_result[gsea_result["pval"] < pval_cutoff].copy()
    df = df.rename(columns={"score": "ES"}).sort_values("ES")
    if top_flag:
        df = pd.concat([df.head(top_number), df.tail(top_number)])
    df["logP"] = -np.log10(df["pval"].replace(0, 1e-7))

    # 名称换行
    df["pathway_wrapped"] = df["source"].apply(
        lambda s: "\n".join(textwrap.wrap(s, label_wrap))
    )


    # ------------------------
    # 4. 图尺寸
    # ------------------------
    height = max(6, len(df) * 0.45)
    fig, ax = plt.subplots(figsize=(12, height))

    # ------------------------
    # 5. 颜色映射
    # ------------------------
    cmap = sns.color_palette(cmap_name, as_cmap=True)
    norm = plt.Normalize(df["logP"].min(), df["logP"].max())
    df["color"] = df["logP"].apply(lambda x: cmap(norm(x)))

    # ------------------------
    # 6. 中心对称 barplot
    # ------------------------
    ax.barh(
        y=df["pathway_wrapped"],
        width=df["ES"],
        color=df["color"],
        edgecolor="black",
        linewidth=0.5
    )

    # ------------------------
    # 7. 中心线（x=0）
    # ------------------------
    ax.axvline(0, color="black", linewidth=1.2)

    # ------------------------
    # 8. 坐标轴与标签
    # ------------------------
    ax.set_xlabel("Enrichment Score (ES)")
    ax.set_ylabel("Pathway")

    ax.set_title(
        "GSEA Enrichment Barplot\n"
        "Left: Negative ES | Right: Positive ES",
        fontsize=14,
        pad=15
    )

    # ------------------------
    # 9. Colorbar (p-value 显著性)
    # ------------------------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("-log10(P-value)")

    # ------------------------
    # 10. 美化
    # ------------------------
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    ax.grid(False, axis="y")

    plt.tight_layout()
    plt.savefig(outplot, dpi=300, bbox_inches="tight")
    plt.close()

# -----------------------------
# 函数3：堆叠柱状图（支持 x 轴彩色标签 + 样本分组图例）
# -----------------------------
def plot_stacking_bar(
    df_counts: pd.DataFrame,
    xlabels: List[str] = None,
    groups: List[str] = None,
    group_colors: Dict[str, str] = None,
    title: str = "Mutation Distribution (Proportion)",
    xlabel: str = "Sample",
    ylabel: str = "Proportion",
    legend_title_type: str = "Mutation Type",
    legend_title_group: str = "Sample Group",
    save_path: Union[str, Path] = None,
    legend_width: float = 0.25,
    figsize: tuple = (12, 6),
    legend_fontsize: int = 13,
    legend_title_fontsize: int = 16,
    rotation: int = 45,
    colormap: str = "tab20"
):
    """
    绘制突变分布的堆叠柱状图，支持比例转化、样本分组上色以及双图例显示。

    Args:
        df_counts (pd.DataFrame): 输入数据，行名为突变类型，列名为样本名。
        xlabels (List[str], optional): X轴刻度标签。默认为 DataFrame 的列名。
        groups (List[str], optional): 样本对应的分组信息，长度需与样本数一致。
        group_colors (Dict[str, str], optional): 分组名到颜色值的映射字典。
        title (str, optional): 图表标题。
        xlabel (str, optional): X轴标题。
        ylabel (str, optional): Y轴标题。
        legend_title_type (str, optional): 突变类型图例的标题。
        legend_title_group (str, optional): 分组图例的标题。
        save_path (Union[str, Path], optional): 图片保存路径。
        legend_width (float, optional): 右侧预留给图例的宽度比例 (0-1)。
        figsize (tuple, optional): 画布尺寸。
        legend_fontsize (int, optional): 图例字体大小。
        legend_title_fontsize (int, optional): 图例标题及轴标签字体大小。
        rotation (int, optional): X轴刻度标签旋转角度。
        colormap (str, optional): 柱状图使用的颜色映射。
    """
    
    # 1. 比例转化：处理单列或多列情况
    # 确保 sum 不为 0 以免除以 0
    df_prop = df_counts.div(df_counts.sum(axis=0).replace(0, 1), axis=1)

    fig, ax = plt.subplots(figsize=figsize)

    # ======== 右侧预留空间给图例 ========
    fig.subplots_adjust(right=1 - legend_width)

    # 2. 绘制堆叠柱状图 (转置后行为样本，列为类型)
    df_prop.T.plot(
        kind="bar",
        stacked=True,
        colormap=colormap,
        width=0.8,
        ax=ax,
        legend=False
    )

    n_samples = df_counts.shape[1]

    # -------------------------------
    # X 轴标签设置
    # -------------------------------
    if xlabels is None:
        xlabels = df_counts.columns.tolist()

    if len(xlabels) != n_samples:
        raise ValueError("xlabels 长度必须与样本数量 (columns) 一致")

    # 设置刻度标签（针对单组样本，ax.set_xticklabels 前需确保有 ticks）
    ax.set_xticks(range(n_samples))
    ax.set_xticklabels(xlabels, rotation=rotation, ha='right')

    # -------------------------------
    # 根据分组上色 (X轴标签颜色)
    # -------------------------------
    xlabel_colors = ["black"] * n_samples
    if groups is not None:
        if len(groups) != n_samples:
            raise ValueError("groups 长度必须与样本数量一致")

        if group_colors is None:
            unique_groups = list(dict.fromkeys(groups))
            cmap_group = plt.get_cmap("tab10")
            group_colors = {g: cmap_group(i) for i, g in enumerate(unique_groups)}

        xlabel_colors = [group_colors[g] for g in groups]

    for label, c in zip(ax.get_xticklabels(), xlabel_colors):
        label.set_color(c)

    # -------------------------------
    # 第一个图例：类型
    # -------------------------------
    legend_types = df_prop.index.tolist()
    logger.info(legend_types)
    cmap_types = plt.get_cmap(colormap)
    types_patches = [
        mpatches.Patch(color=cmap_types(i / max(len(legend_types)-1, 1)), label=legend_types[i])
        for i in range(len(legend_types))
    ]

    legend1 = ax.legend(
        handles=types_patches,
        title=legend_title_type,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=legend_fontsize,
        title_fontsize=legend_title_fontsize,
        frameon=False
    )
    ax.add_artist(legend1)

    # -------------------------------
    # 第二个图例：样本分组
    # -------------------------------
    if groups is not None:
        # 保持分组出现的原始顺序
        unique_groups = list(dict.fromkeys(groups))
        group_patches = [
            mpatches.Patch(color=group_colors[g], label=g)
            for g in unique_groups
        ]
        ax.legend(
            handles=group_patches,
            title=legend_title_group,
            bbox_to_anchor=(1.02, 0.3),
            loc="upper left",
            fontsize=legend_fontsize,
            title_fontsize=legend_title_fontsize,
            frameon=False
        )

    # -------------------------------
    # 轴细节优化
    # -------------------------------
    ax.set_xlabel(xlabel, fontsize=legend_title_fontsize)
    ax.set_ylabel(ylabel, fontsize=legend_title_fontsize)
    ax.set_title(title, fontsize=legend_title_fontsize + 2)
    
    # 移除上方和右侧边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300) # bbox_inches='tight' is not allowed
        plt.close(fig)
    else:
        plt.show()


def plot_stacking_bar_count(
    df_counts: pd.DataFrame,
    xlabels: List[str] = None,
    groups: List[str] = None,
    group_colors: Dict[str, str] = None,
    title: str = "Mutation Distribution (Proportion)",
    xlabel: str = "Sample",
    ylabel: str = "Proportion",
    legend_title_type: str = "Mutation Type",
    legend_title_group: str = "Sample Group",
    save_path: Union[str, Path] = None,
    legend_width: float = 0.25,
    figsize: tuple = (12, 6),
    legend_fontsize: int = 13,
    legend_title_fontsize: int = 16,
    rotation: int = 45,
    colormap: str = "tab20",
    # ========= 新增参数 =========
    show_block_counts: bool = False,
    show_total_counts: bool = False,
    block_count_fmt: str = "{count}",
    total_count_fmt: str = "n={total}",
):
    """
    绘制突变分布的堆叠柱状图（比例），支持：
    - 色块内绝对计数标注
    - 柱子顶部总计数标注
    - 双图例（突变类型 / 样本分组）
    """

    # -------------------------------
    # 1. 比例转化
    # -------------------------------
    logger.info(f"\n{df_counts.head()}")
    df_prop = df_counts.div(df_counts.sum(axis=0).replace(0, 1), axis=1)
    logger.info(f"\n{df_prop.head()}")

    fig, ax = plt.subplots(figsize=figsize)

    # 右侧预留空间给图例
    fig.subplots_adjust(right=1 - legend_width)

    # -------------------------------
    # 2. 绘制堆叠柱状图
    # -------------------------------
    df_prop.T.plot(
        kind="bar",
        stacked=True,
        colormap=colormap,
        width=0.8,
        ax=ax,
        legend=False
    )

    n_samples = df_counts.shape[1]

    # -------------------------------
    # 3. X 轴刻度标签
    # -------------------------------
    if xlabels is None:
        xlabels = df_counts.columns.tolist()

    if len(xlabels) != n_samples:
        raise ValueError("xlabels 长度必须与样本数量一致")

    ax.set_xticks(range(n_samples))
    ax.set_xticklabels(xlabels, rotation=rotation, ha="right")

    # -------------------------------
    # 4. 根据分组给 X 轴标签上色
    # -------------------------------
    xlabel_colors = ["black"] * n_samples

    if groups is not None:
        if len(groups) != n_samples:
            raise ValueError("groups 长度必须与样本数量一致")

        if group_colors is None:
            unique_groups = list(dict.fromkeys(groups))
            cmap_group = plt.get_cmap("tab10")
            group_colors = {g: cmap_group(i) for i, g in enumerate(unique_groups)}

        xlabel_colors = [group_colors[g] for g in groups]

    for label, c in zip(ax.get_xticklabels(), xlabel_colors):
        label.set_color(c)

    # -------------------------------
    # 5. 色块计数 / 柱子总数标注
    # -------------------------------
    if show_block_counts or show_total_counts:
        df_counts_T = df_counts.T   # 行：样本，列：类型
        df_prop_T = df_prop.T

        n_types = df_counts.shape[0]

        # -------------------------------
        # 色块内部绝对计数（修正版）
        # -------------------------------
        if show_block_counts:
            df_counts_T = df_counts.T   # 行：sample，列：type
            df_prop_T = df_prop.T

            n_samples = df_counts_T.shape[0]
            n_types = df_counts_T.shape[1]

            patch_idx = 0
            for i_type in range(n_types):
                for i_sample in range(n_samples):
                    patch = ax.patches[patch_idx]
                    patch_idx += 1

                    height = patch.get_height()
                    if height <= 0:
                        continue

                    count = df_counts_T.iloc[i_sample, i_type]
                    prop = df_prop_T.iloc[i_sample, i_type]

                    x = patch.get_x() + patch.get_width() / 2
                    y = patch.get_y() + height / 2

                    ax.text(
                        x,
                        y,
                        block_count_fmt.format(count=count, prop=prop),
                        ha="center",
                        va="center",
                        fontsize=10,
                    )

        # ---- 每根柱子的总计数 ----
        if show_total_counts:
            totals = df_counts.sum(axis=0).values
            for i, total in enumerate(totals):
                ax.text(
                    i,
                    1.02,
                    total_count_fmt.format(total=total),
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    fontweight="bold",
                    transform=ax.get_xaxis_transform(),
                )

    # -------------------------------
    # 6. 图例 1：突变类型
    # -------------------------------
    legend_types = df_prop.index.tolist()
    cmap_types = plt.get_cmap(colormap)

    types_patches = [
        mpatches.Patch(
            color=cmap_types(i / max(len(legend_types) - 1, 1)),
            label=legend_types[i]
        )
        for i in range(len(legend_types))
    ]

    legend1 = ax.legend(
        handles=types_patches,
        title=legend_title_type,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=legend_fontsize,
        title_fontsize=legend_title_fontsize,
        frameon=False
    )
    ax.add_artist(legend1)

    # -------------------------------
    # 7. 图例 2：样本分组
    # -------------------------------
    if groups is not None:
        unique_groups = list(dict.fromkeys(groups))
        group_patches = [
            mpatches.Patch(color=group_colors[g], label=g)
            for g in unique_groups
        ]

        ax.legend(
            handles=group_patches,
            title=legend_title_group,
            bbox_to_anchor=(1.02, 0.3),
            loc="upper left",
            fontsize=legend_fontsize,
            title_fontsize=legend_title_fontsize,
            frameon=False
        )

    # -------------------------------
    # 8. 轴与样式优化
    # -------------------------------
    ax.set_xlabel(xlabel, fontsize=legend_title_fontsize)
    ax.set_ylabel(ylabel, fontsize=legend_title_fontsize)
    ax.set_title(title, fontsize=legend_title_fontsize + 2)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # -------------------------------
    # 9. 保存或显示
    # -------------------------------
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()


def plot_biax_bar(
        local_alignment_file: str,
        global_alignment_file: str,
        out_png: str = "alignment_comparison.png"
    ):
    """
    Plot alignment comparison using dual y-axis (Local vs Global).
    """

    # 读取数据
    local_df = pd.read_csv(local_alignment_file, sep="\t")
    global_df = pd.read_csv(global_alignment_file, sep="\t")

    # 合并数据
    df = pd.merge(local_df, global_df, on="AminoAcid", suffixes=("_Local", "_Global"))

    # 按 global similarity 排序
    df = df.sort_values("Similarity_Global", ascending=False)

    # x轴位置
    x = np.arange(len(df["AminoAcid"]))
    width = 0.35

    # 创建双y轴
    fig, ax1 = plt.subplots(figsize=(10,5))
    ax2 = ax1.twinx()

    # 绘制柱状图
    bars1 = ax1.bar(
        x - width/2,
        df["Similarity_Local"],
        width,
        label="Local Alignment",
        color="steelblue",
        alpha=0.8
    )

    bars2 = ax2.bar(
        x + width/2,
        df["Similarity_Global"],
        width,
        label="Global Alignment",
        color="orange",
        alpha=0.8
    )

    # 轴标签
    ax1.set_xlabel("Amino Acid")
    ax1.set_ylabel("Local Similarity")
    ax2.set_ylabel("Global Similarity")

    # x轴标签
    ax1.set_xticks(x)
    ax1.set_xticklabels(df["AminoAcid"], rotation=45)

    # 保存图
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    
if __name__ == '__main__':
    # DEG_file = "/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell_DEG.tsv"
    # rmsk_file = "/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz"
    # TEfamily(DEG_file,rmsk_file,"a")

    gsea_result = pd.read_csv("/disk5/luosg/scRNAseq/output/combine/Intestine/GSEA/table/combined_groupE2_chang_10XSC3_Intestinal_stem_cell-combined_groupCKO_chang_10XSC3_Intestinal_stem_cell_DEG_gene_gsea.tsv",sep="\t")
    GSEA_BarPlot_Enhanced(
        gsea_result=gsea_result,
        outplot="/disk5/luosg/scRNAseq/output/combine/Intestine/GSEA/plotTop/Plasmacytoid_dendritic_cells_E2_CKO_gsea.png",
        sort_by="pval",
        top_flag=True,
        top_number=5
    )
