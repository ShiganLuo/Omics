import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu,chi2_contingency,fisher_exact
import matplotlib.patches as mpatches
from typing import List, Union, Dict, Optional
import logging
from scipy.interpolate import make_interp_spline
from scipy.ndimage import gaussian_filter1d
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


def plot_sv_type_barplot(
    summary_df: pd.DataFrame,
    outpng: str,
    xlabel: str,
    ylabel: str,
    figsize: tuple = (6, 4),
    pivot_index: str = "svtype",
    pivot_columns: str = "group",
    pivot_values: str = "count"
):
    """
    Plot a grouped bar chart of SV type distribution (wide layout).

    Pivots a long-format DataFrame into wide format and draws a side-by-side
    bar chart comparing structural variant counts across sample groups.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Long-format DataFrame containing columns for pivot index, columns,
        and values.
    outpng : str
        Output file path for the saved PNG image.
    xlabel : str
        Label for the x-axis (e.g. ``"SV Type"``).
    ylabel : str
        Label for the y-axis (e.g. ``"Count"`` or ``"Number of SVs"``).
    figsize : tuple, optional
        Figure size ``(width, height)`` in inches. Default is ``(6, 4)``.
    pivot_index : str, optional
        Column name used as pivot table row index. Default is ``"svtype"``.
    pivot_columns : str, optional
        Column name used as pivot table columns. Default is ``"group"``.
    pivot_values : str, optional
        Column name used as pivot table values. Default is ``"count"``.

    Returns
    -------
    None
        Saves the figure to *outpng* at 300 DPI and closes the canvas.
    """
    os.makedirs(os.path.dirname(outpng), exist_ok=True)
    
    # 将长格式转换为宽格式绘图
    pivot = (
        summary_df.pivot(index=pivot_index, columns=pivot_columns, values=pivot_values)
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=figsize)
    pivot.plot(kind="bar", ax=ax, edgecolor="black")

    # 设置轴标签
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # 视觉美化：移除多余边框线，增加参考网格
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close(fig)


def plot_stacking_bar(
    df_counts: pd.DataFrame,
    xlabels: Optional[List[str]] = None,
    groups: Optional[List[str]] = None,
    group_colors: Optional[Dict[str, str]] = None,
    title: str = "Mutation Distribution (Proportion)",
    xlabel: str = "Sample",
    ylabel: str = "Proportion",
    legend_title_type: str = "Mutation Type",
    legend_title_group: str = "Sample Group",
    save_path: Optional[Union[str, Path]] = None,
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
    Plot a stacked bar chart of mutation distribution (proportions).

    Supports in-block absolute count labels, per-bar total count labels,
    and dual legends for mutation types and sample groups.

    Parameters
    ----------
    df_counts : pd.DataFrame
        Count matrix with mutation types as rows and samples as columns.
    xlabels : list of str, optional
        Custom x-axis tick labels. If ``None``, uses column names of
        *df_counts*.
    groups : list of str, optional
        Sample group assignments (length must equal number of columns).
    group_colors : dict, optional
        Mapping from group name to color. Auto-generated if ``None``.
    title : str, optional
        Plot title. Default is ``"Mutation Distribution (Proportion)"``.
    xlabel : str, optional
        X-axis label. Default is ``"Sample"``.
    ylabel : str, optional
        Y-axis label. Default is ``"Proportion"``.
    legend_title_type : str, optional
        Title for the mutation-type legend. Default is ``"Mutation Type"``.
    legend_title_group : str, optional
        Title for the sample-group legend. Default is ``"Sample Group"``.
    save_path : str or Path, optional
        If provided, figure is saved to this path; otherwise displayed.
    legend_width : float, optional
        Fraction of figure width reserved for legends. Default is ``0.25``.
    figsize : tuple, optional
        Figure size ``(width, height)``. Default is ``(12, 6)``.
    legend_fontsize : int, optional
        Font size for legend labels. Default is ``13``.
    legend_title_fontsize : int, optional
        Font size for legend titles. Default is ``16``.
    rotation : int, optional
        Rotation angle for x-tick labels. Default is ``45``.
    colormap : str, optional
        Matplotlib colormap name. Default is ``"tab20"``.
    show_block_counts : bool, optional
        If ``True``, annotate each block with its absolute count.
    show_total_counts : bool, optional
        If ``True``, annotate each bar top with its total count.
    block_count_fmt : str, optional
        Format string for block count labels. Default is ``"{count}"``.
    total_count_fmt : str, optional
        Format string for total count labels. Default is ``"n={total}"``.

    Returns
    -------
    None
        Displays the plot or saves it to *save_path*.
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


def plot_sv_length_boxplot(
    data_dict: Dict[str, pd.DataFrame],
    svlen_col: str,
    outpng: str,
    ylabel: str,
    figsize: tuple = (4, 5),
    jitter_width: float = 0.08,
    point_alpha: float = 0.6,
    large_sv_threshold: int = 10_000,
):
    """
    Plot grouped boxplots of SV lengths with jittered points and significance brackets.

    Filters SVs below *large_sv_threshold*, applies log10 transformation,
    performs pairwise Mann-Whitney U tests, and annotates significance.

    Parameters
    ----------
    data_dict : dict of str to pd.DataFrame
        Mapping of group labels to DataFrames containing SV records.
    svlen_col : str
        Column name for SV length.
    outpng : str
        Output file path for the saved PNG image.
    ylabel : str
        Label for the y-axis.
    figsize : tuple, optional
        Figure size ``(width, height)``. Default is ``(4, 5)``.
    jitter_width : float, optional
        Standard deviation of jitter noise for strip plot. Default is ``0.08``.
    point_alpha : float, optional
        Transparency of jitter points. Default is ``0.6``.
    large_sv_threshold : int, optional
        Minimum SV length to include (bp). Default is ``10000``.

    Returns
    -------
    None
        Saves the figure to *outpng* at 300 DPI and closes the canvas.
    """
    os.makedirs(os.path.dirname(outpng), exist_ok=True)

    # -------------------------
    # 数据筛选和转换
    # -------------------------
    filtered_data = {}
    log_lengths = {}

    for group, df in data_dict.items():
        filtered_df = df.loc[df[svlen_col] >= large_sv_threshold, :]
        filtered_data[group] = filtered_df
        log_lengths[group] = np.log10(filtered_df[svlen_col] + 1)

    # -------------------------
    # 统计检验（原始长度）
    # -------------------------
    group_labels = list(data_dict.keys())
    p_values = {}

    for i, group1 in enumerate(group_labels):
        for group2 in group_labels[i + 1:]:
            raw1 = filtered_data[group1][svlen_col]
            raw2 = filtered_data[group2][svlen_col]

            if len(raw1) > 0 and len(raw2) > 0:
                _, p_value = mannwhitneyu(raw1, raw2, alternative="two-sided")
                p_values[(group1, group2)] = p_value

    def p_to_label(p):
        if p < 1e-3:
            return "***"
        elif p < 1e-2:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    # -------------------------
    # 作图
    # -------------------------
    fig, ax = plt.subplots(figsize=figsize)

    box = ax.boxplot(
        log_lengths.values(),
        labels=group_labels,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(edgecolor="black"),
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
    )

    colors = plt.cm.tab10.colors[:len(group_labels)]
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # 抖点
    for i, (group, lengths) in enumerate(log_lengths.items(), start=1):
        x_jitter = np.random.normal(i, jitter_width, size=len(lengths))
        ax.scatter(
            x_jitter, lengths,
            s=10, marker=".",
            color=colors[i - 1],
            alpha=point_alpha, linewidths=0
        )

    # -------------------------
    # 统计横盖（bracket）
    # -------------------------
    y_max = max(max(lengths) for lengths in log_lengths.values())
    h = 0.03
    y = y_max + h

    for (group1, group2), p_value in p_values.items():
        x1 = group_labels.index(group1) + 1
        x2 = group_labels.index(group2) + 1

        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color="black")
        ax.text(
            (x1 + x2) / 2,
            y + h * 1.1,
            p_to_label(p_value),
            ha="center",
            va="bottom",
            fontsize=9
        )
        y += h * 2  # 增加高度避免重叠

    # -------------------------
    # 样式（无标题）
    # -------------------------
    ax.set_ylabel(ylabel)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close(fig)




def plot_large_sv_barplot(
    summary_df: pd.DataFrame,
    outpng: str,
    size_threshold: int,
    title: str,
    ylabel: str,
    figsize: tuple = (6, 4),
):
    """
    Plot a grouped bar chart of SV counts filtered by a size threshold.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Long-format DataFrame with columns ``"svtype"``, ``"group"``, and
        ``"count"``.
    outpng : str
        Output file path for the saved PNG image.
    size_threshold : int
        SV size cutoff used in the plot title.
    title : str
        Plot title.
    ylabel : str
        Label for the y-axis.
    figsize : tuple, optional
        Figure size ``(width, height)``. Default is ``(6, 4)``.

    Returns
    -------
    None
        Saves the figure to *outpng* at 300 DPI and closes the canvas.
    """
    os.makedirs(os.path.dirname(outpng), exist_ok=True)
    pivot = (
        summary_df.pivot(index="svtype", columns="group", values="count")
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=figsize)
    pivot.plot(kind="bar", ax=ax, edgecolor="black")

    ax.set_title(title)
    ax.set_xlabel("SV type")
    ax.set_ylabel(ylabel)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close(fig)


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


def plot_multi_smooth_curves(
        datasets: List[Dict],
        outfig: str,
        x_label: str = "SV Length (bp)", 
        y_label: str = "Proportion", # 通常分箱后会是比例
        title: str = "SV Length Distribution",
        smooth_sigma: float = 2.0,  # 增加此值可以变得更平滑
        highlight_x_values: List[float] = None # 新增参数：用于绘制竖直线的 X 值列表
) -> plt.Figure:
    """
    Plot multiple smoothed SV length distribution curves on a log-scale x-axis.

    Uses Gaussian filtering on log-spaced interpolated data to produce smooth
    curves, with optional vertical reference lines at specified x positions.

    Parameters
    ----------
    datasets : list of dict
        Each dict must contain ``'x'`` (list of bin midpoints) and ``'y'``
        (list of proportions). Optional keys: ``'label'``, ``'color'``.
    outfig : str
        Output file path. If ``None``, the figure is only displayed.
    x_label : str, optional
        X-axis label. Default is ``"SV Length (bp)"``.
    y_label : str, optional
        Y-axis label. Default is ``"Proportion"``.
    title : str, optional
        Plot title. Default is ``"SV Length Distribution"``.
    smooth_sigma : float, optional
        Gaussian filter sigma controlling smoothness. Default is ``2.0``.
    highlight_x_values : list of float, optional
        X positions at which to draw vertical dashed reference lines.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object.
    """
    fig = plt.figure(figsize=(10, 6))
    ax = fig.gca() # 获取当前 Axes 对象，方便后面添加注释

    for data in datasets:
        x = np.array(data['x'], dtype=float)
        y = np.array(data['y'], dtype=float)
        
        # 过滤非法值
        mask = (x > 0) & (y >= 0)
        x, y = x[mask], y[mask]
        if len(x) < 2: continue
            
        label = data.get('label', 'Unnamed')
        color = data.get('color', None) 

        # 1. 在 Log 空间对数据进行高精度重采样（生成 1000 个点）
        log_x = np.log10(x)
        # 确保 x 轴范围涵盖所有数据点
        min_log_x = log_x.min() if len(log_x) > 0 else np.log10(1) # 至少从 1 开始
        max_log_x = log_x.max() if len(log_x) > 0 else np.log10(1000000) # 至少到 1Mb
        log_x_new = np.linspace(min_log_x, max_log_x, 1000)
        
        # 2. 首先通过线性插值获取密集点
        y_interp = np.interp(log_x_new, log_x, y)
        
        # 3. 【核心步骤】使用高斯平滑
        y_smooth = gaussian_filter1d(y_interp, sigma=smooth_sigma)
        
        # 4. 确保非负性
        y_smooth = np.maximum(y_smooth, 0)

        # 5. 绘图
        ax.plot(np.power(10, log_x_new), y_smooth, label=label, color=color, linewidth=2.5)

    # --- 新增功能：绘制竖直线 ---
    if highlight_x_values:
        for val in highlight_x_values:
            ax.axvline(x=val, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
            # 可以在这里添加文本标签，但需要手动调整位置
            if val < 1000:
                label_text = f"{val:.0f} bp"
            else:
                # 如果是整数kb则不带小数，否则带1位小数
                label_text = f"{val/1000:.1f} kb".replace(".0 ", " ")
            ax.text(val, ax.get_ylim()[1]*0.95, label_text, 
                    rotation=0, va='top', ha='right', color='gray', fontsize=9)


    # 坐标轴设置
    ax.set_xscale('log') 
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both", linestyle=':', alpha=0.5)
    plt.tight_layout()
    
    if outfig:
        plt.savefig(outfig, dpi=300, bbox_inches='tight')
    
    plt.show()
    return fig