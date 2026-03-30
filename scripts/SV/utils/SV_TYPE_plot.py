from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu,chi2_contingency
import matplotlib.patches as mpatches
from typing import List, Union, Dict
import logging
from scipy.interpolate import make_interp_spline
from scipy.ndimage import gaussian_filter1d
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


def plot_sv_type_barplot(
    summary_df: pd.DataFrame,
    outpng: Path,
    xlabel: str,
    ylabel: str,
    figsize: tuple = (6, 4),
    pivot_index: str = "svtype",
    pivot_columns: str = "group",
    pivot_values: str = "count"
):
    """
    根据汇总数据绘制 SV 类型分布的分组柱状图（宽幅布局）。

    该函数将输入的长格式 (Long-format) DataFrame 转换为宽格式 (Wide-format)，
    并生成一张对比不同样本组 (Group) 间结构变异类型 (SV Type) 数量的条形图。

    功能步骤:
    1. 自动创建输出目录（如果不存在）。
    2. 执行透视表操作 (Pivot)，将变异类型设为行，样本组设为列。
    3. 绘制分组柱状图，并进行精简的视觉美化（去除上、右边框）。
    4. 以 300 DPI 的高分辨率保存图像。

    Args:
        summary_df (pd.DataFrame): 包含统计结果的 DataFrame，需包含透视所需的索引、列和值。
        outpng (Path): 图片保存的完整路径（包含文件名，建议以 .png 结尾）。
        xlabel (str): X 轴标题（通常为 "SV Type"）。
        ylabel (str): Y 轴标题（通常为 "Count" 或 "Number of SVs"）。
        figsize (tuple, optional): 画布尺寸 (宽, 高)。 默认为 (6, 4)。
        pivot_index (str, optional): 用于透视表行索引的列名。 默认为 "svtype"。
        pivot_columns (str, optional): 用于透视表分组列的列名。 默认为 "group"。
        pivot_values (str, optional): 用于透视表统计值的列名。 默认为 "count"。

    Returns:
        None: 该函数直接保存文件到磁盘并关闭画布。
        
    Example:
        >>> summary = pd.DataFrame({
        ...     'svtype': ['DEL', 'INS', 'DEL', 'INS'],
        ...     'group': ['DMSO', 'DMSO', 'PlaB', 'PlaB'],
        ...     'count': [800, 1200, 863, 1475]
        ... })
        >>> plot_sv_type_barplot(summary, Path("reports/sv_dist.png"), "SV Type", "Count")
    """
    outpng.parent.mkdir(parents=True, exist_ok=True)
    
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


def plot_sv_length_boxplot(
    ctrl_df: pd.DataFrame,
    exp_df: pd.DataFrame,
    svlen_col: str,
    group_labels: tuple,
    outpng: Path,
    ylabel: str,
    figsize: tuple = (4, 5),
    jitter_width: float = 0.08,
    point_alpha: float = 0.6,
    large_sv_threshold: int = 10_000,
):
    outpng.parent.mkdir(parents=True,exist_ok=True)
    # -------------------------
    # 数据筛选
    # -------------------------
    ctrl_df = ctrl_df.loc[ctrl_df[svlen_col] >= large_sv_threshold, :]
    exp_df = exp_df.loc[exp_df[svlen_col] >= large_sv_threshold, :]

    ctrl_raw = ctrl_df[svlen_col]
    exp_raw = exp_df[svlen_col]

    ctrl_len = np.log10(ctrl_raw + 1)
    exp_len = np.log10(exp_raw + 1)

    # -------------------------
    # 统计检验（原始长度）
    # -------------------------
    p_value = np.nan
    if len(ctrl_raw) > 0 and len(exp_raw) > 0:
        _, p_value = mannwhitneyu(
            ctrl_raw,
            exp_raw,
            alternative="two-sided"
        )

    def p_to_label(p):
        if p < 1e-3:
            return "***"
        elif p < 1e-2:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    p_label = p_to_label(p_value)

    # -------------------------
    # 作图
    # -------------------------
    fig, ax = plt.subplots(figsize=figsize)

    box = ax.boxplot(
        [ctrl_len, exp_len],
        labels=list(group_labels),
        patch_artist=True,
        showfliers=False,
        boxprops=dict(edgecolor="black"),
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
    )

    colors = ["#4C72B0", "#DD8452"]
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # 抖点
    x_ctrl = np.random.normal(1, jitter_width, size=len(ctrl_len))
    x_exp  = np.random.normal(2, jitter_width, size=len(exp_len))

    ax.scatter(
        x_ctrl, ctrl_len,
        s=10, marker=".",
        color=colors[0],
        alpha=point_alpha, linewidths=0
    )
    ax.scatter(
        x_exp, exp_len,
        s=10, marker=".",
        color=colors[1],
        alpha=point_alpha, linewidths=0
    )

    # -------------------------
    # 统计横盖（bracket）
    # -------------------------
    y_max = max(ctrl_len.max(), exp_len.max())
    h = 0.03
    y = y_max + h

    ax.plot([1, 1, 2, 2], [y, y + h, y + h, y],
            lw=1.2, color="black")

    ax.text(
        1.5,
        y + h * 1.1,
        f"{p_label}",
        ha="center",
        va="bottom",
        fontsize=9
    )

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
    outpng: Path,
    size_threshold: int,
    title: str,
    ylabel: str,
    figsize: tuple = (6, 4),
):
    outpng.parent.mkdir(parents=True,exist_ok=True)
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


def plot_svtype_comparison(
    df,
    out_png,
    group_col="group",
    svtype_col="svtype",
    count_col="count",
    group_order=("Control", "Experiment"),
    svtype_order=("BND", "DEL", "DUP", "INS", "INV"),
    legend_map=None,
    figsize=(9, 5),
    ylabel="SV count",
    dpi=300,
):

    if legend_map is None:
        legend_map = {g: g for g in group_order}

    pivot = (
        df.pivot(index=svtype_col, columns=group_col, values=count_col)
        .reindex(svtype_order)
        .fillna(0)
    )

    def p_to_star(p):
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

    stars = {}
    for sv in pivot.index:
        g1, g2 = group_order
        table = [
            [pivot.loc[sv, g1], pivot[g1].sum() - pivot.loc[sv, g1]],
            [pivot.loc[sv, g2], pivot[g2].sum() - pivot.loc[sv, g2]],
        ]
        _, p, _, _ = chi2_contingency(table)
        stars[sv] = p_to_star(p)

    sig_sv = [sv for sv, s in stars.items() if s != "ns"]
    low_max = (
        pivot.loc[sig_sv].values.max() * 1.15
        if sig_sv else
        np.median(pivot.values) * 1.5
    )

    global_max = pivot.values.max()
    high_min = low_max * 1.1
    high_max = global_max * 1.15   # ⬅️ 预留空间给显著帽子

    x = np.arange(len(pivot.index))
    width = 0.36

    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, sharex=True,
        figsize=figsize,
        gridspec_kw={"height_ratios": [1, 3]},
    )

    colors = {
        group_order[0]: "#4C72B0",
        group_order[1]: "#DD8452",
    }

    for ax in (ax_top, ax_bottom):
        ax.bar(x - width / 2, pivot[group_order[0]], width,
               color=colors[group_order[0]], label=legend_map[group_order[0]])
        ax.bar(x + width / 2, pivot[group_order[1]], width,
               color=colors[group_order[1]], label=legend_map[group_order[1]])

    ax_bottom.set_ylim(0, low_max)
    ax_top.set_ylim(high_min, high_max)

    # ---------- fixed broken axis (LEFT ONLY, SAFE) ----------
    d = 0.008

    # top panel: bottom edge
    ax_top.plot(
        (-d, +d),
        (-d, +d),
        transform=ax_top.transAxes,
        color="black",
        clip_on=False,
    )

    # bottom panel: top edge
    ax_bottom.plot(
        (-d, +d),
        (1 - d, 1 + d),
        transform=ax_bottom.transAxes,
        color="black",
        clip_on=False,
    )


    # ---------- significance (VISUALLY CONSISTENT, SAFE) ----------
    LEG_PT = 8     # 显著腿高度（物理单位）
    TEXT_PT = 3

    fig.canvas.draw()  # 保证 transform 可用

    for i, sv in enumerate(pivot.index):
        y1 = pivot.loc[sv, group_order[0]]
        y2 = pivot.loc[sv, group_order[1]]
        y_base = max(y1, y2)

        ax = ax_bottom if y_base <= low_max else ax_top

        # data -> display
        trans = ax.transData
        inv = ax.transData.inverted()

        _, y_disp = trans.transform((0, y_base))
        _, y_hat = inv.transform((0, y_disp + LEG_PT))
        _, y_text = inv.transform((0, y_disp + LEG_PT + TEXT_PT))

        x1 = x[i] - width / 2
        x2 = x[i] + width / 2

        ax.plot([x1, x1], [y1, y_hat], lw=1.2, c="black")
        ax.plot([x2, x2], [y2, y_hat], lw=1.2, c="black")
        ax.plot([x1, x2], [y_hat, y_hat], lw=1.2, c="black")

        ax.text(
            x[i], y_text, stars[sv],
            ha="center", va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(pivot.index)
    ax_bottom.set_ylabel(ylabel)

    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_top.spines["bottom"].set_visible(False)

    for ax in (ax_top, ax_bottom):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax_top.legend(frameon=False)

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
    极致平滑版本：使用高斯滤波处理数据，解决震荡和负数问题
    新增功能：可在指定 X 轴位置绘制竖直线
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