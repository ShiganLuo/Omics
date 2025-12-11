import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Union, Dict
import glob
import matplotlib.patches as mpatches
from scipy.stats import chi2_contingency
import numpy as np
# -----------------------------
# 函数1：统计单个 CSV 文件 Func.refGene 突变类型
# -----------------------------
def count_mutation_types(file_path: Union[str, Path]) -> pd.Series:
    df = pd.read_csv(file_path)
    if 'Func.refGene' not in df.columns:
        raise ValueError(f"CSV 文件中不存在 'Func.refGene': {file_path}")
    return df['Func.refGene'].value_counts()


# -----------------------------
# 函数2：统计多样本，支持单样本
# -----------------------------
def count_mutation_types_multiple(files: List[Union[str, Path]]) -> pd.DataFrame:
    all_counts = []
    sample_names = []

    for f in files:
        sample_name = Path(f).parent.parent.name
        counts = count_mutation_types(f)
        all_counts.append(counts)
        sample_names.append(sample_name)

    df_counts = pd.concat(all_counts, axis=1)
    df_counts.columns = sample_names

    return df_counts.fillna(0).astype(int)


# -----------------------------
# 函数3：绘图（支持 x 轴彩色标签 + 样本分组图例）
# -----------------------------
def plot_mutation_distribution_multi(
        df_counts: pd.DataFrame,
        xlabels: List[str] = None,
        groups: List[str] = None,
        group_colors: dict = None,
        title: str = "Func.refGene Mutation Distribution (Proportion)",
        comparisons: Dict[str, List[str]] = None,
        p_values: Dict[str,float] = None,
        group_map: Dict[str,List[str]] = None,
        save_path: Union[str, str] = None,
        legend_width: float = 0.25,
        figsize: tuple = (13, 8),
        legend_fontsize: int = 13,
        legend_title_fontsize: int = 15,
        legend2_start:float = 0.16
    ):

    df_prop = df_counts.div(df_counts.sum(axis=0), axis=1)

    fig, ax = plt.subplots(figsize=figsize)

    # ======== 右侧预留空间给图例 ========
    fig.subplots_adjust(right=1 - legend_width,bottom=0.2)

    # 堆叠柱状图
    df_prop.T.plot(
        kind="bar",
        stacked=True,
        colormap="tab20",
        width=0.8,
        ax=ax,
        legend=False,
        fontsize=legend_fontsize
    )

    n_samples = df_counts.shape[1]

    # -------------------------------
    # X 轴标签设置
    # -------------------------------
    if xlabels is None:
        xlabels = df_counts.columns.tolist()

    if len(xlabels) != n_samples:
        raise ValueError("xlabels 长度必须与样本数量一致")

    xlabel_colors = ["black"] * n_samples

    # -------------------------------
    # 根据分组上色
    # -------------------------------
    if groups is not None:
        if len(groups) != n_samples:
            raise ValueError("groups 长度必须与样本数量一致")

        if group_colors is None:
            unique_groups = list(dict.fromkeys(groups))
            cmap = plt.get_cmap("tab10")
            group_colors = {g: cmap(i) for i, g in enumerate(unique_groups)}

        xlabel_colors = [group_colors[g] for g in groups]

    # 设置 x 轴颜色
    for label, c in zip(ax.get_xticklabels(), xlabel_colors):
        label.set_color(c)

    # -------------------------------
    # 第一个图例：突变类型
    # -------------------------------
    mutation_types = df_prop.index.tolist()
    cmap = plt.get_cmap("tab20")
    mutation_patches = [
        mpatches.Patch(color=cmap(i), label=mutation_types[i])
        for i in range(len(mutation_types))
    ]

    legend1 = ax.legend(
        handles=mutation_patches,
        title="Mutation Type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=legend_fontsize,
        title_fontsize=legend_title_fontsize,
        frameon = False,
    )
    ax.add_artist(legend1)

    # -------------------------------
    # 第二个图例：样本分组
    # -------------------------------
    if groups is not None:
        unique_groups = sorted(set(groups))
        group_patches = [
            mpatches.Patch(color=group_colors[g], label=g)
            for g in unique_groups
        ]
        ax.legend(
            handles=group_patches,
            title="Sample Group",
            bbox_to_anchor=(1.02, legend2_start),
            loc="upper left",
            fontsize=legend_fontsize,
            title_fontsize=legend_title_fontsize,
            frameon = False,
        )
    if comparisons is not None and p_values is not None:
        ymax = df_prop.sum(axis=0).max()
        y_offset = 0.06 * ymax  # 每条横线的垂直间隔
        line_height = 0.02 * ymax  # 横线高度

        for idx, (comp_name, items) in enumerate(comparisons.items()):
            if len(items) < 2:
                continue

            # 获取柱子中心位置
            x_pos = []
            for it in items:
                if group_map is not None and it in group_map:
                    # 组内样本均值位置
                    cols = group_map[it]
                    x_center = sum([df_counts.columns.get_loc(c) for c in cols])/len(cols)
                    x_pos.append(x_center)
                else:
                    x_pos.append(df_counts.columns.get_loc(it))

            x_start = min(x_pos)
            x_end = max(x_pos)
            y = 1.05*ymax + idx*y_offset

            # 横线
            ax.plot([x_start, x_end], [y, y], color='black', linewidth=1.2)

            # 竖线短线连接到柱子上
            ax.plot([x_start, x_start], [y - line_height, y], color='black', linewidth=1.2)
            ax.plot([x_end, x_end], [y - line_height, y], color='black', linewidth=1.2)

            # 星号标记
            p = p_values.get(comp_name, 1)
            if p < 0.001:
                sig_label = "***"
            elif p < 0.01:
                sig_label = "**"
            elif p < 0.05:
                sig_label = "*"
            else:
                sig_label = "ns"

            ax.text((x_start + x_end)/2, y + 0.01*ymax, sig_label,
                    ha='center', va='bottom', fontsize=12)
    
    ax.set_xlabel("Sample", fontsize=legend_title_fontsize)
    ax.set_ylabel("Proportion", fontsize=legend_title_fontsize)
    ax.set_title(title, fontsize=legend_title_fontsize)
    ax.set_xticklabels(xlabels, rotation=45, ha='right')

    if save_path:
        plt.savefig(save_path, dpi=300)


def chi2_tests(df_counts: pd.DataFrame, comparisons: Dict[str, List[str]], 
               group_map: Dict[str, List[str]] = None) -> Dict[str,float]:
    """
    对指定的样本组合做卡方检验（基于比例）
    
    Args:
        df_counts: 行为 Func.refGene，列为样本名
        comparisons: 字典，key=比较名, value=样本列名列表或组名列表
        group_map: 当跨组比较时，指定组名到样本名列表的映射，用于计算组总和
    
    Returns:
        每个比较的 p 值字典
    """
    p_values = {}

    df_prop = df_counts
    print(df_prop)
    for comp_name, items in comparisons.items():
        if group_map is not None:
            # items 可能是组名
            samples = []
            for item in items:
                if item in group_map:
                    # 跨组比较：取组内样本平均比例
                    cols = group_map[item]
                    samples.append(df_prop[cols].mean(axis=1))
                else:
                    # 单样本
                    samples.append(df_prop[item])
            table = pd.concat(samples, axis=1)
            table.columns = items
        else:
            table = df_prop[items]

        if table.shape[1] < 2:
            raise ValueError(f"比较 {comp_name} 必须至少包含两个样本/组合")

        # 卡方检验
        chi2, p, _, _ = chi2_contingency(table)
        p_values[comp_name] = p

    return p_values

def proportion_permutation_test(df_counts: pd.DataFrame, 
                                comparisons: dict, 
                                group_map: dict = None, 
                                n_permutations: int = 10000,
                                random_state: int = 42) -> dict:
    """
    对指定样本组合做比例分布差异的非参数检验（Permutation Test）
    
    Args:
        df_counts: 行为 Func.refGene，列为样本名
        comparisons: 字典，key=比较名, value=样本列名列表或组名列表
        group_map: 跨组比较时，用组名映射到样本列表
        n_permutations: 排列次数
        random_state: 随机种子
        
    Returns:
        dict，每个比较的 p 值
    """
    rng = np.random.default_rng(random_state)
    p_values = {}
    
    # 计算比例
    df_prop = df_counts.div(df_counts.sum(axis=0), axis=1)
    
    for comp_name, items in comparisons.items():
        # 提取样本比例
        samples = []
        for item in items:
            if group_map is not None and item in group_map:
                cols = group_map[item]
                samples.append(df_prop[cols].mean(axis=1))
            else:
                samples.append(df_prop[item])
        
        # 构造一个矩阵，行=功能区，列=样本/组合
        table = pd.concat(samples, axis=1)
        
        # 统计差异：这里用每列平均值的差异总和作为 statistic
        obs_stat = 0
        n_cols = table.shape[1]
        for i in range(n_cols):
            for j in range(i+1, n_cols):
                obs_stat += np.sum(np.abs(table.iloc[:, i] - table.iloc[:, j]))
        
        # permutation
        combined = table.values.flatten()
        count = 0
        for _ in range(n_permutations):
            perm = rng.permutation(combined).reshape(table.shape)
            perm_stat = 0
            for i in range(n_cols):
                for j in range(i+1, n_cols):
                    perm_stat += np.sum(np.abs(perm[:, i] - perm[:, j]))
            if perm_stat >= obs_stat:
                count += 1
        p_values[comp_name] = count / n_permutations
    
    return p_values

if __name__ == "__main__":

    files = glob.glob(
        "/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/annotate/**/*multianno.csv",
        recursive=True
    )

    files_sorted = sorted(
        files,
        key=lambda x: (("GRCh38" not in x), ("GRCm39" not in x), x)
    )

    df_counts = count_mutation_types_multiple(files_sorted)

    print(df_counts)
    xlabels = [Path(f).parent.parent.name for f in files_sorted]

    # 假设你根据文件名区分人（GRCh38）和鼠（GRCm39）
    groups = ["Human" if "GRCh38" in f else "Mouse" for f in files_sorted]
    print(groups)
    # 可自定义颜色（也可以不传，会自动生成）
    group_colors = {
        "Human": "red",
        "Mouse": "blue"
    }
    comparisons = {
        "Mouse": ["TLSC","ciTotiSC"],    # 组内比较
        "Human": ["ci8CLC","hTBLC"],            # 组内比较
        "Mouse_vs_Human": ["Mouse","Human"]     # 跨组比较，取组总和
    }
    groups_map = {
        "Mouse": ["TLSC","ciTotiSC"],
        "Human": ["ci8CLC","hTBLC"]
    }
    p_values = proportion_permutation_test(df_counts, comparisons, group_map=groups_map)
    print(p_values)
    plot_mutation_distribution_multi(
        df_counts,
        xlabels=xlabels,
        groups=groups,
        group_colors=group_colors,
        legend_width = 0.25,
        comparisons=comparisons,
        p_values=p_values,
        group_map=groups_map,
        title="",
        save_path="/disk5/luosg/Totipotent20251031/output/SNP/vcf/intersect/plot/mut_type_prop_grouped.png"
    )
