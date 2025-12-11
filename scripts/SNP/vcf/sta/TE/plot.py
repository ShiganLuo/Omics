import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

def plot_te_family_heatmap_clustered(snv_file, group_file, outpng, top_n=30, distance_col="distance"):
    # 读取数据
    snv_df = pd.read_csv(snv_file, sep="\t")
    group_df = pd.read_csv(group_file, sep="\t")

    # 合并 group 信息
    df = snv_df.merge(group_df, on="sample", how="left")

    # 提取 family_id
    df["family"] = df["t_attr"].astype(str).str.extract(r'family_id "([^"]+)"')[0].fillna("NA")

    # 只保留落在TE内的SNV
    df_te = df[df[distance_col] == 0]

    # 取 top N 高频 family
    top_families = df_te["family"].value_counts().head(top_n).index.tolist()
    if not top_families:
        print("No TE families found.")
        return

    df_te = df_te[df_te["family"].isin(top_families)]

    # 构建 pivot 表: family × sample
    pivot = df_te.groupby(["family", "sample"]).size().unstack(fill_value=0)

    # 转换为比例
    pivot_norm = pivot.div(pivot.sum(axis=0).replace(0,1), axis=1)

    # 样本按 group 排序
    sample_groups = group_df.set_index("sample").loc[pivot_norm.columns, "group"]
    sorted_samples = sample_groups.sort_values().index
    pivot_norm = pivot_norm[sorted_samples]

    # family 聚类
    link = linkage(pivot_norm, method='average', metric='correlation')
    row_order = leaves_list(link)
    pivot_norm = pivot_norm.iloc[row_order, :]

    # group 色条
    group_lut = dict(zip(sample_groups.unique(), sns.color_palette("Set2", len(sample_groups.unique()))))
    col_colors = sample_groups.map(group_lut)

    # 绘制热图
    sns.set(font_scale=0.8)
    g = sns.clustermap(
        pivot_norm,
        row_cluster=False,  # family 已聚类
        col_cluster=False,  # 样本按 group 排序
        col_colors=col_colors,
        cmap="rocket_r",
        linewidths=0.5,
        linecolor="gray",
        annot=True,
        fmt=".2f",
        cbar_kws={"label":"Proportion"},
        figsize=(max(6, len(pivot_norm.columns)*0.4), max(4, len(pivot_norm.index)*0.4))
    )

    # 添加 group 图例
    for label in sample_groups.unique():
        g.ax_col_dendrogram.bar(0, 0, color=group_lut[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center", ncol=len(sample_groups.unique()), bbox_to_anchor=(0.5, 1.05))

    if outpng:
        plt.savefig(outpng, bbox_inches="tight")
    plt.show()




def plot_family_group_barstack_from_df(df: pd.DataFrame, outpng: str, top_n: int = 10):
    """
    绘制不同 group 中 TE 内 SNV 前 top_n 个 family 的堆叠条形图，
    保证颜色和图例对应正确。
    """
    # 提取 family
    df["family"] = (
        df["t_attr"].astype(str)
        .str.extract(r'family_id "([^"]+)"')[0]
        .fillna("NA")
    )

    # 统计 top_n family
    top = df[df["distance"] == 0]["family"].value_counts().head(top_n).index.tolist()
    if not top:
        return

    # 筛选 top family 数据
    sub = df[df["family"].isin(top)]
    tab = sub.groupby(["group", "family"]).size().unstack(fill_value=0)

    # 按总数量降序排序列
    tab = tab.loc[:, tab.sum(axis=0).sort_values(ascending=False).index]

    # 使用 seaborn colormap，确保颜色顺序与列一致
    colors = sns.color_palette("tab20", n_colors=len(tab.columns))

    # 绘制堆叠条形图
    ax = tab.plot(
        kind="bar",
        stacked=True,
        figsize=(10, 5),
        width=0.8,
        color=colors
    )

    # 手动设置 legend，确保颜色和 family 对应
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], tab.columns[::-1], title="Family",
              bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

    # 坐标轴和标题美化
    plt.xlabel("Group", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.title(f"Top {top_n} Families in TE-overlapping SNVs by Group", fontsize=14)

    plt.tight_layout()
    plt.savefig(outpng, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/TE/All_samples.distance.tsv",sep = "\t")
    plot_family_group_barstack_from_df(df,"/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/TE/family_group_stack.png")
    # plot_te_family_heatmap_clustered("/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE224794/All_samples.distance.tsv",
    #                     "/disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE224794.tsv",
    #                     outpng="/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE224794/heatmap.png")
    # plot_te_family_heatmap_clustered("/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE204801/All_samples.distance.tsv",
    #                     "/disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE204801.tsv",
    #                     outpng="/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE204801/heatmap.png")
    # plot_te_family_heatmap_clustered("/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE185005/All_samples.distance.tsv",
    #                     "/disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE185005.tsv",
    #                     outpng="/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE185005/heatmap.png")
    # plot_te_family_heatmap_clustered("/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE166216/All_samples.distance.tsv",
    #                     "/disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE166216.tsv",
    #                     outpng="/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE166216/heatmap.png")