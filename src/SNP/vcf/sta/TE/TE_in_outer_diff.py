import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon

def plot_te_inside_outside_comparison(snv_file, group_file, outpng=None, distance_col="distance"):
    # 读取数据
    snv_df = pd.read_csv(snv_file, sep="\t")
    group_df = pd.read_csv(group_file, sep="\t")
    df = snv_df.merge(group_df, on="sample", how="left")

    # 统计每个样本 TE 内外 SNV 数量
    counts = df.groupby(["sample"]).apply(
        lambda x: pd.Series({
            "TE_inside": (x[distance_col] == 0).sum(),
            "TE_outside": (x[distance_col] != 0).sum(),
            "group": x["group"].iloc[0]
        })
    ).reset_index()

    # 将数据转换为长格式
    counts_melt = counts.melt(id_vars=["sample","group"], value_vars=["TE_inside","TE_outside"],
                              var_name="TE_status", value_name="SNV_count")

    # 绘制条形图
    plt.figure(figsize=(8,6))
    ax = sns.barplot(x="TE_status", y="SNV_count", hue="group", data=counts_melt,
                     ci="sd", palette="Set2", edgecolor="black")

    # 绘制差异线（每个样本连线 TE_inside vs TE_outside）
    for _, row in counts.iterrows():
        x1, x2 = 0, 1
        y1, y2 = row["TE_inside"], row["TE_outside"]
        ax.plot([x1, x2], [y1, y2], color="gray", alpha=0.5, linewidth=1)

    # 配对 Wilcoxon 检验
    stat_results = []
    for grp in counts["group"].unique():
        grp_df = counts[counts["group"]==grp]
        stat, p = wilcoxon(grp_df["TE_inside"], grp_df["TE_outside"])
        stat_results.append((grp, stat, p))
        print(f"Group {grp}: Wilcoxon test TE_inside vs TE_outside -> statistic={stat}, p={p:.3e}")

    plt.ylabel("SNV count")
    plt.xlabel("TE status")
    # plt.title("Comparison of SNV counts inside vs outside TE")
    plt.tight_layout()

    if outpng:
        plt.savefig(outpng, bbox_inches="tight")
    plt.show()

    return stat_results

if __name__ == "__main__":
    plot_te_inside_outside_comparison(
        snv_file="/disk5/luosg/Totipotent20251031/output/SNP/vcf/TE/GSE166216/All_samples.distance.tsv",
        group_file="/disk5/luosg/Totipotent20251031/data/Totipotency/assests/GSE166216.tsv",
        outpng="TE_inside_outside_comparison.png"
    )