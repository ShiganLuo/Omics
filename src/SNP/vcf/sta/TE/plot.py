import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

def plot_te_family_heatmap_clustered(snv_file, group_file, outpng, top_n=30, distance_col="distance"):
    snv_df = pd.read_csv(snv_file, sep="\t")
    group_df = pd.read_csv(group_file, sep="\t")

    df = snv_df.merge(group_df, on="sample", how="left")

    df["family"] = df["t_attr"].astype(str).str.extract(r'family_id "([^"]+)"')[0].fillna("NA")

    df_te = df[df[distance_col] == 0]

    top_families = df_te["family"].value_counts().head(top_n).index.tolist()
    if not top_families:
        print("No TE families found.")
        return

    df_te = df_te[df_te["family"].isin(top_families)]

    pivot = df_te.groupby(["family", "sample"]).size().unstack(fill_value=0)

    pivot_norm = pivot.div(pivot.sum(axis=0).replace(0,1), axis=1)

    sample_groups = group_df.set_index("sample").loc[pivot_norm.columns, "group"]
    sorted_samples = sample_groups.sort_values().index
    pivot_norm = pivot_norm[sorted_samples]

    link = linkage(pivot_norm, method='average', metric='correlation')
    row_order = leaves_list(link)
    pivot_norm = pivot_norm.iloc[row_order, :]

    group_lut = dict(zip(sample_groups.unique(), sns.color_palette("Set2", len(sample_groups.unique()))))
    col_colors = sample_groups.map(group_lut)

    sns.set(font_scale=0.8)
    g = sns.clustermap(
        pivot_norm,
        row_cluster=False,
        col_cluster=False,
        col_colors=col_colors,
        cmap="rocket_r",
        linewidths=0.5,
        linecolor="gray",
        annot=True,
        fmt=".2f",
        cbar_kws={"label":"Proportion"},
        figsize=(max(6, len(pivot_norm.columns)*0.4), max(4, len(pivot_norm.index)*0.4))
    )

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
    df["family"] = (
        df["t_attr"].astype(str)
        .str.extract(r'family_id "([^"]+)"')[0]
        .fillna("NA")
    )

    top = df[df["distance"] == 0]["family"].value_counts().head(top_n).index.tolist()
    if not top:
        return

    sub = df[df["family"].isin(top)]
    tab = sub.groupby(["group", "family"]).size().unstack(fill_value=0)

    tab = tab.loc[:, tab.sum(axis=0).sort_values(ascending=False).index]

    colors = sns.color_palette("tab20", n_colors=len(tab.columns))

    ax = tab.plot(
        kind="bar",
        stacked=True,
        figsize=(10, 5),
        width=0.8,
        color=colors
    )

    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], tab.columns[::-1], title="Family",
              bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

    plt.xlabel("Group", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.title(f"Top {top_n} Families in TE-overlapping SNVs by Group", fontsize=14)

    plt.tight_layout()
    plt.savefig(outpng, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/TE/All_samples.distance.tsv",sep = "\t")
    plot_family_group_barstack_from_df(df,"/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/TE/family_group_stack.png")
