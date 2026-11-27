import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import warnings
from typing import Dict

def p_to_star(p):
    if p >= 0.05:
        return "ns"
    elif p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    else:
        return "*"


def plot_all_datasets_compare_groups(
        integrated_burden_file: str,
        group_files: Dict[str, str],
        outpng:str = "all_datasets_violin_compare.png",
        figsize:tuple = (14,6),
        violin_inner:str = "box",
        point_jitter: float = 0.08,
        dpi:int = 300,
        burden_col:str = "mutation_burden_per_Mb",
        ylim_pad:float = 0.12
    ) -> pd.DataFrame:

    df_burden = pd.read_csv(integrated_burden_file, sep="\t", engine="python")
    df_burden.columns = [c.strip() for c in df_burden.columns]
    colmap = {c.lower(): c for c in df_burden.columns}

    if "sample" not in colmap:
        raise ValueError("突变负担文件必须包含 sample 列")
    if burden_col.lower() not in colmap:
        raise ValueError(f"突变负担文件未找到 burden 列: {burden_col}")

    sample_col = colmap["sample"]
    burden_real_col = colmap[burden_col.lower()]

    df_burden = df_burden[[sample_col, burden_real_col]].rename(
        columns={sample_col: "sample", burden_real_col: "burden"}
    )

    merged_list = []
    for ds_name, grp_path in group_files.items():
        dfg = pd.read_csv(grp_path, sep=None, engine="python")
        dfg.columns = [c.strip() for c in dfg.columns]

        gmap = {c.lower(): c for c in dfg.columns}
        if "sample" not in gmap or "group" not in gmap:
            raise ValueError(f"{grp_path} 必须包含 sample 和 group 列")

        dfg = dfg[[gmap["sample"], gmap["group"]]].rename(
            columns={gmap["sample"]: "sample", gmap["group"]: "group"}
        )

        df = df_burden.merge(dfg, on="sample", how="inner")
        df["dataset"] = ds_name

        if df.empty:
            warnings.warn(f"数据集 {ds_name} 与 burden 文件无交集样本")
            continue

        merged_list.append(df)

    if not merged_list:
        raise ValueError("没有有效数据集")

    df_all = pd.concat(merged_list, ignore_index=True)

    df_all["group_unique"] = df_all["dataset"] + "|" + df_all["group"]
    group_order = []

    for ds_name, grp_path in group_files.items():

        dfg = pd.read_csv(grp_path, sep=None, engine="python")
        dfg.columns = [c.strip() for c in dfg.columns]
        gmap = {c.lower(): c for c in dfg.columns}

        group_col = gmap["group"]

        original_group_order = list(dict.fromkeys(dfg[group_col].tolist()))

        for g in original_group_order:
            group_unique = f"{ds_name}|{g}"
            group_order.append(group_unique)

    valid_groups = set(df_all["group_unique"].unique())
    group_order = [g for g in group_order if g in valid_groups]

    x_labels = [gu.split("|", 1)[1] for gu in group_order]

    stats = []
    for ds in sorted(df_all["dataset"].unique()):
        sub = df_all[df_all["dataset"] == ds]
        groups = sub["group"].unique()

        if len(groups) == 2:
            g1, g2 = groups
            v1 = sub[sub["group"] == g1]["burden"].dropna()
            v2 = sub[sub["group"] == g2]["burden"].dropna()

            if len(v1) >= 1 and len(v2) >= 1:
                _, p = mannwhitneyu(v1, v2, alternative="two-sided")
            else:
                p = np.nan

            stats.append((ds, g1, g2, p))
        else:
            stats.append((ds, None, None, np.nan))

    sns.set(style="whitegrid")
    plt.figure(figsize=figsize)

    ax = sns.violinplot(
        data=df_all,
        x="group_unique",
        y="burden",
        hue="dataset",
        dodge=False,
        inner=violin_inner,
        order=group_order
    )

    sns.stripplot(
        data=df_all,
        x="group_unique",
        y="burden",
        hue="dataset",
        dodge=False,
        order=group_order,
        jitter=point_jitter,
        alpha=0.6,
        size=3,
        linewidth=0
    )

    handles, labels = ax.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))

    leg = ax.legend(
        uniq.values(),
        uniq.keys(),
        title="Dataset",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        handletextpad=0.4,
        borderpad=0.3,
        handlelength=1.5
    )

    for lh in leg.legend_handles:
        try:
            lh.set_markersize(12)
            lh.set_alpha(0.8)
        except:
            pass

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(x_labels, rotation=45, ha="right")

    ymax = df_all["burden"].max()
    ymin = df_all["burden"].min()
    yrange = ymax - ymin if ymax != ymin else 1.0

    base = ymax + yrange * 0.02
    step = yrange * 0.06

    for i, (ds, g1, g2, p) in enumerate(stats):
        if g1 is None or np.isnan(p):
            continue

        star = p_to_star(p)

        try:
            x1 = group_order.index(ds + "|" + g1)
            x2 = group_order.index(ds + "|" + g2)
        except ValueError:
            continue

        left = x1
        right = x2
        y = base + i * step

        ax.plot([left, right], [y, y], color="black")
        ax.plot([left, left], [y, y-step*0.25], color="black")
        ax.plot([right, right], [y, y-step*0.25], color="black")

        ax.text((left+right)/2, y + step*0.05, star,
                ha="center", va="bottom", fontsize=13)

    ax.set_ylim(ymin - yrange * 0.05,
                base + len(stats)*step + yrange * ylim_pad)

    plt.xlabel("")
    plt.ylabel(burden_col)
    plt.title("")
    plt.tight_layout()
    plt.savefig(outpng, dpi=dpi)
    plt.close()

    print(f"Saved plot: {outpng}")
    return df_all


def plot_mutation_burden(df, outfile=None, value_col="mutation_burden_per_Mb", title="Mutation burden per Mb"):
    """
    绘制按培养条件自动分组着色的突变负荷柱状图，
    x 轴仅显示样本编号（例如 GSM4777778）
    """

    # ---- 提取培养条件，例如 SL、2iL、A2iL、LCDM、2iL-F ----
    def get_condition(sample):
        return sample.split("_")[-1].split("-")[0]

    # ---- 提取 GSM 样本编号 ----
    def get_sample_id(sample):
        return sample.split("_")[0]

    df = df.copy()
    df["condition"] = df["sample"].apply(get_condition)
    df["sample_id"] = df["sample"].apply(get_sample_id)

    conditions = df["condition"].unique()
    color_map = {cond: plt.cm.tab20(i) for i, cond in enumerate(conditions)}

    plt.figure(figsize=(14, 6))

    bar_colors = df["condition"].map(lambda x: color_map[x])

    plt.bar(df["sample_id"], df[value_col], color=bar_colors)

    plt.xticks(rotation=60, ha='right')
    plt.ylabel(value_col)
    plt.title(title)

    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[c]) for c in conditions]
    plt.legend(handles, conditions, title="Condition")

    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=300)
    else:
        plt.show()



if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/PRJNA663159/SNV/PRJNA663159_exonic_mutation_burden.csv")
    plot_mutation_burden(df,"bar.png")
