import gzip
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


########################################
# 1. 打开 vcf / vcf.gz
########################################

def open_vcf(path: str):
    """
    Open VCF or VCF.GZ file
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


########################################
# 2. 解析 pbsv VCF
########################################

def parse_pbsv_vcf(vcf_file: str) -> pd.DataFrame:
    """
    Parse pbsv SV VCF into DataFrame
    """
    records = []

    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            chrom, pos, _, _, _, _, _, info = fields[:8]
            pos = int(pos)

            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v

            svtype = info_dict.get("SVTYPE", "NA")
            end = int(info_dict.get("END", pos))

            svlen_raw = info_dict.get("SVLEN")
            if svlen_raw is not None:
                try:
                    svlen = abs(int(svlen_raw))
                except ValueError:
                    svlen = abs(end - pos)
            else:
                svlen = abs(end - pos)

            records.append(
                {
                    "chrom": chrom,
                    "start": pos,
                    "end": end,
                    "svtype": svtype,
                    "svlen": svlen,
                }
            )

    return pd.DataFrame(records)


########################################
# 3. 绘图函数（全部显式形参）
########################################

def plot_sv_type_barplot(
    summary_df: pd.DataFrame,
    outpng: Path,
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


########################################
# 4. 主分析流程
########################################

def run_pbsv_diff_analysis(
    ctrl_vcf: str,
    exp_vcf: str,
    outdir: str,
    large_sv_threshold: int,
):
    outdir = Path(outdir)
    table_dir = outdir / "table"
    table_dir.mkdir(parents=True,exist_ok=True)
    plot_dir = outdir / "plot"
    plot_dir.mkdir(parents=True,exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    ctrl_df = parse_pbsv_vcf(ctrl_vcf)
    exp_df = parse_pbsv_vcf(exp_vcf)

    ctrl_df["group"] = "Control"
    exp_df["group"] = "Experiment"

    all_df = pd.concat([ctrl_df, exp_df], ignore_index=True)

    # 保存所有 SV
    all_df.to_csv(table_dir / "all_sv_records.tsv", sep="\t", index=False)

    # SV 类型统计
    type_summary = (
        all_df.groupby(["group", "svtype"])
        .size()
        .reset_index(name="count")
    )

    type_summary.to_csv(
        table_dir / "sv_type_summary.tsv", sep="\t", index=False
    )

    # 图 1：SV 类型
    plot_sv_type_barplot(
        summary_df=type_summary,
        outpng = plot_dir / "sv_type_barplot.png",
        title="SV type comparison (pbsv)",
        ylabel="SV count",
    )

    # 图 2：SV 长度分布
    plot_sv_length_boxplot(
        ctrl_df=ctrl_df,
        exp_df=exp_df,
        svlen_col="svlen",
        group_labels=("Control", "Experiment"),
        outpng = plot_dir / f"large_sv_{large_sv_threshold}bp_boxplot.png",
        ylabel="log10(SV length + 1)",
        large_sv_threshold = large_sv_threshold
    )

    # ≥阈值的大 SV
    large_sv_df = all_df[all_df["svlen"] >= large_sv_threshold]
    large_summary = (
        large_sv_df.groupby(["group", "svtype"])
        .size()
        .reset_index(name="count")
    )

    large_summary.to_csv(
        table_dir / "large_sv_{large_sv_threshold}bp_summary.tsv",
        sep="\t",
        index=False,
    )

    if not large_summary.empty:
        plot_large_sv_barplot(
            summary_df=large_summary,
            outpng = plot_dir / f"large_sv_{large_sv_threshold}bp_barplot.png",
            size_threshold=large_sv_threshold,
            title=f"Large SV comparison (≥{large_sv_threshold} bp)",
            ylabel=f"SV count (≥{large_sv_threshold} bp)",
        )

    print("pbsv SV differential analysis finished.")
    print(f"Results saved in: {outdir}")


########################################
# 5. CLI
########################################

def main():
    parser = argparse.ArgumentParser(
        description="pbsv SV differential analysis (VCF / VCF.GZ)"
    )
    parser.add_argument("--control", required=True, help="Control pbsv VCF(.gz)")
    parser.add_argument("--experiment", required=True, help="Experiment pbsv VCF(.gz)")
    parser.add_argument("--outdir", default="pbsv_sv_diff")
    parser.add_argument("--large-sv-threshold", type=int, default=10_000)

    args = parser.parse_args()

    run_pbsv_diff_analysis(
        ctrl_vcf=args.control,
        exp_vcf=args.experiment,
        outdir=args.outdir,
        large_sv_threshold=args.large_sv_threshold,
    )


if __name__ == "__main__":
    main()
