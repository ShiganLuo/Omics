from pathlib import Path
import pandas as pd
import argparse
import logging
from utils.SV_TYPE_plot import plot_sv_type_barplot,plot_sv_length_boxplot,plot_large_sv_barplot,plot_svtype_comparison
from utils.SV_TYPE import parse_pbsv_vcf
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


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

    logging.info("pbsv SV differential analysis finished.")
    logging.info(f"Results saved in: {outdir}")
    plot_svtype_comparison(
        type_summary,
        "/disk5/luosg/Totipotent20251031/PacBio/DEG/plot/sv_type.png",
        legend_map={
        "Control": "DMSO",
        "Experiment": "PlaB"
    })

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
