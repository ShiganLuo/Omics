import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.LogUtil import setup_logger
import pandas as pd
import argparse
import os
import logging
from typing import Dict, List, Literal, Optional
from utils.SV_TYPE_plot import plot_sv_type_barplot,plot_sv_length_boxplot,plot_large_sv_barplot,plot_svtype_comparison
from utils.SV_TYPE import parse_pbsv_vcf
logger = setup_logger("pbsvDiffAnalysis", level=logging.INFO)

PlotFormat = Literal["png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf", "raw", "rgba"]


def run_pbsv_diff_analysis(
    group_vcf: Dict[str, str],
    out_dir: str,
    large_sv_threshold: int,
    plot_formats: Optional[List[PlotFormat]] = None,
):
    """Run SV differential analysis across multiple sample groups.

    Parse VCF files for each group, compute SV type and length summaries,
    and generate comparison plots (bar plots, box plots).

    Parameters
    ----------
    group_vcf : dict of str to str
        Mapping of group name to VCF file path.
    out_dir : str
        Output directory for tables and plots.
    large_sv_threshold : int
        Minimum SV length (bp) to classify as a large SV.
    plot_formats : list of str, optional
        Image formats for output plots. Defaults to ``["png"]``.

    Returns
    -------
    None
    """
    if plot_formats is None:
        plot_formats = ["png"]

    os.makedirs(f"{out_dir}/table", exist_ok=True)
    os.makedirs(f"{out_dir}/plot", exist_ok=True)
    df_all = {}
    for group, vcf in group_vcf.items():
        df = parse_pbsv_vcf(vcf)
        df["group"] = group
        df_all[group] = df
    all_df = pd.concat(df_all.values(), ignore_index=True)

    # 保存所有 SV
    all_df.to_csv(f"{out_dir}/table/all_sv_records.tsv", sep="\t", index=False)

    # SV 类型统计
    type_summary = (
        all_df.groupby(["group", "svtype"])
        .size()
        .reset_index(name="count")
    )

    type_summary.to_csv(
        f"{out_dir}/table/sv_type_summary.tsv", sep="\t", index=False
    )

    # 图 1：SV 类型
    for fmt in plot_formats:
        plot_sv_type_barplot(
            summary_df=type_summary,
            outpng = f"{out_dir}/plot/sv_type_barplot.{fmt}",
            xlabel="SV type",
            ylabel="SV count",
        )

    # 图 2：SV 长度分布
    for fmt in plot_formats:
        plot_sv_length_boxplot(
            data_dict=df_all,
            svlen_col="svlen",
            outpng = f"{out_dir}/plot/large_sv_{large_sv_threshold}bp_boxplot.{fmt}",
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
        f"{out_dir}/table/large_sv_{large_sv_threshold}bp_summary.tsv",
        sep="\t",
        index=False,
    )

    if not large_summary.empty:
        for fmt in plot_formats:
            plot_large_sv_barplot(
                summary_df=large_summary,
                outpng = f"{out_dir}/plot/large_sv_{large_sv_threshold}bp_barplot.{fmt}",
                size_threshold=large_sv_threshold,
                title=f"Large SV comparison (≥{large_sv_threshold} bp)",
                ylabel=f"SV count (≥{large_sv_threshold} bp)",
            )

    logging.info("pbsv SV differential analysis finished.")
    logging.info(f"Results saved in: {out_dir}")
    for fmt in plot_formats:
        out_svtype_image = f"{out_dir}/plot/sv_type_barplot_comparison.{fmt}"
        logging.info(f"SV type distribution plot: {out_svtype_image}")
        plot_svtype_comparison(
            type_summary,
            out_svtype_image,
            group_order=list(group_vcf.keys()),
        )

def specific_svtype_diff_analysis():
    """Run SV differential analysis on hardcoded comparison groups.

    Compares DMSO and PlaB treatment groups at two time points using
    pre-defined VCF paths and a 10 kb large-SV threshold.

    Returns
    -------
    None
    """
    DMSO20_vs_DMSO06 = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO20_vs_DMSO06/DMSO_only.vcf"
    PlaB20_vs_PlaB06 = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_PlaB06/PlaB_only.vcf"
    PlaB20_vs_DMSO20 = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB_only.vcf"
    PlaB06_vs_DMSO06 = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_only.vcf"
    group_vcf = {
        "DMSO20_vs_DMSO06": DMSO20_vs_DMSO06,
        "PlaB20_vs_PlaB06": PlaB20_vs_PlaB06,
        "PlaB20_vs_DMSO20": PlaB20_vs_DMSO20,
        "PlaB06_vs_DMSO06": PlaB06_vs_DMSO06,
    }
    run_pbsv_diff_analysis(
        group_vcf=group_vcf,
        out_dir="/data/pub/zhousha/Totipotent20251031/PacBio/SV/pbsv_specific_svtype_diff",
        large_sv_threshold=10_000,
    )

def main():
    """Entry point for command-line execution of pbsv SV differential analysis.

    Parse command-line arguments for group VCFs, output directory,
    large-SV threshold, and plot formats, then invoke
    :func:`run_pbsv_diff_analysis`.

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(
        description="pbsv SV differential analysis (VCF / VCF.GZ)"
    )
    parser.add_argument(
        "-g", "--group",
        action="append",
        required=True,
        metavar="NAME:VCF",
        help="Group in format 'name:vcf_path'. Can be specified multiple times. "
             "Example: -g Control:ctrl.vcf -g Experiment:exp.vcf",
    )
    parser.add_argument("-o", "--outdir", default="pbsv_sv_diff")
    parser.add_argument("-s", "--large-sv-threshold", type=int, default=10_000)
    parser.add_argument(
        "-f", "--format",
        action="append",
        dest="formats",
        metavar="FMT",
        help="Image output format (png, pdf, svg, tiff, ...). "
             "Can be specified multiple times. Default: png.",
    )

    args = parser.parse_args()

    group_vcf = {}
    for g in args.group:
        if ":" not in g:
            parser.error(f"Invalid format '{g}', expected 'NAME:VCF'")
        name, vcf = g.split(":", 1)
        group_vcf[name] = vcf

    run_pbsv_diff_analysis(
        group_vcf=group_vcf,
        out_dir=args.outdir,
        large_sv_threshold=args.large_sv_threshold,
        plot_formats=args.formats,
    )


if __name__ == "__main__":
    main()
