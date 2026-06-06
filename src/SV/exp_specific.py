import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.LogUtil import setup_logger
from utils.VEP_SV import VEP_SV
from utils.SV_TYPE import parse_pbsv_vcf,run_sv_stratification,extract_te_candidate_ins,generate_plot_input
from utils.SV_TYPE_plot import plot_stacking_bar,plot_multi_smooth_curves
from utils.repeatmasker_analysis import run_te_annotation_pipeline,RepeatMaskerOutCompare
from utils.repeatmasker_plot import plot_enrichment
from typing import Literal, List, Optional
PlotFormat = Literal["png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf", "raw", "rgba"]
from pathlib import Path
import logging
import argparse
logger = setup_logger("ExpSpecificSVPipeline", level=logging.INFO)


def run_exp_specific_annotation(
    ctrl_vcf:str, 
    exp_vcf:str, 
    outprefix:str, 
    vep_cache:str = "~/.vep",
    species:str = "mus_musculus",
    assembly:str ="GRCm39",
    dist:int = 500,
    min_support:int = 1,
    annotate_format:str = "vcf"
):
    """Identify and annotate treatment-specific SVs using a VEP-based pipeline.

    The pipeline merges control and experiment VCFs via SURVIVOR, extracts
    SVs specific to the experiment (SUPP_VEC='01'), and annotates them with
    VEP. Output format can be VCF or tab-delimited.

    Parameters
    ----------
    ctrl_vcf : str
        Path to the control VCF file (e.g. DMSO).
    exp_vcf : str
        Path to the experiment VCF file (e.g. PlaB).
    outprefix : str
        Prefix for all output files (merged, specific, annotated).
    vep_cache : str, optional
        Path to the VEP cache directory. Default is ``"~/.vep"``.
    species : str, optional
        Species name for VEP annotation. Default is ``"mus_musculus"``.
    assembly : str, optional
        Genome assembly for VEP annotation. Default is ``"GRCm39"``.
    dist : int, optional
        Distance threshold (bp) for SURVIVOR merging. Default is 500.
    min_support : int, optional
        Minimum sample support for SURVIVOR merging. Default is 1.
    annotate_format : str, optional
        Output format for VEP annotation (``"vcf"`` or ``"tab"``).
        Default is ``"vcf"``.

    Returns
    -------
    str
        Path to the annotated file containing experiment-specific SVs
        with VEP annotations.
    """
    
    outdir = os.path.dirname(outprefix)
    os.makedirs(outdir, exist_ok=True)
    
    analysis = VEP_SV(
        vep_cache_dir=vep_cache,
        species=species,
        assembly=assembly
    )
    

    merged_vcf = f"{outprefix}_merged.vcf"
    specific_vcf = f"{outprefix}_only.vcf"
    annotated_file = f"{outprefix}_annotated.{annotate_format}"

    logger.info(">>> Starting Pipeline: Merge -> Extract -> Annotate")
    
    # 合并 (注意列表顺序: DMSO 在前, PlaB 在后，对应 SUPP_VEC='01')
    raw_inputs = [ctrl_vcf, exp_vcf]
    analysis.merge_sv_survivor(raw_inputs, merged_vcf, dist=dist, min_support=min_support)
    
    # 提取 PlaB 独有 (Vector 为 01)
    analysis.extract_specific_sv(merged_vcf, specific_vcf, vec="01")
    
    analysis.annotate_sv_vep(specific_vcf, annotated_file, result_format=annotate_format)
    
    logger.info(f">>> Pipeline Finished. Annotated VCF: {annotated_file}")
    
    return annotated_file

def run_vcf_analysis(
        vcf:str,
        out_dir:str,
        image_formats:Optional[List[PlotFormat]] = None
):
    """Analyze an annotated VCF to produce summary statistics and visualizations.

    Parses the annotated VCF to extract SV type and size information, then
    generates a stacked-bar size distribution plot and a smoothed curve plot
    of SV counts across size bins.

    Parameters
    ----------
    vcf : str
        Path to the annotated VCF file containing experiment-specific SVs.
    out_dir : str
        Directory where output tables and plots will be saved.
        Subdirectories ``"table"`` and ``"plot"`` are created automatically.
    image_formats : list of str, optional
        Output image formats (e.g. ``["png", "pdf"]``). Default is
        ``["png"]``.

    Returns
    -------
    None

    Notes
    -----
    Side effects include:

    - A stacked-bar size distribution plot saved as ``size_bin.<fmt>``
      in the ``plot`` subdirectory.
    - A smoothed curve plot of SV count across size bins saved as
      ``svlen_count.<fmt>`` in the ``plot`` subdirectory.
    - Intermediate tables used for plotting saved in the ``table``
      subdirectory.
    """
    if image_formats is None:
        image_formats = ["png"]
    outdir = Path(out_dir)
    outdir_table = outdir / "table"
    outdir_table.mkdir(parents=True,exist_ok=True)
    outdir_plot = outdir / "plot"
    outdir_plot.mkdir(parents=True,exist_ok=True)
    logger.info("Starting VCF analysis and visualization...")
    matrix,summary = run_sv_stratification(vcf,str(outdir_table))
    for fmt in image_formats:
        outpng_bar = outdir_plot / f"size_bin.{fmt}"
        plot_stacking_bar(matrix,
                          xlabel="bin",
                          ylabel="count",
                          title="",
                          legend_title_type="sv type",
                          legend_width=0.15,
                          save_path=outpng_bar,
                          show_block_counts=True)

    # ##### svlen count
    df_sta = parse_pbsv_vcf(vcf)
    plot_data = generate_plot_input(df_sta)
    for fmt in image_formats:
        outpng_curve = outdir_plot / f"svlen_count.{fmt}"
        plot_multi_smooth_curves(plot_data,
                                 x_label="svlen",
                                 y_label="count",
                                 title="",
                                 outfig=str(outpng_curve),
                                 highlight_x_values=[6000])

def run_exp_enricher(
        vcf_01:str,
        vcf_1x:str,
        out_dir:str,
        min_svlen:int = 2000,
        max_svlen:int = 10000,
        image_formats:Optional[List[PlotFormat]] = None
):
    """Perform repeat-element enrichment analysis of experiment-specific SV insertions.

    Extracts candidate transposable-element (TE) insertions from both
    experiment-specific and control VCFs, annotates them with RepeatMasker,
    and tests for enrichment of repeat subfamilies in the experiment set
    relative to controls.

    Parameters
    ----------
    vcf_01 : str
        Path to the VCF file containing experiment-specific SVs
        (SUPP_VEC='01').
    vcf_1x : str
        Path to the VCF file containing control SVs (SUPP_VEC='1x').
    out_dir : str
        Directory where output tables and plots will be saved.
        Subdirectories ``"fa"`` and ``"repeatmasker"`` are created
        automatically for intermediate files.
    min_svlen : int, optional
        Minimum SV length (bp) for candidate insertion extraction.
        Default is 2000.
    max_svlen : int, optional
        Maximum SV length (bp) for candidate insertion extraction.
        Default is 10000.
    image_formats : list of str, optional
        Output image formats (e.g. ``["png", "pdf"]``). Default is
        ``["png"]``.

    Returns
    -------
    None

    Notes
    -----
    Side effects include:

    - FASTA files of candidate insertions in the ``fa`` subdirectory.
    - RepeatMasker annotation outputs in the ``repeatmasker`` subdirectory.
    - An enrichment results CSV file saved as ``Enrichment.csv``
      in the ``repeatmasker`` subdirectory.
    - An enrichment plot for significant subfamilies (FDR < 0.05) saved
      as ``Enrichment.<fmt>`` in the ``repeatmasker`` subdirectory.
    """
    if image_formats is None:
        image_formats = ["png"]
    logger.info("Starting PlaB specific insertion enrichment analysis...")
    outdir = Path(out_dir)
    outdir.mkdir(parents=True,exist_ok=True)
    logger.info(f"Output directory: {outdir}")

    logger.info(f"Extracting TE candidate insertions from VCFs...: {vcf_01}, {vcf_1x}")
    fa_01 = f"{outdir}/fa/01/INS_{min_svlen/1000}-{max_svlen/1000}kb_01.fa"
    extract_te_candidate_ins(vcf_01,fa_01,min_len=min_svlen,max_len=max_svlen)
    out_01 = run_te_annotation_pipeline(fa_01,f"{outdir}/repeatmasker/01",species="mus musculus")
    fa_1x = f"{outdir}/fa/1x/INS_{min_svlen/1000}-{max_svlen/1000}_1x.fa"
    extract_te_candidate_ins(vcf_1x,fa_1x,min_len=min_svlen,max_len=max_svlen)
    out_1x = run_te_annotation_pipeline(fa_1x,f"{outdir}/repeatmasker/1x",species="mus musculus")

    logger.info("Performing enrichment test between PlaB only and DMSO SV insertions...")
    repeatMaskerOutCompare = RepeatMaskerOutCompare(bg_out=str(out_1x),fg_out=str(out_01))
    df = repeatMaskerOutCompare.enrichment_test(level="subfamily",max_div=3)
    out_enrich_csv = f"{outdir}/repeatmasker/Enrichment.csv"
    df.to_csv(out_enrich_csv,sep="\t",index=False)
    logger.info(f"Saving enrichment results to: {out_enrich_csv}")

    logger.info("Generating enrichment plot for significant subfamilies (FDR < 0.05)...")
    df = df[df["fdr"] < 0.05]
    for fmt in image_formats:
        out_enrich_png = f"{outdir}/repeatmasker/Enrichment.{fmt}"
        plot_enrichment(df,out_enrich_png)
        logger.info(f"Generating enrichment plot: {out_enrich_png}")

def main():
    """Entry point for the experiment-specific SV analysis pipeline.

    Parses command-line arguments and orchestrates the full pipeline:
    merging VCFs, extracting experiment-specific SVs, annotating with VEP,
    generating summary plots, and performing repeat-element enrichment.

    Returns
    -------
    None
    """
    paraser = argparse.ArgumentParser(description="PlaB specific SV analysis pipeline")
    paraser.add_argument("-c","--ctrl_vcf",type=str,required=True,help="control sample VCF path")
    paraser.add_argument("-e","--exp_vcf",type=str,required=True,help="experiment sample VCF path")
    paraser.add_argument("-o","--outprefix",type=str,required=True,help="out prefix for merged, specific, and annotated files")
    paraser.add_argument("-d","--dist",type=int,default=500,help="SURVIVOR merge distance threshold")
    paraser.add_argument("--vep_cache",type=str,default="~/.vep",help="VEP cache directory")
    paraser.add_argument("--species",type=str,default="mus_musculus",help="species name for VEP annotation")
    paraser.add_argument("--assembly",type=str,default="GRCm39",help="genome assembly for VEP annotation")
    paraser.add_argument("--annotate_format",type=str,choices=["vcf","tab"],default="tab",help="output format for VEP annotation")
    paraser.add_argument("-f","--format",action="append",dest="formats",metavar="FMT",help="Image output format (png, pdf, svg, ...). Can be specified multiple times. Default: png.")
    args = paraser.parse_args()

    annotated_file = run_exp_specific_annotation(
        ctrl_vcf=args.ctrl_vcf,
        exp_vcf=args.exp_vcf,
        outprefix=args.outprefix,
        dist=args.dist,
        annotate_format=args.annotate_format,
    )

    outdir = os.path.dirname(args.outprefix)
    run_vcf_analysis(
        vcf=f"{args.outprefix}_only.vcf",
        out_dir=outdir,
        image_formats=args.formats
    )
    run_exp_enricher(
        vcf_01=f"{args.outprefix}_only.vcf",
        vcf_1x=args.ctrl_vcf,
        out_dir=f"{outdir}/enrichment",
        image_formats=args.formats
    )

if __name__ == "__main__":
    main()
