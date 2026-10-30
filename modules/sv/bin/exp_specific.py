import os
import sys
from typing import Literal, List, Optional
from pathlib import Path
import logging
import argparse
try:
    from .utils.common import setup_logger
    from .utils.VEP_SV import VEP_SV
    from .utils.SV_TYPE import parse_pbsv_vcf,run_sv_stratification,extract_te_candidate_ins,generate_plot_input
    from .utils.SV_TYPE_plot import plot_stacking_bar,plot_multi_smooth_curves
    from .utils.repeatmasker_analysis import run_te_annotation_pipeline,RepeatMaskerOutCompare
    from .utils.repeatmasker_plot import plot_enrichment
except ImportError:
    from utils.common import setup_logger
    from utils.VEP_SV import VEP_SV
    from utils.SV_TYPE import parse_pbsv_vcf,run_sv_stratification,extract_te_candidate_ins,generate_plot_input
    from utils.SV_TYPE_plot import plot_stacking_bar,plot_multi_smooth_curves
    from utils.repeatmasker_analysis import run_te_annotation_pipeline,RepeatMaskerOutCompare
    from utils.repeatmasker_plot import plot_enrichment
logger = setup_logger("exp_specific", level=logging.INFO)
PlotFormat = Literal["png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf", "raw", "rgba"]


def run_exp_specific_annotation(
    ctrl_vcfs:List[str],
    exp_vcfs:List[str],
    outprefix:str,
    vep_cache: str = "~/.vep",
    species: str = "mus_musculus",
    assembly: str = "GRCm39",
    custom_gtf: Optional[str] = None,
    custom_fasta: Optional[str] = None,
    dist:int = 500,
    min_support:int = 1,
    annotate_format:str = "vcf"
):
    """Identify and annotate experiment-specific SVs using a VEP-based pipeline.

    Merges control and experiment VCFs via SURVIVOR (controls first,
    experiments second), extracts SVs supported by EVERY experiment sample
    and no control sample (SUPP_VEC = ``"0" * n_ctrl + "1" * n_exp``), and
    annotates them with VEP. Output format can be VCF or tab-delimited.

    Parameters
    ----------
    ctrl_vcfs : list of str
        Control sample VCF paths (one per sample; replicate groups allowed).
    exp_vcfs : list of str
        Experiment sample VCF paths (one per sample; replicate groups allowed).
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
    tuple of (str, str)
        Path to the annotated file containing experiment-specific SVs with
        VEP annotations, and path to the control callset VCF (the TE
        enrichment background): the single control VCF when one control is
        given, otherwise a SURVIVOR merge of the control VCFs
        (``<outprefix>_control_merged.vcf``).
    """

    outdir = os.path.dirname(outprefix)
    os.makedirs(outdir, exist_ok=True)

    analysis = VEP_SV(
        vep_cache_dir=vep_cache,
        species=species,
        assembly=assembly,
        custom_gtf=custom_gtf,
        custom_fasta=custom_fasta
    )


    merged_vcf = f"{outprefix}_merged.vcf"
    specific_vcf = f"{outprefix}_only.vcf"
    annotated_file = f"{outprefix}_annotated.{annotate_format}"

    logger.info(">>> Starting Pipeline: Merge -> Extract -> Annotate")

    # Merge controls first, experiments second: SUPP_VEC slots follow this order
    raw_inputs = list(ctrl_vcfs) + list(exp_vcfs)
    analysis.merge_sv_survivor(raw_inputs, merged_vcf, dist=dist, min_support=min_support)

    # Experiment-specific SVs: present in every experiment sample and
    # absent from every control sample (SUPP_VEC = 00..11).
    vec = "0" * len(ctrl_vcfs) + "1" * len(exp_vcfs)
    analysis.extract_specific_sv(merged_vcf, specific_vcf, vec=vec)

    analysis.annotate_sv_vep(specific_vcf, annotated_file, result_format=annotate_format)

    # Control callset used as the TE-enrichment background
    if len(ctrl_vcfs) == 1:
        control_callset = ctrl_vcfs[0]
    else:
        control_callset = f"{outprefix}_control_merged.vcf"
        analysis.merge_sv_survivor(list(ctrl_vcfs), control_callset, dist=dist, min_support=min_support)

    logger.info(f">>> Pipeline Finished. Annotated VCF: {annotated_file}")

    return annotated_file, control_callset

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
        te_lib: Optional[str] = None,
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
    out_01 = run_te_annotation_pipeline(fa_01,f"{outdir}/repeatmasker/01",species="mus musculus",lib=te_lib)
    fa_1x = f"{outdir}/fa/1x/INS_{min_svlen/1000}-{max_svlen/1000}_1x.fa"
    extract_te_candidate_ins(vcf_1x,fa_1x,min_len=min_svlen,max_len=max_svlen)
    out_1x = run_te_annotation_pipeline(fa_1x,f"{outdir}/repeatmasker/1x",species="mus musculus",lib=te_lib)
    if out_01 is None or out_1x is None:
        raise RuntimeError(
            "RepeatMasker TE annotation failed (see log above); "
            "cannot run the enrichment test"
        )

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
    paraser.add_argument("-c","--ctrl_vcf",type=str,nargs="+",required=True,help="control sample VCF path(s)")
    paraser.add_argument("-e","--exp_vcf",type=str,nargs="+",required=True,help="experiment sample VCF path(s)")
    paraser.add_argument("-o","--outprefix",type=str,required=True,help="out prefix for merged, specific, and annotated files")
    paraser.add_argument("-d","--dist",type=int,default=500,help="SURVIVOR merge distance threshold")
    paraser.add_argument("--vep_cache",type=str,default="~/.vep",help="VEP cache directory")
    paraser.add_argument("--species",type=str,default="mus_musculus",help="species name for VEP annotation")
    paraser.add_argument("--assembly",type=str,default="GRCm39",help="genome assembly for VEP annotation")
    paraser.add_argument("--annotate_format",type=str,choices=["vcf","tab"],default="tab",help="output format for VEP annotation")
    paraser.add_argument("--custom_gtf",type=str,default=None,help="sorted bgzipped tabix-indexed GTF/GFF for custom genomes (bypasses the Ensembl cache)")
    paraser.add_argument("--custom_fasta",type=str,default=None,help="genome FASTA matching --custom_gtf (required with it)")
    paraser.add_argument("--te_lib",type=str,default=None,help="RepeatMasker library FASTA for TE annotation (for species whose Dfam partition is absent from the container, e.g. mouse)")
    paraser.add_argument("-f","--format",action="append",dest="formats",metavar="FMT",help="Image output format (png, pdf, svg, ...). Can be specified multiple times. Default: png.")
    args = paraser.parse_args()

    annotated_file, control_callset = run_exp_specific_annotation(
        ctrl_vcfs=args.ctrl_vcf,
        exp_vcfs=args.exp_vcf,
        outprefix=args.outprefix,
        vep_cache=args.vep_cache,
        species=args.species,
        assembly=args.assembly,
        custom_gtf=args.custom_gtf,
        custom_fasta=args.custom_fasta,
        dist=args.dist,
        annotate_format=args.annotate_format,
    )

    outdir = os.path.dirname(args.outprefix)
    run_vcf_analysis(
        vcf=f"{args.outprefix}_only.vcf",
        out_dir=outdir,
        image_formats=args.formats
    )
    # NOTE: TE annotation writes to {outdir}/repeatmasker, deliberately NOT
    # under {outdir}/enrichment — sv_exp_enrichment declares
    # directory(enrichment/) as its Snakemake output, and Snakemake removes
    # the whole directory before that rule runs, which would silently wipe
    # these results (observed 2026-10-29).
    run_exp_enricher(
        vcf_01=f"{args.outprefix}_only.vcf",
        vcf_1x=control_callset,
        out_dir=f"{outdir}/repeatmasker",
        te_lib=args.te_lib,
        image_formats=args.formats
    )

if __name__ == "__main__":
    main()
