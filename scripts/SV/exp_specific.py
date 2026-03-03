import os

from utils.VEP_SV import VEP_SV
from utils.SV_TYPE import parse_pbsv_vcf,run_sv_stratification,extract_te_candidate_ins,generate_plot_input
from utils.SV_TYPE_plot import plot_stacking_bar,plot_multi_smooth_curves
from utils.repeatmasker_analysis import run_te_annotation_pipeline,RepeatMaskerOutCompare
from utils.repeatmasker_plot import plot_enrichment
from pathlib import Path
import logging
import pandas as pd
import argparse
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


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
    """
    Function: Run a pipeline to identify and annotate treatment-specific SVs using VEP.
    Steps:
    1. Merge control and experiment VCFs using SURVIVOR with specified distance and support thresholds.
    2. Extract SVs specific to the experiment (SUPP_VEC='01') using bcftools.
    3. Annotate the extracted SVs with VEP, using the provided cache, species, and assembly parameters. Output format can be VCF or tab
    Parameters:
    - ctrl_vcf: path to control VCF (e.g. DMSO)
    - exp_vcf: path to experiment VCF (e.g. PlaB)
    - outprefix: prefix for output files
    - vep_cache: path to VEP cache directory
    - species: species name for VEP annotation (e.g. "mus_musculus")
    - assembly: genome assembly for VEP annotation (e.g. "GRCm39")
    - dist: distance threshold for SURVIVOR merging (default 500 bp)
    - min_support: minimum support for SURVIVOR merging (default 1)
    - annotate_format: output format for VEP annotation ("vcf" or "tab", default "vcf")
    Returns:
    - Path to the annotated VCF or tab file containing experiment-specific SVs with VEP annotations.

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
        outdir:str
):
    """
    Function: Analyze the annotated VCF to generate summary statistics and visualizations.
    Steps:
    1. Parse the annotated VCF to extract SV type and size information.
    2. Generate a size distribution plot (stacked bar) for different SV types.
    3. Generate a curve plot showing SV count across size bins, highlighting specific size thresholds (e.g. 6kb).
    Parameters:
    - vcf: path to the annotated VCF file containing experiment-specific SVs with VEP annotations.
    - outdir: directory where output tables and plots will be saved. Subdirectories "table" and "plot" will be created within this directory.
    Returns:
    - None. Side effects include:
        - A size distribution plot saved as "size_bin.png" in the "plot" subdirectory.
        - A curve plot of SV count across size bins saved as "svlen_count.png" in the "plot" subdirectory.
        - Intermediate tables used for plotting saved in the "table" subdirectory.

    """
    outdir = Path(outdir)
    outdir_table = outdir / "table"
    outdir_table.mkdir(parents=True,exist_ok=True)
    outdir_plot = outdir / "plot"
    outdir_plot.mkdir(parents=True,exist_ok=True)
    logger.info("Starting VCF analysis and visualization...")
    matrix,summary = run_sv_stratification(vcf,str(outdir_table))
    outpng_bar = outdir_plot / "size_bin.png"
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
    outpng_curve = outdir_plot / "svlen_count.png"
    plot_multi_smooth_curves(plot_data,
                             x_label="svlen",
                             y_label="count",
                             title="",
                             outfig=outpng_curve,
                             highlight_x_values=[6000])

def run_exp_enricher(
        vcf_01:str,
        vcf_1x:str,
        outdir:str,
        min_svlen:int = 2000,
        max_svlen:int = 10000
):
    """
    Function: Perform enrichment analysis of repeat elements in experiment-specific SV insertions compared to control SV insertions.
    Steps:
    1. Extract candidate TE insertions from both experiment-specific (SUPP_VEC='01') and control (SUPP_VEC='1x') VCFs, filtering by specified SV length range.
    2. Annotate the extracted candidate insertions using RepeatMasker to identify repeat element composition.
    3. Perform enrichment analysis comparing the repeat element composition of experiment-specific insertions to control insertions, focusing on subfamily level and applying a divergence filter (e.g. max_div=3).
    4. Generate a plot visualizing significantly enriched subfamilies (FDR < 0.05) in the experiment-specific insertions.
    Parameters:
    - vcf_01: path to the VCF file containing experiment-specific SVs (SUPP_VEC='01').
    - vcf_1x: path to the VCF file containing control SVs (SUPP_VEC='1x').
    - outdir: directory where output tables and plots will be saved. Subdirectories "fa" and "repeatmasker" will be created within this directory for intermediate files.
    - min_svlen: minimum SV length (in bp) for candidate insertion extraction (default 2000 bp).
    - max_svlen: maximum SV length (in bp) for candidate insertion extraction (default 10000 bp).
    Returns:
    - None. Side effects include:
        - FASTA files of candidate insertions saved in the "fa" subdirectory.
        - RepeatMasker annotation outputs saved in the "repeatmasker" subdirectory.
        - An enrichment results CSV file saved as "PlaBOnlyEnrichment.csv" in the "repeatmasker" subdirectory.
        - An enrichment plot saved as "PlaBOnlyEnrichment.png" in the "repeatmasker" subdirectory.

    """
    logger.info("Starting PlaB specific insertion enrichment analysis...")
    outdir = Path(outdir)
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
    repeatMaskerOutCompare = RepeatMaskerOutCompare(bg_out=out_1x,fg_out=out_01)
    df = repeatMaskerOutCompare.enrichment_test(level="subfamily",max_div=3)
    out_enrich_csv = f"{outdir}/repeatmasker/PlaBOnlyEnrichment.csv"
    df.to_csv(out_enrich_csv,sep="\t",index=False)
    logger.info(f"Saving enrichment results to: {out_enrich_csv}")

    logger.info("Generating enrichment plot for significant subfamilies (FDR < 0.05)...")
    df = df[df["fdr"] < 0.05]
    out_enrich_png = f"{outdir}/repeatmasker/PlaBOnlyEnrichment.png"
    plot_enrichment(df,out_enrich_png)
    logger.info(f"Generating enrichment plot: {out_enrich_png}")

def main():
    paraser = argparse.ArgumentParser(description="PlaB specific SV analysis pipeline")
    paraser.add_argument("-c","--ctrl_vcf",type=str,required=True,help="control sample VCF path")
    paraser.add_argument("-e","--exp_vcf",type=str,required=True,help="experiment sample VCF path")
    paraser.add_argument("-o","--outprefix",type=str,required=True,help="out prefix for merged, specific, and annotated files")
    paraser.add_argument("-d","--dist",type=int,default=500,help="SURVIVOR merge distance threshold")
    paraser.add_argument("--vep_cache",type=str,default="~/.vep",help="VEP cache directory")
    paraser.add_argument("--species",type=str,default="mus_musculus",help="species name for VEP annotation")
    paraser.add_argument("--assembly",type=str,default="GRCm39",help="genome assembly for VEP annotation")
    paraser.add_argument("--annotate_format",type=str,choices=["vcf","tab"],default="tab",help="output format for VEP annotation")
    args = paraser.parse_args()

    annotated_file = run_exp_specific_annotation(
        ctrl_vcf=args.ctrl_vcf,
        exp_vcf=args.exp_vcf,
        outprefix=args.outprefix,
        dist=args.dist,
        annotate_format=args.annotate_format
    )

    outdir = os.path.dirname(args.outprefix)
    run_vcf_analysis(
        vcf=f"{args.outprefix}_only.vcf",
        outdir=outdir
    )
    run_exp_enricher(
        vcf_01=f"{args.outprefix}_only.vcf",
        vcf_1x=args.ctrl_vcf,
        outdir=f"{outdir}/enrichment"
    )

if __name__ == "__main__":
    main()
