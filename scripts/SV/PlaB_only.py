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


def run_PlaB_only_annotation(
    dmso_vcf:str, 
    plab_vcf:str, 
    work_dir:str, 
    vep_cache:str = "~/.vep",
    species:str = "mus_musculus",
    assembly:str ="GRCm39",
    dist:int = 500,
    min_support:int = 1
):
    """
    参数化的 PlaB 独有变异注释流程
    :param dmso_vcf: DMSO 样本 VCF 路径 (Control)
    :param plab_vcf: PlaB 样本 VCF 路径 (Treatment)
    :param work_dir: 分析结果输出目录
    :param vep_cache: VEP 缓存目录
    :param species: 物种名称
    :param assembly: 基因组版本
    :param dist: SURVIVOR 合并距离阈值
    :param min_support: SURVIVOR 最小支持样本数
    """
    
    work_path = Path(work_dir).absolute()
    work_path.mkdir(parents=True, exist_ok=True)
    
    analysis = VEP_SV(
        vep_cache_dir=vep_cache,
        species=species,
        assembly=assembly
    )
    

    merged_vcf = work_path / "merged_sv.vcf"
    specific_vcf = work_path / "PlaB_only.vcf"
    annotated_vcf = work_path / "PlaB_only_annotated.vcf"

    logger.info(">>> Starting Pipeline: Merge -> Extract -> Annotate")
    
    # 合并 (注意列表顺序: DMSO 在前, PlaB 在后，对应 SUPP_VEC='01')
    raw_inputs = [dmso_vcf, plab_vcf]
    analysis.merge_sv_survivor(raw_inputs, str(merged_vcf), dist=dist, min_support=min_support)
    
    # 提取 PlaB 独有 (Vector 为 01)
    analysis.extract_specific_sv(str(merged_vcf), str(specific_vcf), vec="01")
    
    analysis.annotate_sv_vep(str(specific_vcf), str(annotated_vcf))
    
    logger.info(f">>> Pipeline Finished. Annotated VCF: {annotated_vcf}")
    
    return str(annotated_vcf)

def run_PlaB_specific_vcf_analysis(
        vcf:str,
        outdir:str
):
    """
    Docstring for run_PlaB_specific_vcf_analysis
    
    :param vcf: PlaB specific SV VCF path
    :param outdir: Output directory for analysis results
    """
    outdir = Path(outdir)
    outdir_table = outdir / "table"
    outdir_table.mkdir(parents=True,exist_ok=True)
    outdir_plot = outdir / "plot"
    outdir_plot.mkdir(parents=True,exist_ok=True)
    logger.info("Starting PlaB specific VCF analysis and visualization...")
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

def run_PlaB_enricher(
        vcf_01:str,
        vcf_1x:str,
        outdir:str,
        min_svlen:int = 2000,
        max_svlen:int = 10000
):
    """
    Docstring for run_PlaB_enricher
    (DMSO,PlaB),01: only PlaB have, 1x DMSO have
    :param vcf_01: PlaB only SV VCF path
    :param vcf_1x: DMSO SV VCF path
    :param outdir: Output directory for enrichment analysis results
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
    paraser.add_argument("-c","--dmso_vcf",type=str,required=True,help="DMSO sample VCF path")
    paraser.add_argument("-e","--plab_vcf",type=str,required=True,help="PlaB sample VCF path")
    paraser.add_argument("-o","--work_dir",type=str,required=True,help="Working directory for intermediate and final results")
    paraser.add_argument("-d","--dist",type=int,default=500,help="SURVIVOR merge distance threshold")
    args = paraser.parse_args()
    annotated_vcf = run_PlaB_only_annotation(
        dmso_vcf=args.dmso_vcf,
        plab_vcf=args.plab_vcf,
        work_dir=args.work_dir,
        dist=args.dist
    )
    run_PlaB_specific_vcf_analysis(
        vcf=annotated_vcf,
        outdir=f"{args.work_dir}/PlaB"
    )
    run_PlaB_enricher(
        vcf_01=f"{args.work_dir}/PlaB_only.vcf",
        vcf_1x=args.dmso_vcf,
        outdir=f"{args.work_dir}/PlaB_enrichment"
    )

if __name__ == "__main__":
    main()
    # config = {
    #     "dmso_vcf": "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/DMSO.sv.vcf.gz",
    #     "plab_vcf": "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/PlaB.sv.vcf.gz",
    #     "work_dir": "/disk5/luosg/Totipotent20251031/PacBio/SV",
    #     "dist": 500 
    # }

    # final_output = run_PlaB_only_annotation(**config)

    # config = {
    #     "vcf": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB_only.vcf",
    #     "outdir": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB"
    # }
    # run_PlaB_specific_vcf_analysis(**config)

    # config = {
    #     "vcf_01": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB_only.vcf",
    #     "vcf_1x": "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/DMSO.sv.vcf.gz",
    #     "outdir": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB"
    # }
    # run_PlaB_enricher(**config)