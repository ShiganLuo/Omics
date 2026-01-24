from utils.VEP_SV import VEP_SV
from utils.SV_TYPE import parse_pbsv_vcf,run_sv_stratification,extract_te_candidate_ins,generate_plot_input
from utils.SV_TYPE_plot import plot_stacking_bar,plot_multi_smooth_curves
from utils.repeatmasker_analysis import run_te_annotation_pipeline,RepeatMaskerOutCompare
from utils.repeatmasker_plot import plot_enrichment
from pathlib import Path
import logging
import pandas as pd
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
logger = logging.getLogger(__name__)


def run_PlaB_only_annotation(
    dmso_vcf, 
    plab_vcf, 
    work_dir, 
    vep_cache="/home/luosg/.vep",
    species="mus_musculus",
    assembly="GRCm39",
    dist=500,
    min_support=1
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
    
    analysis = VEP_SV(
        vep_cache_dir=vep_cache,
        species=species,
        assembly=assembly
    )
    
    work_path = Path(work_dir).absolute()
    work_path.mkdir(parents=True, exist_ok=True)
    
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
    outdir = Path(outdir)
    outdir_table = outdir / "table"

    ######### stacking bar
    matrix,summary = run_sv_stratification(vcf,str(outdir_table))
    outpng = outdir / "plot/size_bin.png"
    plot_stacking_bar(matrix,xlabel="bin",ylabel="count",title="",legend_title_type="sv type",legend_width=0.15,save_path=outpng)

    ##### svlen count
    df_sta = parse_pbsv_vcf(vcf)
    plot_data = generate_plot_input(df_sta)
    outfig = outdir / "plot/svlen_count.png"
    print(plot_data)
    plot_multi_smooth_curves(plot_data,x_label="svlen",y_label="count",title="",outfig=outfig,highlight_x_values=[6000])

def run_PlaB_enricher(
        vcf_01:str,
        vcf_1x:str,
        outdir:str
):
    """
    Docstring for run_PlaB_enricher
    (DMSO,PlaB),01: only PlaB have, 1x DMSO have
    """
    # fa_01 = f"{outdir}/fa/01/INS_2-10kb_01.fa"
    # extract_te_candidate_ins(vcf_01,fa_01)
    # out_01 = run_te_annotation_pipeline(fa_01,f"{outdir}/repeatmasker/01",species="mus musculus")
    # fa_1x = f"{outdir}/fa/1x/INS_2-10kb_1x.fa"
    # extract_te_candidate_ins(vcf_1x,fa_1x)
    # out_1x = run_te_annotation_pipeline(fa_1x,f"{outdir}/repeatmasker/1x",species="mus musculus")

    # repeatMaskerOutCompare = RepeatMaskerOutCompare(bg_out=out_1x,fg_out=out_01)
    repeatMaskerOutCompare = RepeatMaskerOutCompare(bg_out="/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB/repeatmasker/1x/INS_2-10kb_1x.fa.out",
                                                    fg_out="/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB/repeatmasker/01/INS_2-10kb_01.fa.out")
    df = repeatMaskerOutCompare.enrichment_test(level="subfamily",max_div=3)
    df.to_csv("/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB/repeatmasker/PlaBOnlyEnrichment.csv",sep="\t",index=False)
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB/repeatmasker/PlaBOnlyEnrichment.csv",sep="\t")
    df = df[df["fdr"] < 0.05]
    plot_enrichment(df,"/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB/repeatmasker/PlaBOnlyEnrichment.png")


if __name__ == "__main__":
    config = {
        "dmso_vcf": "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/DMSO.sv.vcf.gz",
        "plab_vcf": "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/PlaB.sv.vcf.gz",
        "work_dir": "/disk5/luosg/Totipotent20251031/PacBio/SV",
        "dist": 500 
    }

    # final_output = run_PlaB_only_annotation(**config)

    # config = {
    #     "vcf": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB_only.vcf",
    #     "outdir": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB"
    # }
    # run_PlaB_specific_vcf_analysis(**config)

    config = {
        "vcf_01": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB_only.vcf",
        "vcf_1x": "/disk5/luosg/Totipotent20251031/data/Pacbio/unphased/DMSO.sv.vcf.gz",
        "outdir": "/disk5/luosg/Totipotent20251031/PacBio/SV/PlaB"
    }
    run_PlaB_enricher(**config)