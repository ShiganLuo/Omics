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
    Overall purpose:
    - Merge structural variant (SV) calls from two samples (control and treatment),
        extract variants that are specific to the treatment (PlaB), and annotate those
        treatment-specific SVs using VEP. The function writes intermediate files into
        the provided working directory and returns the path to the annotated VCF.

    Use cases:
    - Comparative SV analysis between a control (e.g. DMSO) and a treatment (PlaB)
        to discover treatment-specific structural variants followed by functional annotation.

    Data flow (input -> processing -> output):
    - Input: `dmso_vcf`, `plab_vcf`, and parameters controlling merge/annotation behavior.
    - Processing: call SURVIVOR merge via `analysis.merge_sv_survivor` -> extract treatment-only
        entries by support vector -> annotate extracted VCF via VEP (`analysis.annotate_sv_vep`).
    - Output: annotated VCF file path string.

    Parameters:
    - dmso_vcf (str): required. Path to the control sample VCF (.vcf or .vcf.gz). Must be readable.
        Constraint: VCF should contain fields compatible with SURVIVOR/pbsv (e.g. SUPP_VEC) if used upstream.
    - plab_vcf (str): required. Path to the treatment (PlaB) sample VCF. Same constraints as `dmso_vcf`.
    - work_dir (str): required. Directory where intermediate and result files are written. The directory
        will be created if not present. Caller must ensure write permission and sufficient disk space.
    - vep_cache (str): optional (default "~/.vep"). Path to VEP cache directory used to accelerate VEP annotation.
        Constraint: if provided, it should be accessible and contain the expected VEP cache files.
    - species (str): required (default "mus_musculus"). Species identifier for VEP.
    - assembly (str): required (default "GRCm39"). Genome assembly identifier for VEP.
    - dist (int): required (default 500). Merge distance threshold for SURVIVOR (in base pairs). Must be >= 0.
    - min_support (int): required (default 1). Minimum support samples for SURVIVOR merging. Must be >= 1.

    Returns:
    - str: absolute path to the annotated VCF file (e.g. /path/to/PlaB_only_annotated.vcf).
        The file is a standard VCF and may include VEP annotation fields (e.g. CSQ).
        On failure the function raises an exception and does not return a sentinel value.

    Exceptions & error handling:
    - Possible exceptions: FileNotFoundError, PermissionError, OSError, ValueError, RuntimeError.
        - FileNotFoundError: when input VCFs or required resources (e.g. VEP cache) are missing.
        - PermissionError/OSError: when unable to create or write into `work_dir` or temporary files.
        - ValueError: invalid parameter values (e.g. negative `dist`).
        - RuntimeError: underlying tools (SURVIVOR, VEP) fail during execution.
    - Caller should catch and handle these exceptions; for long-running steps (VEP) consider adding
        external retry, timeout, or job scheduling.

    Implementation notes:
    - Algorithmic idea: perform proximity-based SV merging using SURVIVOR, filter merged records by
        support vector (vec="01" to indicate presence only in the second input), then annotate the
        filtered VCF using VEP via the `VEP_SV` helper class.
    - Complexity: relies on external tools; this function does not maintain persistent internal state.
    - Dependencies: `VEP_SV` and its methods `merge_sv_survivor`, `extract_specific_sv`, and `annotate_sv_vep`.

    Performance:
    - Time complexity: dominated by VEP annotation and merge operations; roughly O(N) to O(N log N)
        where N is the number of input SV records. VEP typically dominates runtime and scales at least linearly.
    - Space complexity: disk I/O heavy; memory usage depends on VEP and SURVIVOR internals.
    - For large inputs (tens of thousands of SVs) run on nodes with sufficient CPU/memory and fast disk.

    Example:
    - annotated = run_exp_specific_annotation(
                dmso_vcf="/data/DMSO.sv.vcf.gz",
                plab_vcf="/data/PlaB.sv.vcf.gz",
                work_dir="/results/PlaB_SV",
                vep_cache="/opt/vep",
                species="mus_musculus",
                assembly="GRCm39",
                dist=500,
                min_support=1
        )

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

def run_vcf_analysis(
        vcf:str,
        outdir:str
):
    """
    Overall purpose:
    - Produce stratified statistics and visualizations for a SV VCF.
        This includes generation of count tables and two types of plots: a stacked bar
        chart for size bins and a smoothed curve of SV length counts.

    Use cases:
    - Post-annotation exploration of treatment-specific SVs to summarize distributions
        by SV type and size for reporting or downstream analyses.

    Data flow (input -> processing -> output):
    - Input: a single VCF file path and an output directory.
    - Processing: call `run_sv_stratification` to build stratified matrices and summaries,
        then call plotting helpers to render PNG files.
    - Output: files are written to `outdir/table` and `outdir/plot` (tables and images).

    Parameters:
    - vcf (str): required. Path to the PlaB-specific VCF. Must exist and be parseable by `parse_pbsv_vcf`.
    - outdir (str): required. Output directory where `table/` and `plot/` subdirectories will be created.
        Caller must ensure write permission.

    Returns:
    - None. Side effects include writing the following artifacts:
        - `outdir/table/*` : stratification tables produced by `run_sv_stratification`.
        - `outdir/plot/size_bin.png` : stacked bar chart of size bins.
        - `outdir/plot/svlen_count.png` : smoothed SV length count curve.

    Exceptions & error handling:
    - Possible exceptions: FileNotFoundError, PermissionError, OSError, RuntimeError.
        - FileNotFoundError: if the input VCF is missing or required intermediate files are absent.
        - PermissionError/OSError: if the function cannot create directories or write files.
        - RuntimeError: if downstream stratification or plotting functions fail.
    - Caller should validate the input VCF contains fields required for stratification (e.g. sv type, svlen).

    Implementation notes:
    - Algorithmic idea: build stratified count matrices using `run_sv_stratification`, then
        convert parsed VCF data into plotting-friendly structures (`generate_plot_input`) and
        call plotting routines.
    - Complexity: dominated by upstream parsing and stratification routines; function itself
        orchestrates those steps and calls plotting helpers.
    - Dependencies: `run_sv_stratification`, `parse_pbsv_vcf`, `generate_plot_input`,
        `plot_stacking_bar`, `plot_multi_smooth_curves`.

    Performance:
    - Time complexity: roughly O(N) where N is the number of SV records; memory usage
        grows with the size of intermediate data frames. For very large VCFs, consider
        sampling or running on machines with more memory.

    Example:
    - run_exp_specific_vcf_analysis(vcf="/results/PlaB_only_annotated.vcf", outdir="/results/PlaB")

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
    Overall purpose:
    - Extract insertion candidate sequences from foreground (PlaB-only) and background (DMSO)
        VCFs within a specified length range, annotate these sequences with RepeatMasker, and
        perform a comparative enrichment test at the repeat subfamily level.

    Use cases:
    - Investigate whether treatment-specific insertions are enriched for particular
        transposable element families/subfamilies as biological evidence supporting downstream experiments.

    Data flow (input -> processing -> output):
    - Input: two VCFs (foreground and background) and a length interval (min/max).
    - Processing: extract candidate insertion sequences -> run RepeatMasker for fg/bg ->
        use `RepeatMaskerOutCompare` to compute enrichment statistics -> write CSV and plot results.
    - Output: RepeatMasker outputs under `outdir/repeatmasker`, enrichment CSV and a PNG figure
        containing significant subfamilies (FDR < 0.05).

    Parameters:
    - vcf_01 (str): required. Path to the PlaB-only VCF (foreground). INS records will be extracted
        and filtered by length.
    - vcf_1x (str): required. Path to the DMSO VCF (background) used to construct the background annotation.
    - outdir (str): required. Root output directory where `fa/` and `repeatmasker/` subdirectories will be created.
    - min_svlen (int): required (default 2000). Minimum insertion length (bp). Must be >= 0 and <= max_svlen.
    - max_svlen (int): required (default 10000). Maximum insertion length (bp). Must be >= min_svlen.

    Returns:
    - None. Side effects include:
        - `outdir/repeatmasker/PlaBOnlyEnrichment.csv` : enrichment test results (subfamily, fg/bg counts, p-value, fdr, etc.).
        - `outdir/repeatmasker/PlaBOnlyEnrichment.png` : visualization for significant subfamilies (FDR < 0.05).

    Exceptions & error handling:
    - Possible exceptions: FileNotFoundError, PermissionError, OSError, ValueError, RuntimeError.
        - FileNotFoundError: when input VCFs or expected RepeatMasker outputs are missing.
        - ValueError: when `min_svlen` > `max_svlen` or invalid numeric arguments are provided.
        - RuntimeError: when RepeatMasker annotation or enrichment calculation fails.
    - Caller must ensure RepeatMasker and its databases are installed and available in the runtime environment.
        The pipeline creates intermediate files and therefore requires sufficient disk space.

    Implementation notes:
    - Algorithmic idea: use `extract_te_candidate_ins` to generate FASTA sequences for candidate insertions,
        run `run_te_annotation_pipeline` (RepeatMasker) on foreground and background FASTAs, then use
        `RepeatMaskerOutCompare.enrichment_test` to compute enrichment statistics at the subfamily level.
    - Complexity: the most expensive step is RepeatMasker annotation, which scales with the number
        and length of sequences and the size of the repeat database.
    - Dependencies: `extract_te_candidate_ins`, `run_te_annotation_pipeline`, `RepeatMaskerOutCompare`, `plot_enrichment`.

    Performance:
    - Time complexity: dominated by RepeatMasker runtime; expected to increase approximately linearly
        with the number of base pairs annotated but can be substantial for large datasets.
    - Space complexity: significant disk usage for FASTA and RepeatMasker outputs; memory requirements
        depend on RepeatMasker implementation and dataset size.
    - For very large candidate sets, consider batching or running on an HPC environment. Increasing
        `min_svlen` reduces candidate set size as a practical control.

    Example:
    - run_exp_enricher(
                vcf_01="/results/PlaB_only.vcf",
                vcf_1x="/data/DMSO.sv.vcf.gz",
                outdir="/results/PlaB_enrichment",
                min_svlen=2000,
                max_svlen=10000
        )

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
    paraser.add_argument("-o","--out_dir",type=str,required=True,help="out directory for intermediate and final results")
    paraser.add_argument("-d","--dist",type=int,default=500,help="SURVIVOR merge distance threshold")
    args = paraser.parse_args()

    annotated_vcf = run_exp_specific_annotation(
        dmso_vcf=args.dmso_vcf,
        plab_vcf=args.plab_vcf,
        work_dir=args.out_dir,
        dist=args.dist
    )

    run_vcf_analysis(
        vcf=annotated_vcf,
        outdir=f"{args.out_dir}/PlaB"
    )
    run_exp_enricher(
        vcf_01=f"{args.out_dir}/PlaB_only.vcf",
        vcf_1x=args.dmso_vcf,
        outdir=f"{args.out_dir}/PlaB_enrichment"
    )

if __name__ == "__main__":
    # main()
    config1 = {
        "vcf": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz",
        "outdir": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO_P6"
    }
    config2 = {
        "vcf": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz",
        "outdir": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB_P6"
    }
    config3 = {
        "vcf": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf",
        "outdir": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/DMSO_P20"
    }
    config4 = {
        "vcf": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf",
        "outdir": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB_P20"
    }
    run_vcf_analysis(**config1)
    run_vcf_analysis(**config2)
    run_vcf_analysis(**config3)
    run_vcf_analysis(**config4)