import argparse
from copy import deepcopy
import json
import os
import re
from src.common.MetaUtil import MetadataUtils
from src.common.LogUtil import setup_logger
from src.common.CmdUtil import _run_cmd, _run_cmds_parallel
from src.common.type import DesignPair
from src.common.SchemaValidator import SchemaValidator
import logging
from typing import Dict, Any, Optional, List
logger = logging.getLogger(__name__)


def smart_cast(val):
    """尝试将字符串转换为 int/float/bool，否则原样返回；list 逐元素转换"""
    if isinstance(val, list):
        return [smart_cast(v) for v in val]
    if isinstance(val, str):
        if val.lower() in {"true", "false"}:
            return val.lower() == "true"
        try:
            if val.startswith("0") and len(val) > 1 and not val.startswith("0."):
                return val  # 避免八进制等
            return int(val)
        except Exception:
            pass
        try:
            return float(val)
        except Exception:
            pass
    return val

def dict_set_by_path(d, keys, value):
    """递归设置嵌套字典的值，keys为key列表，自动类型转换"""
    for k in keys[:-1]:
        if k not in d or not isinstance(d[k], dict):
            d[k] = {}
        d = d[k]
    d[keys[-1]] = smart_cast(value)

def parse_dot_args(extra_args):
    """从extra_args中提取点号语法参数，返回{(k1,k2,...):v}"""
    dot_args = {}
    for k, v in list(extra_args.items()):
        if '.' in k:
            dot_args[tuple(k.split('.'))] = v
    return dot_args


def _load_model_json(model_json_file: str) -> Dict[str, Any]:
    """Load model JSON template from disk."""
    with open(model_json_file, 'r', encoding='utf-8') as f:
        return json.load(f)
    
def runCoCulture(
    datajson: Dict[str,Any],
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,

):
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    outfiles = []
    paired_samples = []
    single_samples = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            outfiles.append(f"{outdir}/common/2_trimmed_fastq/{sample_id}/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/common/2_trimmed_fastq/{sample_id}/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCh38/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCm39/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount_name.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount_name.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/common/2_trimmed_fastq/{sample_id}/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCm39/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/common/3_raw_bam/GRCh38/{sample_id}/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount_name.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount_name.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/disambiguate/disambiguate_qc.tsv")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runMERIP(
    datajson: Dict[str, Any],
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
):
    """
    Function: Prepare input JSON for MERIP workflow based on the provided model JSON template and sample information.
    Parameters:
    - input_json: Path to the model JSON template file.
    - samples_info_dict: A dictionary containing sample information, where keys are sample IDs and values are objects with attributes 'layout' and 'design'.
    - indir: Input directory containing raw data (e.g., FASTQ files).
    - outdir: Output directory where results will be stored.
    Returns:
    - instance_json: Path to the generated input JSON file that will be used for the MERIP workflow.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    paired_samples = []
    single_samples = []
    ip_samples = []
    input_samples = []
    treated_ip_samples = []
    treated_input_samples = []
    outfiles = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_1.fq.gz")
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_2.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/{sample_id}/{sample_id}.dedup.bam")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}.single.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/{sample_id}/{sample_id}.dedup.bam")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
        
        if sample_info.design == "ip":
            ip_samples.append(sample_id)
        elif sample_info.design == "input":
            input_samples.append(sample_id)
        elif sample_info.design == "treated_ip":
            treated_ip_samples.append(sample_id)
        elif sample_info.design == "treated_input":
            treated_input_samples.append(sample_id)
        else:
            logger.error(f"Unknown design type for sample {sample_id}: {sample_info.design}")
    outfiles.append(f"{outdir}/exomePeak/sig_diff_peak_gene_names.xls")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["ip_samples"] = ip_samples
    datajson["input_samples"] = input_samples
    datajson["treated_ip_samples"] = treated_ip_samples
    datajson["treated_input_samples"] = treated_input_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runRNAseq(
    datajson: Dict[str, Any],
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
):
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    outfiles = []
    paired_samples = []
    single_samples = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            outfiles.append(f"{outdir}/stringtie/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt")
            outfiles.append(f"{outdir}/fusion/arriba/{sample_id}/{sample_id}_passed_fusions.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/stringtie/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt")
            outfiles.append(f"{outdir}/fusion/arriba/{sample_id}/{sample_id}_passed_fusions.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    # outfiles.append(f"{outdir}/TEtranscripts/TEcount/all_TEcount.tsv")
    outfiles.append(f"{outdir}/fusion/arriba/../arriba_report/arriba_fusion_report.html")
    outfiles.append(f"{outdir}/stringtie/stringtie_merged.gtf")
    outfiles.append(f"{outdir}/stringtie/result/TE_chimeric/TE_chimeric_group_stacked.png")
    outfiles.append(f"{outdir}/stringtie/result/TE_chimeric/TE_chimeric_te_type_top.png")
    outfiles.append(f"{outdir}/stringtie/result/TE_chimeric/TE_chimeric_te_type_by_group.png")
    outfiles.append(f"{outdir}/stringtie/result/TE_chimeric/TE_chimeric_sample_summary.tsv")
    outfiles.append(f"{outdir}/stringtie/result/TE_chimeric/TE_chimeric_group_summary.tsv")
    outfiles.append(f"{outdir}/stringtie/result/TE_chimeric/TE_chimeric_te_type_counts.tsv")
    sample_groups = {}
    for sample_id in samples_info_dict:
        group_key = sample_id.split('-', 1)[0]
        sample_groups.setdefault(group_key, []).append(sample_id)
    datajson["sample_groups"] = sample_groups
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runCLIP(
    datajson: Dict[str, Any],
    samples_info_dict:Dict[str, Any],
    indir:str,
    outdir: str,
):
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    outfiles = []
    paired_samples = []
    single_samples = []
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}_2.fq.gz")
            if datajson["aligner"] == "star":
                outfiles.append(f"{outdir}/star/{sample_id}/{sample_id}.bam")
            elif datajson["aligner"] == "hisat2":
                outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/fastqc/raw/{sample_id}/fastqc.raw.txt")
            outfiles.append(f"{outdir}/fastqc/trimmed/{sample_id}/fastqc.trimmed.txt")
            outfiles.append(f"{outdir}/PureCLIP/{sample_id}.pureclip.sites.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.plus.bw")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.minus.bw")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}.single.fq.gz")
            if datajson["aligner"] == "star":
                outfiles.append(f"{outdir}/star/{sample_id}/{sample_id}.bam")
            elif datajson["aligner"] == "hisat2":
                outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            outfiles.append(f"{outdir}/fastqc/raw/{sample_id}/fastqc.raw.txt")
            outfiles.append(f"{outdir}/fastqc/trimmed/{sample_id}/fastqc.trimmed.txt")
            outfiles.append(f"{outdir}/PureCLIP/{sample_id}.pureclip.sites.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.bed")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.plus.bw")
            outfiles.append(f"{outdir}/bedtools/{sample_id}/{sample_id}.minus.bw")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/track/igv_track_iclip.html")
    outfiles.append(f"{outdir}/track/ucsc_track_iclip.txt")
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    # parameters suggest by https://doi.org/10.1016/j.ymeth.2019.11.008
    datajson["Params"]["bamCoverage"]["offset"] = "-1"
    datajson["Params"]["bamCoverage"]["binSize"] = 1
    datajson["Params"]["bamCoverage"]["normalizeUsing"] = "CPM"
    datajson["Params"]["bamCoverage"]["extendReads"] = 1
    datajson["Params"]["STAR"]["alignEndsType"] = "Extend5pOfRead1"
    datajson["Params"]["STAR"]["outFilterMismatchNoverReadLmax"] = 0.04
    datajson["Params"]["STAR"]["outFilterMismatchNmax"] = 999
    datajson["Params"]["STAR"]["outFilterMultimapNmax"] = 999
    datajson["Params"]["igv"]["js"] = "/data/pub/zhousha/Reference/igv.min.js"
    datajson["Params"]["igv"]["publicPathMap"] = {
        "/data/pub/zhousha/": "/data/",
        "/data/pub/zhousha/Reference/": "/ref/"
    }
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runPacVar(
    datajson: Dict[str, Any],
    samples_info_dict: Dict[str, Any],
    indir: str,
    outdir: str,
):
    """Prepare input JSON for PacVar (PacBio variant calling) workflow."""
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    outfiles = []
    samples = []
    skip_snp = datajson.get("Params", {}).get("skip_snp", False)
    skip_sv = datajson.get("Params", {}).get("skip_sv", False)
    skip_phase = datajson.get("Params", {}).get("skip_phase", False)
    skip_repeat = datajson.get("Params", {}).get("skip_repeat", False)
    snv_caller = datajson.get("Params", {}).get("snv_caller", "deepvariant")

    for sample_id, sample_info in samples_info_dict.items():
        samples.append(sample_id)
        # SNP calling
        if not skip_snp:
            if snv_caller == "deepvariant":
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.vcf.gz")
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.vcf.gz.csi")
            elif snv_caller == "gatk4":
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.filtered.vcf.gz")
                outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.filtered.vcf.gz.csi")
        # SV calling
        if not skip_sv:
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.vcf.gz.csi")
        # phasing
        if not skip_phase and not skip_snp and not skip_sv:
            outfiles.append(f"{outdir}/variation/germline_snv_indel/{sample_id}/{sample_id}.phased.vcf.gz")
            outfiles.append(f"{outdir}/variation/germline_sv/{sample_id}/{sample_id}.sv.phased.vcf.gz")
        # repeat characterization
        if not skip_repeat and datajson["genome"]["repeat_bed"]:
            outfiles.append(f"{outdir}/repeat/trgt/genotype/{sample_id}/{sample_id}.trgt.vcf.gz")
            outfiles.append(f"{outdir}/repeat/trgt/plot/{sample_id}/{sample_id}.trgt.repeat.png")

    # telomere & centromere analysis
    skip_telomere = datajson.get("Params", {}).get("skip_telomere", False)
    if not skip_telomere:
        for sample_id in samples:
            # Telogator2 (per-chromosome-arm)
            outfiles.append(f"{outdir}/repeat/telomere/{sample_id}/telogator2/tlens_by_allele.tsv")
            outfiles.append(f"{outdir}/repeat/telomere/{sample_id}/telogator2/all_final_alleles.png")
            outfiles.append(f"{outdir}/repeat/telomere/{sample_id}/telogator2/violin_atl.png")
            # Approach A: Assembly contig end scanning
            outfiles.append(f"{outdir}/repeat/telomere/{sample_id}/assembly_scan/{sample_id}_assembly_telomere_stats.txt")
            # Approach B: Read-level k-mer density
            outfiles.append(f"{outdir}/repeat/telomere/{sample_id}/read_density/{sample_id}_read_telomere_stats.txt")
            # Approach C: tidk
            outfiles.append(f"{outdir}/repeat/telomere/{sample_id}/tidk/{sample_id}_tidk_telomeres.tsv")
            # Centromere
            outfiles.append(f"{outdir}/repeat/centromere/{sample_id}/{sample_id}.centromere_stats.txt")
    gatk_tmp_dir = os.path.join(outdir, "tmp")
    os.makedirs(gatk_tmp_dir, exist_ok=True)
    datajson["Params"]["gatk"]["tmp-dir"] = gatk_tmp_dir
    datajson["samples"] = samples
    datajson["outfiles"] = outfiles
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runMutation(
    datajson: Dict[str, Any],
    samples_info_dict:Dict[str, Any],
    designPairs: List[DesignPair],
    indir:str,
    outdir: str,
):
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    outfiles = []
    paired_samples = []
    single_samples = []
    mutect2_samples = []
    sample_somatic_vcf_dict = {}
    sample_group_dict = {}
    for designPair in designPairs:
        somatic_file = f"{outdir}/mutation/gatk/somatic/mutect2-vcf/{designPair.ctr_sample_id}_vs_{designPair.exp_sample_id}/{designPair.ctr_sample_id}_vs_{designPair.exp_sample_id}.vcf.gz"
        sample_somatic_vcf_dict[designPair.exp_sample_id] = somatic_file
        sample_group_dict[designPair.exp_sample_id] = designPair.exp_group
        outfiles.append(somatic_file)
        mutect2_samples.append(designPair.ctr_sample_id)
        mutect2_samples.append(designPair.exp_sample_id)
    
    for sample_id, sample_info in samples_info_dict.items():
        
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            if sample_id in mutect2_samples:
                logger.info(f"Sample {sample_id} is involved in mutect2 analysis, skipping germline workflow for this sample.")
                continue
            outfiles.append(f"{outdir}/mutation/gatk/germline/{sample_id}/{sample_id}.filtered.vcf.gz")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            if sample_id in mutect2_samples:
                logger.info(f"Sample {sample_id} is involved in mutect2 analysis, skipping germline workflow for this sample.")
                continue
            outfiles.append(f"{outdir}/mutation/gatk/germline/{sample_id}/{sample_id}.filtered.vcf.gz")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    outfiles.append(f"{outdir}/mutation/spectrum/somatic_spectrum_stacked_bar.png")
    
    all_samples = paired_samples + single_samples
    
    # Fragment size analysis outputs (unless skipped)
    skip_fragment_size = datajson.get("Params", {}).get("skip_fragment_size", False)
    if not skip_fragment_size:
        outfiles.append(f"{outdir}/results/fragment_size/fragment/FragmentSize.txt")
        outfiles.append(f"{outdir}/results/fragment_size/fragment/FragmentSize.png")
    
    # SV detection with Manta outputs (unless skipped)
    skip_sv = datajson.get("Params", {}).get("skip_sv", False)
    if not skip_sv:
        for sample_id in all_samples:
            outfiles.append(f"{outdir}/mutation/sv/manta/{sample_id}/results/variants/candidateSV.vcf.gz")
    
    # CNV detection with CNVkit outputs (unless skipped)
    skip_cnv = datajson.get("Params", {}).get("skip_cnv", False)
    if not skip_cnv:
        for sample_id in all_samples:
            outfiles.append(f"{outdir}/mutation/cnv/cnvkit/cnv/{sample_id}.cnr")
            outfiles.append(f"{outdir}/mutation/cnv/cnvkit/cnv/{sample_id}.cns")
    
    datajson["Params"]["somatic_spectrum"]["sample_somatic_vcf_dict"] = sample_somatic_vcf_dict
    datajson["Params"]["somatic_spectrum"]["sample_group_dict"] = sample_group_dict
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runKARRseq(
    datajson: Dict[str, Any],
    samples_info_dict: Dict[str, Any],
    indir: str,
    outdir: str,
):
    """Prepare input JSON for KARRseq (Kethoxal-Assisted RNA-RNA interaction sequencing) workflow."""
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    paired_samples = []
    single_samples = []
    outfiles = []

    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            # Final ligation pairs
            outfiles.append(f"{outdir}/chimeric/{sample_id}/{sample_id}.dedup.ligation.pairs.gz")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/chimeric/{sample_id}/{sample_id}.dedup.ligation.pairs.gz")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")

    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runPeakCalling(
    datajson: Dict[str, Any],
    samples_info_dict: Dict[str, Any],
    design_pairs:List[DesignPair],
    indir: str,
    outdir: str,
):
    """Prepare input JSON for PeakCalling (ChIP-seq/DIP-seq peak calling) workflow.
    
    Workflow steps:
    1. FastQC (raw)
    2. Trimming (trim_galore)
    3. FastQC (trimmed)
    4. Bowtie2 alignment
    5. AddReadsGroup + MarkDuplicates (GATK4)
    6. BigWig tracks (bamCoverage)
    7. MACS3 peak calling
    8. FRiP score
    9. HOMER peak annotation
    
    Supports both ChIP-seq and DIP-seq experiments.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    paired_samples = []
    single_samples = []
    ip_samples = []
    input_samples = []
    sample_ip_input_map = {}
    outfiles = []
    for design_pair in design_pairs:
        sample_ip_input_map[design_pair.exp_sample_id] = design_pair.ctr_sample_id
        ip_samples.append(design_pair.exp_sample_id)
        input_samples.append(design_pair.ctr_sample_id)

    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")

    # Build IP -> Input mapping
    # Match each IP sample with an Input sample (if available)
    # Strategy: use the first available Input sample as control for all IPs
    # More sophisticated matching can be implemented based on metadata
    if input_samples:
        default_input = input_samples[0]
        for ip_sample in ip_samples:
            sample_ip_input_map[ip_sample] = default_input
            # Step 5: AddReadsGroup + MarkDuplicates (GATK4)
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            # Step 6: BigWig tracks
            outfiles.append(f"{outdir}/tracks/{ip_sample}/{ip_sample}.bigwig")
            # Step 7: MACS3 peak calling
            outfiles.append(f"{outdir}/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            # Step 8: FRiP score
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            # Step 9: HOMER annotation
            outfiles.append(f"{outdir}/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
        # Also add markdup for input samples
        for input_sample in input_samples:
            outfiles.append(f"{outdir}/common/4_markdup_bam/{input_sample}/{input_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/tracks/{input_sample}/{input_sample}.bigwig")
    else:
        logger.warning("No Input samples found. MACS3 will run without control.")
        for ip_sample in ip_samples:
            sample_ip_input_map[ip_sample] = None
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/tracks/{ip_sample}/{ip_sample}.bigwig")
            outfiles.append(f"{outdir}/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            outfiles.append(f"{outdir}/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
    outfiles.append(f"{outdir}/tracks/ucsc_track.txt")
    outfiles.append(f"{outdir}/tracks/igv_track.html")
    outfiles.append(f"{outdir}/PeakCalling_report.pptx")
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["samples"] = paired_samples + single_samples
    datajson["ip_samples"] = ip_samples
    datajson["input_samples"] = input_samples
    datajson["sample_ip_input_map"] = sample_ip_input_map
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runQuantMS(
    datajson: Dict[str, Any],
    samples_info_dict: Dict[str, Any],
    indir: str,
    outdir: str,
):
    """Prepare input JSON for QuantMS (quantitative proteomics) workflow.
    
    Workflow steps:
    1. Decoy database generation
    2. Database search engines (Comet, MSGF+, Sage)
    3. PSM rescoring (Percolator)
    4. PSM FDR control
    5. Protein inference (EpiFany)
    6. Protein quantification (ProteomicsLFQ or ProteinQuantifier)
    7. Statistical analysis (MSstats)
    
    Supports TMT, LFQ, and DIA quantification methods.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    
    samples = []
    mzml_files = []
    outfiles = []
    
    for sample_id, sample_info in samples_info_dict.items():
        samples.append(sample_id)
        # For proteomics, we expect mzML files in the input directory
        mzml_file = os.path.join(indir, f"{sample_id}.mzML")
        if not os.path.exists(mzml_file):
            # Try with .mzML.gz extension
            mzml_file_gz = os.path.join(indir, f"{sample_id}.mzML.gz")
            if os.path.exists(mzml_file_gz):
                mzml_file = mzml_file_gz
            else:
                logger.warning(f"mzML file not found for sample {sample_id}: {mzml_file}")
                continue
        mzml_files.append(mzml_file)
    
    # Build outfiles based on quantification method
    quantification_method = datajson.get("quantification_method", "lfq")
    
    # Decoy database
    outfiles.append(f"{outdir}/decoy_database/{os.path.basename(datajson['genome']['fasta'])}_decoy.fasta")
    
    # Database search results
    for sample_id in samples:
        outfiles.append(f"{outdir}/search_engine/{sample_id}/{sample_id}.idXML")
    
    # PSM rescoring results
    for sample_id in samples:
        outfiles.append(f"{outdir}/psm_rescoring/{sample_id}/{sample_id}_scored.idXML")
    
    # PSM FDR control results
    for sample_id in samples:
        outfiles.append(f"{outdir}/psm_fdr/{sample_id}/{sample_id}_filtered.idXML")
    
    # Protein inference results
    for sample_id in samples:
        outfiles.append(f"{outdir}/protein_inference/{sample_id}/{sample_id}_protein.idXML")
    
    # Quantification results
    if quantification_method == "tmt":
        outfiles.append(f"{outdir}/quantification/tmt_quantification.mzTab")
    elif quantification_method == "lfq":
        outfiles.append(f"{outdir}/quantification/lfq_quantification.mzTab")
    elif quantification_method == "dia":
        outfiles.append(f"{outdir}/quantification/dia_quantification.mzTab")
    
    # MSstats results
    if not datajson.get("Params", {}).get("skip_post_msstats", False):
        outfiles.append(f"{outdir}/msstats/msstats_results.csv")
    
    datajson["samples"] = samples
    datajson["mzml_files"] = mzml_files
    datajson["outfiles"] = outfiles
    
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runtRNAseq(
    datajson: Dict[str, Any],
    samples_info_dict: Dict[str, Any],
    indir: str,
    outdir: str,
    meta: str,
):
    """Prepare input JSON for tRNAseq (mim-tRNAseq) workflow.

    mim-tRNAseq is an all-in-one pipeline for tRNA sequencing analysis:
    tRNA clustering, GSNAP alignment, coverage analysis, modification
    quantification, CCA analysis, and DESeq2 differential expression.

    The pipeline processes all samples together via a sample data sheet.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    if not datajson.get("Params", {}).get("mimseq", {}).get("data_dir"):
        datajson.setdefault("Params", {}).setdefault("mimseq", {})["data_dir"] = os.path.join(os.path.dirname(__file__), "modules", "mimseq", "mimseq", "data")
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    datajson["meta"] = meta
    samples = []
    for sample_id in samples_info_dict:
        samples.append(sample_id)

    outfiles = [f"{outdir}/mimseq/mimseq.done"]

    datajson["samples"] = samples
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runncRNAseq(
    datajson: Dict[str, Any],
    samples_info_dict: Dict[str, Any],
    indir: str,
    outdir: str,
):
    """Prepare input JSON for ncRNAseq (small/non-coding RNA-seq) workflow.

    Pipeline: cutadapt trim -> hisat2 (no-spliced-alignment) or star -> featureCounts.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    outfiles = []
    paired_samples = []
    single_samples = []

    for sample_id, sample_info in samples_info_dict.items():
        outfiles.append(f"{outdir}/common/3_raw_bam/{sample_id}/{sample_id}.bam")
        outfiles.append(f"{outdir}/results/tailer/{sample_id}/{sample_id}_tail.csv")
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")

    all_samples = paired_samples + single_samples
    datajson["samples"] = all_samples
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def parse_args():
    parser = argparse.ArgumentParser(description="workflow")
    parser.add_argument('-m','--meta', type=str, default=None, help='meta input file or data dir which condatain fastq file')
    parser.add_argument('-w','--workflow_name', type=str, nargs='+',
        choices=["CoCulture", "MERIP", "RNAseq", "ncRNAseq", "CLIP", "Mutation", "PacVar", "KARRseq", "PeakCalling", "QuantMS", "tRNAseq"],
        default=['CoCulture'], help='workflow name(s), multiple for parallel execution')
    parser.add_argument('-o','--output_dir', type=str, default=None, help='output dir')
    parser.add_argument('-t','--threads', type=int, default=10, help='threads')
    parser.add_argument('--dry-run', action='store_true', help='dry run')
    parser.add_argument('--test', type=str, nargs='?', const='all', metavar='WORKFLOW',
        help='run dry-run test for a workflow (or "all"). Auto-sets meta/output/dry-run from test/ directory')
    parser.add_argument('--log', type=str, default='workflow.log', help='log file')
    parser.add_argument('--conda-prefix', type=str, default='/data/pub/zhousha/env/mutation_0.1', help='conda prefix for snakemake')
    parser.add_argument(
        '--rerun-trigger', '--rerun-triggers',
        dest='rerun_trigger',
        nargs='+',
        default=["input"],
        choices=["code", "input", "mtime", "params", "software-env"],
        help='snakemake rerun-triggers, e.g. code input mtime params software-env'
    )
    parser.add_argument('--conda-frontend', type=str, choices=["conda", "mamba"], default="mamba", help='conda frontend for snakemake')
    parser.add_argument(
        '--snakemake-args',
        nargs=argparse.REMAINDER,
        default=[],
        help='additional arguments forwarded to snakemake; place them after this flag'
    )
    
    # 支持 --key=value、--key value、--key v1 v2 v3 三种形式的额外参数
    # 多值参数在碰到下一个 --key 或到达末尾时停止收集
    args, unknown = parser.parse_known_args()
    extra_args = {}
    i = 0
    while i < len(unknown):
        arg = unknown[i]
        if arg.startswith('--'):
            key = arg[2:]
            if '=' in key:
                k, v = key.split('=', 1)
                extra_args[k] = v
            else:
                values = []
                while i + 1 < len(unknown) and not unknown[i + 1].startswith('--'):
                    values.append(unknown[i + 1])
                    i += 1
                if not values:
                    extra_args[key] = True
                elif len(values) == 1:
                    extra_args[key] = values[0]
                else:
                    extra_args[key] = values
        i += 1
    args.extra_args = extra_args
    return args


def build_snakemake_cmd(root_dir, smk, input_json, threads, conda_prefix, rerun_trigger, dry_run, conda_frontend, snakemake_args):
    cmd = [
        "snakemake",
        "-s",
        f"{root_dir}/subworkflow/{smk}",
        "--configfile",
        input_json,
        "--cores",
        str(threads),
        "--conda-prefix",
        conda_prefix,
        "--rerun-triggers",
        *rerun_trigger,
        "--use-conda",
        "--conda-frontend",
        conda_frontend,
    ]
    if dry_run:
        cmd.append("--dry-run")
    if snakemake_args:
        cmd.extend(snakemake_args)
    return cmd


WORKFLOW_DISPATCH = {
    "CoCulture":  lambda cfg, sid, dp, indir, outdir, meta: ("CoCulture.smk", runCoCulture(cfg, sid, indir, outdir)),
    "MERIP":      lambda cfg, sid, dp, indir, outdir, meta: ("MERIP.smk",     runMERIP(cfg, sid, indir, outdir)),
    "RNAseq":     lambda cfg, sid, dp, indir, outdir, meta: ("RNAseq.smk",    runRNAseq(cfg, sid, indir, outdir)),
    "ncRNAseq":   lambda cfg, sid, dp, indir, outdir, meta: ("ncRNAseq.smk",  runncRNAseq(cfg, sid, indir, outdir)),
    "CLIP":       lambda cfg, sid, dp, indir, outdir, meta: ("CLIP.smk",      runCLIP(cfg, sid, indir, outdir)),
    "Mutation":   lambda cfg, sid, dp, indir, outdir, meta: ("Mutation.smk",  runMutation(cfg, sid, dp, indir, outdir)),
    "PacVar":     lambda cfg, sid, dp, indir, outdir, meta: ("PacVar.smk",    runPacVar(cfg, sid, indir, outdir)),
    "KARRseq":    lambda cfg, sid, dp, indir, outdir, meta: ("KARRseq.smk",   runKARRseq(cfg, sid, indir, outdir)),
    "PeakCalling":lambda cfg, sid, dp, indir, outdir, meta: ("PeakCalling.smk",runPeakCalling(cfg, sid,dp, indir, outdir)),
    "QuantMS":    lambda cfg, sid, dp, indir, outdir, meta: ("QuantMS.smk",   runQuantMS(cfg, sid, indir, outdir)),
    "tRNAseq":    lambda cfg, sid, dp, indir, outdir, meta: ("tRNAseq.smk",   runtRNAseq(cfg, sid, indir, outdir, meta)),
}


# ============================================================
# Test path generation
# ============================================================


def setup_test_args(args, root_dir: str):
    """Configure args for --test mode.

    Resolves workflow names, output directory, meta files, and test paths.
    Modifies args in-place and returns it.
    """
    import shutil

    TEST_DIR = os.path.join(root_dir, "assests", "test")

    # All registered workflows
    ALL_WORKFLOWS = list(WORKFLOW_DISPATCH.keys())
    if args.test == "all":
        args.workflow_name = ALL_WORKFLOWS
    elif args.test in ALL_WORKFLOWS:
        args.workflow_name = [args.test]
    else:
        logger.info(f"Unknown workflow: {args.test}")
        logger.info(f" Available: {ALL_WORKFLOWS} or 'all'")
        exit(1)

    # Output to {cwd or --output-dir}/test
    base_out = args.output_dir if args.output_dir else os.getcwd()
    args.output_dir = os.path.join(base_out, "test")
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    logger.info(f"Output: {args.output_dir}")

    args.dry_run = True
    args.log = os.path.join(args.output_dir, "test.log")
    os.makedirs(args.output_dir, exist_ok=True)

    # Build per-workflow meta map dynamically
    args._test_meta_map = {}
    for wf in args.workflow_name:
        meta_path = os.path.join(TEST_DIR, f"meta_{wf}.tsv")
        if os.path.isfile(meta_path):
            args._test_meta_map[wf] = meta_path
        else:
            logger.info(f"Warning: no meta for {wf} at {meta_path}")

    args.meta = None  # will be resolved per-workflow via _get_meta()

    # Store schema validator and test data for per-workflow path injection
    test_data = os.path.join(TEST_DIR, "data")
    GENOME = "GRCm39"
    schema_validator = SchemaValidator()
    schema_validator._schema_dir = os.path.join(root_dir, "config")
    args._test_base_paths = schema_validator.generate_test_paths(test_data, GENOME)
    args._test_genome = GENOME

    # Local conda-prefix (avoid permission issues)
    args.conda_prefix = os.path.join(args.output_dir, ".conda")
    os.makedirs(args.conda_prefix, exist_ok=True)

    return args


def setup_normal_args(args):
    """Validate args for normal (non-test) mode."""
    if not args.meta:
        logger.info("Error: -m/--meta is required (unless --test is used)")
        exit(1)
    if not args.output_dir:
        logger.info("Error: -o/--output_dir is required (unless --test is used)")
        exit(1)
    args._test_meta_map = None
    return args


def print_test_summary(test_results: Dict[str, tuple]):
    """Print summary of test workflow results."""
    logger.info(f"[Results ({len(test_results)} workflows)")
    passed = [k for k, (ok, _) in test_results.items() if ok]
    failed = [k for k, (ok, _) in test_results.items() if not ok]
    for wf in sorted(test_results.keys()):
        ok, err = test_results[wf]
        status = "PASS" if ok else "FAIL"
        logger.info(f"  [{status}] {wf}")
        if not ok and err:
            logger.info(f"         {err.splitlines()[0] if err else 'unknown error'}")
    logger.info(f"\n  Passed: {len(passed)}, Failed: {len(failed)}, Total: {len(test_results)}")
    if failed:
        logger.error(f"  Failed workflows: {', '.join(failed)}")
        exit(1)


def execute_workflows(args, root_dir: str, logger):
    """Execute all configured workflows.

    In test mode (args._test_meta_map is set), runs each workflow and collects results.
    In normal mode, runs workflows directly.
    """
    workflow_names = args.workflow_name
    n_workflows = len(workflow_names)

    def _get_meta(wf_name):
        if args._test_meta_map:
            return args._test_meta_map.get(wf_name, args.meta)
        return args.meta

    # Use first workflow's output dir for metadata (or a shared parent if multi)
    ref_outdir = os.path.join(args.output_dir, workflow_names[0])
    abs_ref_outdir = os.path.abspath(ref_outdir)
    first_meta = _get_meta(workflow_names[0])
    if first_meta and os.path.isfile(first_meta):
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, meta=first_meta)
    else:
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, fastq_dir=first_meta)
    samples_info_dict, designPair, raw_fastq_dir = metadataUtil.run()

    # Thread allocation: user-specified total threads split across workflows
    threads_per_workflow = max(1, args.threads // n_workflows)
    if n_workflows > 1:
        logger.info(f"Parallel mode: {n_workflows} workflows, "
                    f"{args.threads} total threads -> {threads_per_workflow} per workflow")

    # Prepare each workflow
    test_results = {}  # wf_name -> (passed: bool, error: str)
    smk_cmds: list[tuple[list[str], str]] = []  # (cmd, cwd) pairs
    for wf_name in workflow_names:
        abs_outdir = os.path.abspath(os.path.join(args.output_dir, wf_name))
        os.makedirs(abs_outdir, exist_ok=True)

        try:
            # In test mode, reload metadata per workflow (different meta files)
            if args._test_meta_map:
                wf_meta = _get_meta(wf_name)
                if wf_meta and os.path.isfile(wf_meta):
                    metadataUtil = MetadataUtils(outdir=abs_outdir, meta=wf_meta)
                else:
                    metadataUtil = MetadataUtils(outdir=abs_outdir, fastq_dir=wf_meta)
                samples_info_dict, designPair, raw_fastq_dir = metadataUtil.run()

            model_json = os.path.join(root_dir, f"config/{wf_name}.json")
            workflow_config = _load_model_json(model_json)

            # In test mode, inject test paths for ALL path-like fields in config
            if hasattr(args, '_test_base_paths'):
                test_data = os.path.join(os.path.join(root_dir, "assests", "test"), "data")
                base_paths = args._test_base_paths

                def _is_path(val):
                    if val is None: return True
                    if not isinstance(val, str): return False
                    if "/" in val: return True
                    return False

                def _make_test_path(key, test_data, genome):
                    from pathlib import Path
                    ref = Path(test_data) / "ref"
                    ref.mkdir(parents=True, exist_ok=True)
                    for name in ("smallrna", "rRNA", "access", "repeat", "decoy"):
                        if name in key: return str(ref / name)
                    if key.endswith(("_dir", "_index")):
                        d = Path(test_data) / "index" / key.replace("_dir", "").replace("_index", "")
                        d.mkdir(parents=True, exist_ok=True)
                        return str(d)
                    if key.endswith("_prefix"):
                        d = Path(test_data) / "index" / key.replace("_index_prefix", "").replace("_prefix", "")
                        d.mkdir(parents=True, exist_ok=True)
                        return str(d / genome)
                    return str(ref / genome)

                def _inject(cfg, prefix, wf_extra):
                    for field, val in cfg.items():
                        dotted = f"{prefix}.{field}" if prefix else field
                        if isinstance(val, dict):
                            _inject(val, dotted, wf_extra)
                        elif _is_path(val):
                            wf_extra[dotted] = base_paths.get(dotted, _make_test_path(field, test_data, args._test_genome))

                wf_extra = dict(args.extra_args) if hasattr(args, 'extra_args') else {}
                _inject(workflow_config, "", wf_extra)

                flat_args = {k: v for k, v in wf_extra.items() if '.' not in k}
                workflow_config.update(flat_args)
                dot_args = parse_dot_args(wf_extra)
                for key_tuple, v in dot_args.items():
                    dict_set_by_path(workflow_config, list(key_tuple), v)
            else:
                flat_args = {k: v for k, v in args.extra_args.items() if '.' not in k}
                workflow_config.update(flat_args)
                dot_args = parse_dot_args(args.extra_args)
                for key_tuple, v in dot_args.items():
                    dict_set_by_path(workflow_config, list(key_tuple), v)

            if wf_name not in WORKFLOW_DISPATCH:
                raise ValueError(f"Unknown workflow name: {wf_name}")

            smk, input_json = WORKFLOW_DISPATCH[wf_name](
                deepcopy(workflow_config), samples_info_dict, designPair,
                raw_fastq_dir, abs_outdir, _get_meta(wf_name)
            )

            cmd = build_snakemake_cmd(
                root_dir, smk, input_json, threads_per_workflow,
                args.conda_prefix, args.rerun_trigger, args.dry_run,
                args.conda_frontend, args.snakemake_args,
            )
            logger.info(f"[{wf_name}] {cmd}")
            smk_cmds.append((cmd, abs_outdir))

        except Exception as e:
            if args._test_meta_map:
                test_results[wf_name] = (False, str(e))
                logger.error(f"[{wf_name}] Config/build failed: {e}")
            else:
                raise

    # Execute snakemake commands
    if args._test_meta_map:
        # Test mode: run each workflow, catch errors
        for cmd, cwd in smk_cmds:
            wf = os.path.basename(cwd)
            try:
                _run_cmd(cmd, cwd=cwd)
                test_results[wf] = (True, "")
            except Exception as e:
                test_results[wf] = (False, str(e)[:200])

        print_test_summary(test_results)
    else:
        if n_workflows == 1:
            _run_cmd(smk_cmds[0][0], cwd=smk_cmds[0][1])
        else:
            logger.info(f"Launching {n_workflows} snakemake processes in parallel...")
            _run_cmds_parallel(smk_cmds)
            logger.info("All workflows completed.")


if __name__ == "__main__":
    args = parse_args()
    ROOT_DIR = os.path.dirname(__file__)

    # Configure args based on mode
    if args.test is not None:
        setup_test_args(args, ROOT_DIR)
    else:
        setup_normal_args(args)

    logger = setup_logger("root", level=logging.INFO, log_file=args.log)
    execute_workflows(args, ROOT_DIR, logger)