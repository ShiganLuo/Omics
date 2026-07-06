import argparse
from copy import deepcopy
import json
import os
from src.common.MetaUtil import MetadataUtils, DesignPair
from src.common.LogUtil import setup_logger
from src.common.CmdUtil import _run_cmd, _run_cmds_parallel
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
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount_name.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount_name.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
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
            outfiles.append(f"{outdir}/igv/dedup/{sample_id}.dedup.bam")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}/{sample_id}.single.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/dedup/{sample_id}.dedup.bam")
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
        outfiles.append(f"{outdir}/mutation/fragment_size/fragment/FragmentSize.txt")
        outfiles.append(f"{outdir}/mutation/fragment_size/fragment/FragmentSize.png")
    
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
    indir: str,
    outdir: str,
):
    """Prepare input JSON for PeakCalling (ChIP-seq/DIP-seq peak calling) workflow.
    
    Workflow steps:
    1. Trimming (trim_galore)
    2. Bowtie2 index
    3. Bowtie2 align
    4. MACS3 peak calling
    
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

    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")

        if sample_info.design == "ip":
            ip_samples.append(sample_id)
        elif sample_info.design == "input":
            input_samples.append(sample_id)
        else:
            logger.warning(f"Unknown design type for sample {sample_id}: {sample_info.design}, treating as IP")
            ip_samples.append(sample_id)

    # Build IP -> Input mapping
    # Match each IP sample with an Input sample (if available)
    # Strategy: use the first available Input sample as control for all IPs
    # More sophisticated matching can be implemented based on metadata
    if input_samples:
        default_input = input_samples[0]
        for ip_sample in ip_samples:
            sample_ip_input_map[ip_sample] = default_input
            # Add trimming and alignment outputs
            outfiles.append(f"{outdir}/cutadapt/{ip_sample}/{ip_sample}_1.fq.gz")
            outfiles.append(f"{outdir}/cutadapt/{ip_sample}/{ip_sample}_2.fq.gz")
            outfiles.append(f"{outdir}/bowtie2/mm/{ip_sample}/{ip_sample}.bam")
            # Add peak calling output
            outfiles.append(f"{outdir}/macs3/mm/{ip_sample}/{ip_sample}_peaks.narrowPeak")
        # Also add trimming/alignment for input samples
        for input_sample in input_samples:
            outfiles.append(f"{outdir}/cutadapt/{input_sample}/{input_sample}_1.fq.gz")
            outfiles.append(f"{outdir}/cutadapt/{input_sample}/{input_sample}_2.fq.gz")
            outfiles.append(f"{outdir}/bowtie2/mm/{input_sample}/{input_sample}.bam")
    else:
        logger.warning("No Input samples found. MACS3 will run without control.")
        for ip_sample in ip_samples:
            sample_ip_input_map[ip_sample] = None
            outfiles.append(f"{outdir}/cutadapt/{ip_sample}/{ip_sample}_1.fq.gz")
            outfiles.append(f"{outdir}/cutadapt/{ip_sample}/{ip_sample}_2.fq.gz")
            outfiles.append(f"{outdir}/bowtie2/mm/{ip_sample}/{ip_sample}.bam")
            outfiles.append(f"{outdir}/macs3/mm/{ip_sample}/{ip_sample}_peaks.narrowPeak")

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

def parse_args():
    parser = argparse.ArgumentParser(description="workflow")
    parser.add_argument('-m','--meta', type=str, required=True, help='meta input file or data dir which condatain fastq file')
    parser.add_argument('-w','--workflow_name', type=str, nargs='+',
        choices=["CoCulture", "MERIP", "RNAseq", "CLIP", "Mutation", "PacVar", "KARRseq", "PeakCalling", "QuantMS", "tRNAseq"],
        default=['CoCulture'], help='workflow name(s), multiple for parallel execution')
    parser.add_argument('-o','--output_dir', type=str, required=True, help='output dir')
    parser.add_argument('-t','--threads', type=int, default=10, help='threads')
    parser.add_argument('--dry-run', action='store_true', help='dry run')
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
    "CLIP":       lambda cfg, sid, dp, indir, outdir, meta: ("CLIP.smk",      runCLIP(cfg, sid, indir, outdir)),
    "Mutation":   lambda cfg, sid, dp, indir, outdir, meta: ("Mutation.smk",  runMutation(cfg, sid, dp, indir, outdir)),
    "PacVar":     lambda cfg, sid, dp, indir, outdir, meta: ("PacVar.smk",    runPacVar(cfg, sid, indir, outdir)),
    "KARRseq":    lambda cfg, sid, dp, indir, outdir, meta: ("KARRseq.smk",   runKARRseq(cfg, sid, indir, outdir)),
    "PeakCalling":lambda cfg, sid, dp, indir, outdir, meta: ("PeakCalling.smk",runPeakCalling(cfg, sid, indir, outdir)),
    "QuantMS":    lambda cfg, sid, dp, indir, outdir, meta: ("QuantMS.smk",   runQuantMS(cfg, sid, indir, outdir)),
    "tRNAseq":    lambda cfg, sid, dp, indir, outdir, meta: ("tRNAseq.smk",   runtRNAseq(cfg, sid, indir, outdir, meta)),
}


if __name__ == "__main__":
    args = parse_args()
    logger = setup_logger("root",level=logging.INFO, log_file=args.log)
    ROOT_DIR = os.path.dirname(__file__)

    workflow_names = args.workflow_name
    n_workflows = len(workflow_names)

    # Use first workflow's output dir for metadata (or a shared parent if multi)
    # Metadata is shared across workflows — only parsed once
    ref_outdir = os.path.join(args.output_dir, workflow_names[0])
    abs_ref_outdir = os.path.abspath(ref_outdir)
    if os.path.isfile(args.meta):
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, meta=args.meta)
    else:
        metadataUtil = MetadataUtils(outdir=abs_ref_outdir, fastq_dir=args.meta)
    samples_info_dict, designPair, raw_fastq_dir = metadataUtil.run()

    # Thread allocation: user-specified total threads split across workflows
    threads_per_workflow = max(1, args.threads // n_workflows)
    if n_workflows > 1:
        logger.info(f"Parallel mode: {n_workflows} workflows, "
                    f"{args.threads} total threads -> {threads_per_workflow} per workflow")

    # Prepare each workflow
    smk_cmds: list[tuple[list[str], str]] = []  # (cmd, cwd) pairs
    for wf_name in workflow_names:
        abs_outdir = os.path.abspath(os.path.join(args.output_dir, wf_name))
        os.makedirs(abs_outdir, exist_ok=True)

        model_json = os.path.join(ROOT_DIR, f"config/{wf_name}.json")
        workflow_config = _load_model_json(model_json)
        flat_args = {k: v for k, v in args.extra_args.items() if '.' not in k}
        workflow_config.update(flat_args)
        dot_args = parse_dot_args(args.extra_args)
        for key_tuple, v in dot_args.items():
            logger.info(f"Setting config parameter {'.'.join(key_tuple)} to {v}")
            dict_set_by_path(workflow_config, list(key_tuple), v)

        if wf_name not in WORKFLOW_DISPATCH:
            logger.error(f"Unknown workflow name: {wf_name}")
            exit(1)

        smk, input_json = WORKFLOW_DISPATCH[wf_name](
            deepcopy(workflow_config), samples_info_dict, designPair,
            raw_fastq_dir, abs_outdir,args.meta
        )

        cmd = build_snakemake_cmd(
            ROOT_DIR, smk, input_json, threads_per_workflow,
            args.conda_prefix, args.rerun_trigger, args.dry_run,
            args.conda_frontend, args.snakemake_args,
        )
        logger.info(f"[{wf_name}] {cmd}")
        # Each workflow runs from its own output dir to isolate .snakemake/
        smk_cmds.append((cmd, abs_outdir))

    if n_workflows == 1:
        _run_cmd(smk_cmds[0][0], cwd=smk_cmds[0][1])
    else:
        logger.info(f"Launching {n_workflows} snakemake processes in parallel...")
        _run_cmds_parallel(smk_cmds)
        logger.info("All workflows completed.")