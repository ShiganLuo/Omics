import os
from typing import Dict, Any, List
import logging
import json
from src.common.util.type import DesignPair, CompareGroupPair, SampleInfo
from src.common.util.LogUtil import setup_logger
logger = setup_logger(__name__, level=logging.DEBUG)
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
    samples_info_dict:Dict[str, SampleInfo],
    group_pairs: List[CompareGroupPair],
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
    sample_groups = {}
    for group_pair in group_pairs:
        datajson["Params"]["DESeq2"]["group_pairs"].setdefault(f"{group_pair.ctr_group_token}_vs_{group_pair.exp_group_token}", {
            "control_group_name": group_pair.ctr_group_name,
            "experimental_group_name": group_pair.exp_group_name,
            "control_samples": group_pair.ctr_sample_ids,
            "experimental_samples": group_pair.exp_sample_ids
        })
        outfiles.append(f"{outdir}/diff_expression/{group_pair.ctr_group_token}_vs_{group_pair.exp_group_token}/DESeq2.done")
        sample_groups.setdefault(group_pair.ctr_group_name, []).extend(group_pair.ctr_sample_ids)
        sample_groups.setdefault(group_pair.exp_group_name, []).extend(group_pair.exp_sample_ids)

        if datajson.get("Params", {}).get("function", {}).get("enabled", False):
            pair_dir = f"{outdir}/function/{group_pair.ctr_group_token}_vs_{group_pair.exp_group_token}"
            outfiles.append(f"{pair_dir}/go_back_to_back.png")
            outfiles.append(f"{pair_dir}/kegg_back_to_back.png")
            outfiles.append(f"{pair_dir}/GSEA/TEcount_Gene_GSEA.jpeg")
    datajson["Params"]["StringTie"]["sample_groups"] = sample_groups
    Organisms = set()
    for sample_id, sample_info in samples_info_dict.items():
        Organisms.add(sample_info.organism)
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            outfiles.append(f"{outdir}/transcripts/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt")
            outfiles.append(f"{outdir}/fusion/{sample_id}/{sample_id}_passed_fusions.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            outfiles.append(f"{outdir}/transcripts/raw/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt")
            outfiles.append(f"{outdir}/fusion/{sample_id}/{sample_id}_passed_fusions.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    if len(Organisms) != 1:
        raise ValueError(f"meta don't support multiple organsim temporarily, please check your meta file, found: {Organisms}")
    if Organisms.pop() in ["Homo sapiens", "human"]:
        datajson["genome"]["default"] = "GRCh38"
        datajson["Params"]["report"]["genome"] = "GRCh38"
        datajson["Params"]["function"]["species"] = "human"
    elif Organisms.pop() in ["Mus musculus", "mouse"]:
        datajson["genome"]["default"] = "GRCm39"
        datajson["Params"]["report"]["genome"] = "GRCm39"
        datajson["Params"]["function"]["species"] = "mouse"
    else:
        raise ValueError(f"pipeline don't support {Organisms.pop()}, only support human or mouse(Homo sapiens or Mus musculus)")
    # outfiles.append(f"{outdir}/TEtranscripts/TEcount/all_TEcount.tsv")
    outfiles.append(f"{outdir}/fusion/arriba_report/arriba_fusion_report.html")
    outfiles.append(f"{outdir}/transcripts/stringtie_merged.gtf")
    outfiles.append(f"{outdir}/transcripts/TE_chimeric/TE_chimeric_group_stacked.png")
    outfiles.append(f"{outdir}/transcripts/TE_chimeric/TE_chimeric_te_type_top.png")
    outfiles.append(f"{outdir}/transcripts/TE_chimeric/TE_chimeric_te_type_by_group.png")
    outfiles.append(f"{outdir}/transcripts/TE_chimeric/TE_chimeric_sample_summary.tsv")
    outfiles.append(f"{outdir}/transcripts/TE_chimeric/TE_chimeric_group_summary.tsv")
    outfiles.append(f"{outdir}/transcripts/TE_chimeric/TE_chimeric_te_type_counts.tsv")
    outfiles.append(f"{outdir}/RNAseq_report.pptx")
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
    datajson["Params"]["star"]["alignEndsType"] = "Extend5pOfRead1"
    datajson["Params"]["star"]["outFilterMismatchNoverReadLmax"] = 0.04
    datajson["Params"]["star"]["outFilterMismatchNmax"] = 999
    datajson["Params"]["star"]["outFilterMultimapNmax"] = 999
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

    Pipeline: jla-demultiplexer -> trim_galore -> subsample -> STAR (star / star_3pass / star_3pass_gene) -> featureCounts + Tailer.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    aligner = datajson.get("Procedure", {}).get("aligner") or "star_3pass"
    if aligner == "star_3pass":
        bam_subdir = "common/3_raw_bam"
    elif aligner == "star_3pass_gene":
        bam_subdir = "common/4_per_gene_bam"
    else:
        bam_subdir = "common/3_raw_bam"

    outfiles = []
    paired_samples = []
    single_samples = []

    for sample_id, sample_info in samples_info_dict.items():
        outfiles.append(f"{outdir}/{bam_subdir}/{sample_id}/{sample_id}.bam")
        outfiles.append(f"{outdir}/{bam_subdir}/{sample_id}/{sample_id}_tail.csv")
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")

    all_samples = paired_samples + single_samples
    outfiles.append(f"{outdir}/ncRNAseq_report.pptx")
    datajson["samples"] = all_samples
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json
