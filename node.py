import os
import time
from typing import Dict, Any, List
import logging
import json
from src.common.util.type import DesignPair, CompareGroupPair, SampleInfo, CellrangerInput
from src.common.util.LogUtil import setup_logger
logger = setup_logger(__name__, level=logging.DEBUG)

def runCoCulture(
        datajson: Dict[str,Any],
        samples_info_dict:Dict[str, SampleInfo],
        indir:str,
        outdir: str,
        raw_files: List[str],
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
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runMERIP(
        datajson: Dict[str, Any],
        samples_info_dict:Dict[str, SampleInfo],
        indir:str,
        outdir: str,
        raw_files: List[str],
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
    datajson["raw_files"] = raw_files
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


def runCLIP(
        datajson: Dict[str, Any],
        samples_info_dict:Dict[str, SampleInfo],
        indir:str,
        outdir: str,
        raw_files: List[str],
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
    datajson["raw_files"] = raw_files
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
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
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
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runMutation(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        designPairs: List[DesignPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
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
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runKARRseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
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
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runPeakCalling(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        design_pairs:List[DesignPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
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
    9. deeptools enrichment heatmap
    10. HOMER peak annotation
    
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

    # Auto-detect organism and set genome default (RNAseq pattern)
    organisms = set()
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.organism:
            organisms.add(sample_info.organism)
    if len(organisms) == 1:
        organism = next(iter(organisms))
        if organism in ["Homo sapiens", "human"]:
            datajson["genome"]["default"] = "GRCh38"
            datajson["Params"]["report"]["genome"] = "GRCh38"
        elif organism in ["Mus musculus", "mouse"]:
            datajson["genome"]["default"] = "GRCm39"
            datajson["Params"]["report"]["genome"] = "GRCm39"
        else:
            logger.warning(f"Unknown organism '{organism}', using config default genome")
    elif len(organisms) > 1:
        logger.warning(f"Multiple organisms detected: {organisms}, using config default genome")
    datajson["Params"]["report"]["date"] = time.strftime("%Y-%m-%d", time.localtime())

    # Heatmap samples: from computeMatrix.samples, fallback to all ip_samples
    cm_params = datajson.get("Params", {}).get("computeMatrix", {})
    cm_heatmap_samples = cm_params.get("samples", None) or ip_samples
    # Heatmap output suffix based on regions mode
    cm_regions = cm_params.get("regions", "tss")
    if cm_regions == "tss":
        hm_suffix = "tss"
    elif cm_regions == "peaks":
        hm_suffix = "peaks"
    elif isinstance(cm_regions, dict) and "genes" in cm_regions:
        hm_suffix = "per_gene" if cm_regions.get("per_gene") else "genes"
    else:
        hm_suffix = "heatmap"
    # Gene names for per_gene mode
    cm_gene_names = cm_regions.get("genes", []) if isinstance(cm_regions, dict) else []

    if input_samples:
        for ip_sample in ip_samples:
            outfiles.append(f"{outdir}/QC/1_raw_fastqc/{ip_sample}/{ip_sample}.fastqc.raw.txt")
            outfiles.append(f"{outdir}/QC/2_trimmed_fastqc/{ip_sample}/{ip_sample}.fastqc.trimmed.txt")
            # Step 5: AddReadsGroup + MarkDuplicates (GATK4)
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            # Step 6: BigWig tracks
            outfiles.append(f"{outdir}/results/tracks/{ip_sample}/{ip_sample}.bigwig")
            # Step 7: MACS3 peak calling
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_broad_peaks.broadPeak")
            # Step 8: FRiP score
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            # Step 9: deeptools enrichment heatmap
            if ip_sample in cm_heatmap_samples:
                if hm_suffix == "per_gene" and cm_gene_names:
                    for gene_name in cm_gene_names:
                        outfiles.append(f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{gene_name}_heatmap.png")
                else:
                    outfiles.append(f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{hm_suffix}_heatmap.png")
            # Step 10: HOMER annotation
            outfiles.append(f"{outdir}/results/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
            # Step 11: Peak-TE overlap
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_te_subfamily_overlap.tsv")
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_peak_centric_te.tsv")
        outfiles.append(f"{outdir}/results/peaks/cutoff_analysis.png")
        outfiles.append(f"{outdir}/results/te_overlap/te_subfamily_overlap_combined.png")
        for input_sample in input_samples:
            outfiles.append(f"{outdir}/common/4_markdup_bam/{input_sample}/{input_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/results/tracks/{input_sample}/{input_sample}.bigwig")
    else:
        logger.warning("No Input samples found. MACS3 will run without control.")
        for ip_sample in ip_samples:
            sample_ip_input_map[ip_sample] = None
            outfiles.append(f"{outdir}/common/4_markdup_bam/{ip_sample}/{ip_sample}.sorted_markdup.bam")
            outfiles.append(f"{outdir}/results/tracks/{ip_sample}/{ip_sample}.bigwig")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
            outfiles.append(f"{outdir}/results/peaks/{ip_sample}/{ip_sample}_broad_peaks.broadPeak")
            outfiles.append(f"{outdir}/QC/3_frip_score/{ip_sample}/{ip_sample}.FRiP.txt")
            # Step 9: deeptools enrichment heatmap
            if ip_sample in cm_heatmap_samples:
                if hm_suffix == "per_gene" and cm_gene_names:
                    for gene_name in cm_gene_names:
                        outfiles.append(f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{gene_name}_heatmap.png")
                else:
                    outfiles.append(f"{outdir}/results/heatmap/{ip_sample}/{ip_sample}_{hm_suffix}_heatmap.png")
            # Step 10: HOMER annotation
            outfiles.append(f"{outdir}/results/annotation/{ip_sample}/{ip_sample}_peaks.annotatePeaks.txt")
            # Step 11: Peak-TE overlap
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_te_subfamily_overlap.tsv")
            outfiles.append(f"{outdir}/results/te_overlap/{ip_sample}/{ip_sample}_peak_centric_te.tsv")
    outfiles.append(f"{outdir}/results/te_overlap/te_subfamily_overlap_combined.png")
    outfiles.append(f"{outdir}/results/tracks/ucsc_track.txt")
    outfiles.append(f"{outdir}/results/tracks/igv_track.html")
    outfiles.append(f"{outdir}/PeakCalling_report.pptx")
    outfiles.append(f"{outdir}/PeakCalling_report.xlsx")
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["samples"] = paired_samples + single_samples
    datajson["ip_samples"] = ip_samples
    datajson["input_samples"] = input_samples
    datajson["sample_ip_input_map"] = sample_ip_input_map
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runQuantMS(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ):
    """Prepare input JSON for QuantMS (quantitative proteomics) workflow.
    
    Workflow steps:
    1. Raw conversion (optional, when raw_manifest/raw inputs are provided)
    2. Decoy database generation
    3. Database search engines (Comet, MSGF+, Sage)
    4. PSM rescoring (Percolator)
    5. PSM FDR control
    6. Protein inference (EpiFany)
    7. Protein quantification (ProteomicsLFQ or ProteinQuantifier)
    8. Statistical analysis (MSstats)
    
    Supports TMT, LFQ, DIA, and raw-to-mzML entry points.
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    samples: List[str] = []
    outfiles: List[str] = []

    quantification_method = datajson.get("quantification_method", "lfq")
    search_engine = datajson["Params"]["search_engine"]["engine"]
    for sample_id in samples_info_dict:
        samples.append(sample_id)
        if search_engine== "comet":
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}_comet.idXML")
        elif search_engine == "msgf":
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}_msgf.idXML")
        elif search_engine == "sage":
            outfiles.append(f"{outdir}/common/3_search_engine/{sample_id}/{sample_id}_sage.idXML")
        else:
            raise ValueError(f"Unknown search engine: {search_engine}")
        outfiles.append(f"{outdir}/common/4_psm_rescoring/{sample_id}/{sample_id}_scored.idXML")
        outfiles.append(f"{outdir}/common/5_psm_fdr/{sample_id}/{sample_id}_filtered.idXML")
        outfiles.append(f"{outdir}/common/6_protein_inference/{sample_id}/{sample_id}_protein.idXML")
    if quantification_method == "tmt":
        outfiles.append(f"{outdir}/common/7_quantification/tmt_quantification.mzTab")
    elif quantification_method == "lfq":
        outfiles.append(f"{outdir}/common/7_quantification/lfq_quantification.mzTab")
    elif quantification_method == "dia":
        outfiles.append(f"{outdir}/common/7_quantification/dia_quantification.mzTab")

    if not datajson.get("Params", {}).get("skip_post_msstats", False):
        outfiles.append(f"{outdir}/common/8_msstats/msstats_results.csv")

    datajson["samples"] = samples
    datajson["raw_files"] = raw_files

    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        group_pairs: List[CompareGroupPair],
        indir:str,
        outdir: str,
        raw_files: List[str],
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
    organism = next(iter(Organisms))
    if organism in ["Homo sapiens", "human"]:
        datajson["genome"]["default"] = "GRCh38"
        datajson["Params"]["report"]["genome"] = "GRCh38"
        datajson["Params"]["function"]["species"] = "human"
    elif organism in ["Mus musculus", "mouse"]:
        datajson["genome"]["default"] = "GRCm39"
        datajson["Params"]["report"]["genome"] = "GRCm39"
        datajson["Params"]["function"]["species"] = "mouse"
    else:
        raise ValueError(f"pipeline don't support {organism}, only support human or mouse(Homo sapiens or Mus musculus)")
    datajson["Params"]["report"]["date"] = time.strftime("%Y-%m-%d", time.localtime())
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
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runtRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        meta: str,
        raw_files: List[str],
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
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles

    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runncRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        design_pairs:List[DesignPair],
        indir: str,
        outdir: str,
        raw_files: List[str],
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
    sample_ip_input_map = {}
    ip_samples = []
    input_samples = []
    outfiles = []
    paired_samples = []
    single_samples = []
    layouts = set()
    organisms = set()
    for design_pair in design_pairs:
        sample_ip_input_map[design_pair.exp_sample_id] = design_pair.ctr_sample_id
        ip_samples.append(design_pair.exp_sample_id)
        input_samples.append(design_pair.ctr_sample_id)
    # for ip_sample in ip_samples:
    #     outfiles.append(f"{outdir}/peaks/{ip_sample}/{ip_sample}_peaks.narrowPeak")
    #     outfiles.append(f"{outdir}/peaks/{ip_sample}/{ip_sample}_peaks.xls")

    aligner = datajson.get("Procedure", {}).get("aligner") or "star_3pass"
    if aligner == "star_3pass":
        bam_subdir = "common/3_raw_bam"
    elif aligner == "star_3pass_gene":
        bam_subdir = "common/4_per_gene_bam"
        outfiles.append(f"{outdir}/ncRNAseq_report.pptx")
    else:
        bam_subdir = "common/3_raw_bam"

    for sample_id, sample_info in samples_info_dict.items():
        organisms.add(sample_info.organism)
        outfiles.append(f"{outdir}/{bam_subdir}/{sample_id}/{sample_id}.bam")
        outfiles.append(f"{outdir}/{bam_subdir}/{sample_id}/{sample_id}_tail.csv")
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            layouts.add("PE")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            layouts.add("SE")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    if len(organisms) != 1:
        raise ValueError(f"meta don't support multiple organsim temporarily, please check your meta file, found: {organisms}")
    for organism in organisms:
        if organism in ["Homo sapiens", "human"]:
            datajson["Params"]["track"]["default"] = "GRCh38"
        elif organism in ["Mus musculus", "mouse"]:
            datajson["Params"]["track"]["default"] = "GRCm39"
        else:
            raise ValueError(f"pipeline don't support {organisms.pop()}, only support human or mouse(Homo sapiens or Mus musculus)")
    for layout in layouts:
        if layout == "PE":
            outfiles.append(f"{outdir}/counts/all_paired_featureCounts.tsv")
        elif layout == "SE":
            outfiles.append(f"{outdir}/counts/all_single_featureCounts.tsv")
    all_samples = paired_samples + single_samples
    outfiles.append(f"{outdir}/tracks/igv_track.html")
    outfiles.append(f"{outdir}/tracks/ucsc_track.txt")
    datajson["samples"] = all_samples
    datajson["raw_files"] = raw_files
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    datajson["outfiles"] = outfiles
    datajson["ip_samples"] = ip_samples
    datajson["input_samples"] = input_samples
    datajson["sample_ip_input_map"] = sample_ip_input_map
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json

def runscRNAseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
        cellranger_input_dict: Dict[str, CellrangerInput]
    ):
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir
    datajson["raw_files"] = raw_files
    datajson["cellranger_input_dict"] = {k: v.__dict__ for k, v in cellranger_input_dict.items()}
    outfiles = []
    organisms = set()
    counters = datajson["counters"]
    paired_samples = []
    single_samples = []
    for sample_id, sample_info in samples_info_dict.items():
        organisms.add(sample_info.organism)
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            if "scTE" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_scTE.h5ad")
            elif "cellranger" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_cellranger.h5ad")
            else:
                logger.error(f"Unknown counter type for sample {sample_id}: {counters}")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            if "scTE" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_scTE.h5ad")
            elif "cellranger" in counters:
                outfiles.append(f"{outdir}/common/3_raw_h5ad/{sample_id}/{sample_id}_cellranger.h5ad")
            else:
                logger.error(f"Unknown counter type for sample {sample_id}: {counters}")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    datajson["paired_samples"] = paired_samples
    datajson["single_samples"] = single_samples
    if len(organisms) != 1:
        raise ValueError(f"meta don't support multiple organsim temporarily, please check your meta file, found: {organisms}")
    organism = next(iter(organisms))
    if organism in ["Homo sapiens", "human"]:
        datajson["genome"]["default"] = "GRCh38"
        datajson["Params"]["scTE"]["genome"] = "GRCh38"
        datajson["Params"]["cellranger"]["mkref"]["genome_name"] = "GRCh38"
        datajson["Params"]["cellranger"]["mkref"]["version"] = time.strftime("%Y-%m-%d", time.localtime())
    elif organism in ["Mus musculus", "mouse"]:
        datajson["genome"]["default"] = "GRCm39"
        datajson["Params"]["scTE"]["genome"] = "GRCm39"
        datajson["Params"]["cellranger"]["mkref"]["genome_name"] = "GRCm39"
        datajson["Params"]["cellranger"]["mkref"]["version"] = time.strftime("%Y-%m-%d", time.localtime())
    elif organism in ["Macaca mulatta", "rhesus macaque", "mulatta"]:
        datajson["genome"]["default"] = "Mmul_10"
        datajson["Params"]["scTE"]["genome"] = "Mmul_10"
        datajson["Params"]["cellranger"]["mkref"]["genome_name"] = "Mmul_10"
        datajson["Params"]["cellranger"]["mkref"]["version"] = time.strftime("%Y-%m-%d", time.localtime())
    else:
        raise ValueError(f"pipeline don't support {organism}, only support human, mouse, or macaque(Homo sapiens, Mus musculus, or Macaca mulatta)")
    counters = datajson["counters"]
    aligner = datajson["aligner"]
    if "scTE" in counters:
        # TEs are far less numerous than genes — median ~25 TE/cell
        datajson["Params"]["scanpy"]["scTE"]["qc"]["min_genes"] = 10
        datajson["Params"]["scanpy"]["scTE"]["qc"]["max_genes"] = 3000
        datajson["Params"]["scanpy"]["scTE"]["qc"]["scrublet"] = False
        if aligner == "star":
            datajson["Params"]["scTE"]["cb_tag"] = "CR"
            datajson["Params"]["scTE"]["umi_tag"] = "UR"
        elif aligner == "cellranger":
            datajson["Params"]["scTE"]["cb_tag"] = "CB"
            datajson["Params"]["scTE"]["umi_tag"] = "UB"
        else:
            raise ValueError(f"Unsupported aligner for scTE counter: {aligner}, must be star or cellranger")
    elif "cellranger" in counters:
        # cellranger: keep defaults (min_genes=200, max_genes=6000, scrublet=True)
        datajson["Params"]["scanpy"]["cellranger"]["qc"]["min_genes"] = 200
        datajson["Params"]["scanpy"]["cellranger"]["qc"]["max_genes"] = 6000
        datajson["Params"]["scanpy"]["cellranger"]["qc"]["scrublet"] = True
    else:
        raise ValueError(f"Unsupported counter type: {counters}, must be scTE or cellranger")
    # Build tissue_samples for scanpy downstream
    tissue_samples = {}
    for sid in paired_samples + single_samples:
        tissue = getattr(samples_info_dict.get(sid), "tissue", None) or "unknown"
        tissue_samples.setdefault(tissue, []).append(sid)
    for tissue in tissue_samples.keys():
        for counter in counters:
            outfiles.append(f"{outdir}/common/5_combine_h5ad/{tissue}/{tissue}_{counter}_advanced.h5ad")
    datajson["Params"]["scanpy"]["tissue_samples"] = tissue_samples
    # outfiles.append(f"{outdir}/scRNAseq_report.pptx")
    datajson["outfiles"] = outfiles
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json



def runFiberseq(
        datajson: Dict[str, Any],
        samples_info_dict: Dict[str, SampleInfo],
        indir: str,
        outdir: str,
        raw_files: List[str],
    ) -> str:
    """Prepare input JSON for Fiber-seq workflow.

    Fiber-seq: single-molecule chromatin accessibility sequencing.
    Reference: Stergachis et al., 2020, Science (DOI: 10.1126/science.aaz1646).
    Guide: https://fiberseq.github.io/
    """
    datajson["ROOT_DIR"] = os.path.dirname(__file__)
    datajson["indir"] = indir
    datajson["outdir"] = outdir
    logdir = os.path.join(outdir, "log")
    os.makedirs(logdir, exist_ok=True)
    datajson["logdir"] = logdir

    outfiles = []
    samples = []
    for sample_id, sample_info in samples_info_dict.items():
        samples.append(sample_id)
        # Final outputs: FIRE BAM + extracted BED files
        outfiles.append(f"{outdir}/fiberseq/3_fire/{sample_id}/{sample_id}.fiberseq.fire.bam")
        outfiles.append(f"{outdir}/fiberseq/4_extract/{sample_id}/{sample_id}.m6a.bed.gz")
        outfiles.append(f"{outdir}/fiberseq/4_extract/{sample_id}/{sample_id}.nuc.bed.gz")
        outfiles.append(f"{outdir}/fiberseq/4_extract/{sample_id}/{sample_id}.msp.bed.gz")
        outfiles.append(f"{outdir}/fiberseq/4_extract/{sample_id}/{sample_id}.fire.bed.gz")

    datajson["samples"] = samples
    datajson["raw_files"] = raw_files
    datajson["outfiles"] = outfiles
    instance_json = os.path.join(outdir, "raw.json")
    with open(instance_json, 'w', encoding='utf-8') as wf:
        json.dump(datajson, wf, indent=2, ensure_ascii=False)
    return instance_json
