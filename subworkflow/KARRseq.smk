"""
KARRseq Subworkflow

Kethoxal-Assisted RNA-RNA interaction sequencing (KARR-seq) workflow.
Detects RNA-RNA interactions via chimeric read analysis.

Workflow steps:
1. SeqPrep: Merge paired-end reads (PE only)
2. STAR: Align reads (merged PE or raw SE) + output chimeric reads
3. KARRseq: Extract chimeric pairs, remove duplicates, process ligation events

Usage:
    snakemake -s subworkflow/KARRseq.smk --configfile config/KARRseq.json --use-conda --cores 30

Reference: https://github.com/ouyang-lab/KARR-seq.git
"""
shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os

indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
samples = paired_samples + single_samples
ROOT_DIR = config.get("ROOT_DIR", os.path.dirname(os.path.dirname(__file__)))


rule all:
    input:
        outfiles


# =================
# Module 1: SeqPrep (merge paired-end reads)
# =================

if paired_samples:
    seqprep_config = {
        "indir": indir,
        "outdir": f"{outdir}/merge",
        "logdir": logdir,
        "samples": paired_samples,
        "Procedure": {
            "SeqPrep": config.get("Procedure", {}).get("SeqPrep") or "SeqPrep",
            "seqkit": config.get("Procedure", {}).get("seqkit") or "seqkit"
        }
    }
    module seqprep:
        snakefile: "../modules/seqprep/seqprep.smk"
        config: seqprep_config
    logger.info(f"seqprep_config: {seqprep_config}")
    use rule seqprep_merge from seqprep as KARRseq_seqprep_merge


# =================
# Module 2: STAR Alignment
# =================

# STAR alignment with KARRseq-specific parameters
star_index_dir = config.get("genome", {}).get("star_index_dir")
star_params = config.get("Params", {}).get("STAR", {})


def get_alignment_input(wildcards):
    """Dynamically determine alignment input: merged PE or raw SE."""
    if wildcards.sample_id in paired_samples:
        merged = f"{outdir}/merge/{wildcards.sample_id}/{wildcards.sample_id}_merge.fastq.gz"
        return [merged]
    elif wildcards.sample_id in single_samples:
        single = f"{indir}/{wildcards.sample_id}.fastq.gz"
        return [single]
    else:
        raise ValueError(f"Sample {wildcards.sample_id} not in paired_samples or single_samples")


rule KARRseq_star_align:
    input:
        fastq = get_alignment_input
    output:
        aligned = outdir + "/bam/{sample_id}/{sample_id}_Aligned.sortedByCoord.out.bam",
        chimeric = outdir + "/bam/{sample_id}/{sample_id}_Chimeric.out.sam"
    log:
        logdir + "/{sample_id}/STAR.log"
    threads: 15
    params:
        STAR = config.get("Procedure", {}).get("STAR") or "STAR",
        index = star_index_dir,
        prefix = outdir + "/bam/{sample_id}/{sample_id}",
        outFilterMultimapNmax = star_params.get("outFilterMultimapNmax", 100),
        outSAMattributes = star_params.get("outSAMattributes", "All"),
        alignIntronMin = star_params.get("alignIntronMin", 1),
        scoreGapNoncan = star_params.get("scoreGapNoncan", -4),
        scoreGapATAC = star_params.get("scoreGapATAC", -4),
        chimSegmentMin = star_params.get("chimSegmentMin", 15),
        chimJunctionOverhangMin = star_params.get("chimJunctionOverhangMin", 15),
        limitOutSJcollapsed = star_params.get("limitOutSJcollapsed", 10000000),
        limitIObufferSize = star_params.get("limitIObufferSize", 1500000000),
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        import time
        current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        logger.info(f"Start STAR alignment for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/STAR_align_{current_time}.sh")
        input_params = " ".join(input.fastq)
        cmd = [
            params.STAR,
            "--runMode", "alignReads",
            "--runThreadN", str(threads),
            "--genomeDir", params.index,
            "--readFilesIn", input_params,
            "--readFilesCommand", "zcat",
            "--outFileNamePrefix", f"{params.prefix}_",
            "--outReadsUnmapped", "Fastq",
            "--outFilterMultimapNmax", str(params.outFilterMultimapNmax),
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMattributes", params.outSAMattributes,
            "--alignIntronMin", str(params.alignIntronMin),
            "--scoreGapNoncan", str(params.scoreGapNoncan),
            "--scoreGapATAC", str(params.scoreGapATAC),
            "--chimSegmentMin", str(params.chimSegmentMin),
            "--chimJunctionOverhangMin", str(params.chimJunctionOverhangMin),
            "--limitOutSJcollapsed", str(params.limitOutSJcollapsed),
            "--limitIObufferSize", str(params.limitIObufferSize),
        ]
        cmd_index = [params.samtools, "index", output.aligned]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
            f.write(" ".join(cmd_index) + "\n")
        shell(f"bash {script} > {log} 2>&1")


# =================
# Module 3: KARRseq Chimeric Processing
# =================

karrseq_config = {
    "indir": f"{outdir}/bam",
    "outdir": f"{outdir}/chimeric",
    "logdir": logdir,
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools") or "samtools"
    },
    "Params": {
        "karrseq": config.get("Params", {}).get("karrseq", {})
    }
}
module karrseq:
    snakefile: "../modules/karrseq/karrseq.smk"
    config: karrseq_config
logger.info(f"karrseq_config: {karrseq_config}")
use rule karrseq_chimeric_to_pairs from karrseq as KARRseq_chimeric_to_pairs
use rule karrseq_remove_duplicates from karrseq as KARRseq_remove_duplicates
use rule karrseq_ligation from karrseq as KARRseq_ligation
