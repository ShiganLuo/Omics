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


def rule_res(name):
    """Resolve per-rule resources from config['resources'] (common.smk
    semantics; KARRseq.smk does not include common.smk, so this is the
    local copy). Missing keys fall back through 'default'; no section at
    all yields the caller's .get() defaults."""
    res_cfg = config.get("resources") or {}
    merged = dict(res_cfg.get("default") or {})
    merged.update(res_cfg.get(name) or {})
    return merged

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
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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
star_params = config.get("Params", {}).get("star", {})


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



# ============================================
# Sex check (coverage on sex chromosomes)
# ============================================
# Wired unconditionally over the alignment BAMs; node.py registers the
# sexcheck outputs (subworkflow.md QC wiring).
sexcheck_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "sample_bams": {
        s: f"{outdir}/bam/{s}/{s}_Aligned.sortedByCoord.out.bam"
        for s in paired_samples + single_samples
    },
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "sex_chroms": (config.get("genome", {}) or {}).get("references", {})
            .get(config.get("genome", {}).get("default"), {})
            .get("sex_chroms"),
        # RNA-mode marker lookup (XIST / Y-linked genes) reads the annotation
        "gtf": (config.get("genome", {}) or {}).get("references", {})
            .get(config.get("genome", {}).get("default"), {})
            .get("gtf")
    },
    "Params": {
        "sexcheck": config.get("Params", {}).get("sexcheck", {})
    },
}
module sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: sexcheck_config
logger.info(f"sexcheck_config: {sexcheck_config}")
use rule sexcheck_idxstats from sexcheck as KARRseq_sexcheck_idxstats
use rule sexcheck from sexcheck as KARRseq_sexcheck


# ============================================
# QC Summary (cross-sample metric aggregation)
# ============================================
# Wired unconditionally; node.py registers the qc_summary outputs (see
# subworkflow.md QC wiring spec). Only stages with qc_summary parsers are
# declared; unregistered stages stay best-effort scans.
qc_summary_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "qc_dirs": {
        "sexcheck": sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "sexcheck": "sexcheck.tsv",
    },
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as KARRseq_qc_summary


rule KARRseq_star_align:
    input:
        fastq = get_alignment_input
    output:
        aligned = outdir + "/bam/{sample_id}/{sample_id}_Aligned.sortedByCoord.out.bam",
        chimeric = outdir + "/bam/{sample_id}/{sample_id}_Chimeric.out.sam"
    log:
        logdir + "/{sample_id}/STAR.log"
    threads: rule_res("KARRseq_star_align").get("threads", 15)
    resources:
        mem_mb = rule_res("KARRseq_star_align").get("mem_mb", 8000),
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
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
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
