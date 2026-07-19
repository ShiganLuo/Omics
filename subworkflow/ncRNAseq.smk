shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
all_samples = config.get("samples", [])
outfiles = config.get("outfiles", [])
aligner = config.get("Procedure", {}).get("aligner") or "star"

rule all:
    input:
        outfiles

# ── 1. Trim (cutadapt / trim_galore) ─────────────────────────────────────────
cutadapt_config = {
    "indir": indir,
    "outdir": f"{outdir}/ncRNAseq/cutadapt",
    "logdir": logdir,
    "Procedure": {
        "trim_galore": config.get("Procedure", {}).get("trim_galore")
    }
}
logger.info(f"cutadapt_config: {cutadapt_config}")
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
use rule trimming_Paired from cutadapt as ncRNAseq_trimming_Paired
use rule trimming_Single from cutadapt as ncRNAseq_trimming_Single

# ── 2. Align (hisat2-ncRNAseq or star) ───────────────────────────────────────
if aligner == "hisat2":
    hisat2_config = {
        "indir": cutadapt_config["outdir"],
        "outdir": f"{outdir}/ncRNAseq/bam",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "hisat2": config.get("Procedure", {}).get("hisat2"),
            "hisat2-build": config.get("Procedure", {}).get("hisat2-build")
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta"),
            "hisat2_index_prefix": config.get("genome", {}).get("hisat2_index_prefix")
        }
    }
    logger.info(f"hisat2_config: {hisat2_config}")
    module hisat2:
        snakefile: "../modules/hisat2/ncRNAseq/hisat2.smk"
        config: hisat2_config
    use rule hisat2_align_ncRNAseq_single from hisat2 as ncRNAseq_hisat2_align
elif aligner == "star":
    star_config = {
        "indir": cutadapt_config["outdir"],
        "outdir": f"{outdir}/ncRNAseq/bam",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "STAR": config.get("Procedure", {}).get("STAR")
        },
        "Params": {
            "STAR": {
                "genomeLoad": config.get("Params", {}).get("STAR", {}).get("genomeLoad") or "LoadAndRemove",
                "limitBAMsortRAM": config.get("Params", {}).get("STAR", {}).get("limitBAMsortRAM") or 20000000000,
                "outReadsUnmapped": config.get("Params", {}).get("STAR", {}).get("outReadsUnmapped") or "Fastx",
                "outFilterMultimapNmax": config.get("Params", {}).get("STAR", {}).get("outFilterMultimapNmax") or 99999,
                "outFilterMismatchNoverLmax": config.get("Params", {}).get("STAR", {}).get("outFilterMismatchNoverLmax") or 0.1,
                "outFilterMatchNminOverLread": config.get("Params", {}).get("STAR", {}).get("outFilterMatchNminOverLread") or 0.66,
                "alignSJoverhangMin": config.get("Params", {}).get("STAR", {}).get("alignSJoverhangMin") or 999,
                "alignSJDBoverhangMin": config.get("Params", {}).get("STAR", {}).get("alignSJDBoverhangMin") or 999
            }
        },
        "genome": {
            "fasta": config.get("genome", {}).get("fasta"),
            "gtf": config.get("genome", {}).get("gtf"),
            "index_dir": config.get("genome", {}).get("star_index_dir")
        }
    }
    logger.info(f"star_config: {star_config}")
    module star:
        snakefile: "../modules/star/star.smk"
        config: star_config
    use rule star_align from star as ncRNAseq_star_align
else:
    raise ValueError(f"Unsupported aligner: {aligner}. Please choose 'hisat2' or 'star'.")

# ── 3. Quantify (featureCounts) ──────────────────────────────────────────────
featureCounts_config = {
    "indir": hisat2_config["outdir"] if aligner == "hisat2" else star_config["outdir"],
    "outdir": f"{outdir}/ncRNAseq/counts",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "featureCounts": config.get("Procedure", {}).get("featureCounts")
    },
    "genome": {
        "gtf": config.get("genome", {}).get("gtf")
    }
}
logger.info(f"featureCounts_config: {featureCounts_config}")
module featureCounts:
    snakefile: "../modules/featureCounts/featureCounts.smk"
    config: featureCounts_config
use rule featureCounts_single_noMultiple from featureCounts as ncRNAseq_featureCounts_single
use rule featureCounts_paired_noMultiple from featureCounts as ncRNAseq_featureCounts_paired
use rule featureCounts_result from featureCounts as ncRNAseq_featureCounts_result

