shell.prefix("set -x; set -e;")
from snakemake.logging import logger
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])

rule all:
    input:
        outfiles

fastqc_raw_config = {
        "indir": indir,
        "outdir":  f"{outdir}/fastqc/raw",
        "logdir": logdir,
        "log_suffix": "raw.txt",
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "fastqc": config.get("Procedure", {}).get("fastqc") or "fastqc"
        }
    }
module fastqc_raw:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_raw_config
logger.info(f"fastqc_raw_config: {fastqc_raw_config}")
use rule fastqc from fastqc_raw as WES_fastqc_raw

cutadapt_config = {
        "indir": indir,
        "outdir": f"{outdir}/cutadapt",
        "logdir": logdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": {
                "quality": config.get('Params',{}).get("trim_galore", {}).get('quality')
            }
        },
    }

module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
logger.info(f"Cutadapt parameters: {cutadapt_config}")
use rule trimming_Paired from cutadapt as WES_trimming_Paired

fastqc_trimmed_config = {
        "indir": cutadapt_config["outdir"],
        "outdir":  f"{outdir}/fastqc/trimmed",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "log_suffix": "trimmed.txt",
        "Procedure": {
            "fastqc": config.get("Procedure", {}).get("fastqc")
        }
    }
module fastqc_trimmed:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_trimmed_config
logger.info(f"fastqc_trimmed_config: {fastqc_trimmed_config}")
use rule fastqc from fastqc_trimmed as WES_fastqc_trimmed

bwa_mem2_confg = {
    "indir": cutadapt_config["outdir"],
    "outdir":  f"{outdir}/bwa_mem2",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "bwaMem2": config.get("Procedure",{}).get("bwaMem2"),
        "samtools": config.get("Procedure",{}).get("samtools")
    },
    "genome": {
        "fasta": config.get("genome",{}).get("fasta"),
        "index_prefix": config.get("genome",{}).get("bwaMem2_index_prefix")
    }
}

module bwa_mem2:
    snakefile: "../modules/bwa-mem2/bwa-mem2.smk"
    config: bwa_mem2_confg
logger.info(f"BWA-MEM2 parameters: {bwa_mem2_confg}")
use rule bwaMem2_index from bwa_mem2 as WES_bwaMem2_index
use rule bwaMem2_alignment from bwa_mem2 as WES_bwaMem2_alignment

samtools_config = {
    "indir": bwa_mem2_confg["outdir"],
    "outdir": bwa_mem2_confg["outdir"],
    "logdir": logdir
    "Procedure": {
        "samtools": config.get("Procedure",{}).get("samtools")
    },
    "Params": {
        "samtools": {
            "onlykeep_properpair": config.get("Params",{}).get("samtools",{}).get("onlykeep_properpair")
        }
    }
}

module samtools:
    snakefile: "../modules/samtools/samtools.smk"
    config: samtools_config
logger.info(f"Samtools parameters: {samtools_config}")
use rule bam_sort from samtools as WES_bam_sort


gatk_prepare_config = {
    "indir": bwa_mem2_confg["outdir"],
    "outdir": f"{outdir}/gatk",
    "logdir": logdir
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
        "samtools": config.get("Procedure", {}).get("samtools"),
    },
    "addReadsGroup": {
        "RGLB": config.get("addReadsGroup", {}).get("RGLB"),
        "RGPL": config.get("addReadsGroup", {}).get("RGPL"),
        "RGPU": config.get("addReadsGroup", {}).get("RGPU")
    },

    "fasta": fasta,
    "tmp_dir": config.get("tmp_dir")
}
module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
use rule gatk_index from gatk_prepare as WES_gatk_index
use rule addReadsGroup from gatk_prepare as WES_addReadsGroup
use rule MarkDuplicates from gatk_prepare as WES_MarkDuplicates

gatk_somatic_config = {
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk"),
    },
    "outdir": outdir,
    "fasta": fasta,
    "tmp_dir": config.get("tmp_dir"),
    "mutect2_parameters": config.get("mutect2_parameters")
}
module gatk_somatic:
    snakefile: "../modules/gatk/gatk_somatic/gatk_somatic.smk"
    config: gatk_somatic_config
use rule somaticMutect2 from gatk_somatic as WES_somaticMutect2
use rule gatk_index from gatk_prepare as WES_gatk_index


