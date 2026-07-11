shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
samples = config.get("samples", [])
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})

rule all:
    input:
        outfiles

fastqc_raw_config = {
        "indir": indir,
        "outdir":  f"{outdir}/QC/1_raw_fastqc",
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
use rule fastqc from fastqc_raw as PeakCalling_fastqc_raw

trim_galore_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": indir,
        "outdir": f"{outdir}/common/2_trimmed_fastq",
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

module trim_galore:
    snakefile: "../modules/trim-galore/trim-galore.smk"
    config: trim_galore_config
logger.info(f"TrimGalore parameters: {trim_galore_config}")
use rule trimming_Paired from trim_galore as PeakCalling_trimming_Paired
use rule trimming_Single from trim_galore as PeakCalling_trimming_Single

fastqc_trimmed_config = {
        "indir": trim_galore_config["outdir"],
        "outdir":  f"{outdir}/QC/2_trimmed_fastqc",
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
use rule fastqc from fastqc_trimmed as PeakCalling_fastqc_trimmed


bowtie2_config = {
    "indir": trim_galore_config["outdir"],
    "outdir": f"{outdir}/common/3_raw_bam",
    "logdir": logdir,
    "Procedure": {
        "bowtie2-build": config.get("Procedure", {}).get("bowtie2-build"),
        "bowtie2": config.get("Procedure", {}).get("bowtie2")
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta"),
        "index_prefix": config.get("genome", {}).get("bowtie2_index_prefix")
    }
}
module bowtie2:
    snakefile: "../modules/bowtie2/bowtie2.smk"
    config: bowtie2_config
logger.info(f"bowtie2_config: {bowtie2_config}")
use rule bowtie2_index from bowtie2 as PeakCalling_bowtie2_index
use rule bowtie2_align_paired from bowtie2 as PeakCalling_bowtie2_align_paired
use rule bowtie2_align_single from bowtie2 as PeakCalling_bowtie2_align_single

macs3_config = {
    "indir": bowtie2_config["outdir"],
    "outdir": f"{outdir}/peaks",
    "logdir": logdir,
    "samples": ip_samples,
    "ip_samples": ip_samples,
    "input_samples": input_samples,
    "sample_ip_input_map": sample_ip_input_map,
    "Procedure": {
        "macs3": config.get("Procedure", {}).get("macs3")
    },
    "Params": {
        "macs3": config.get("Params", {}).get("macs3", {})
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module macs3:
    snakefile: "../modules/macs3/macs3.smk"
    config: macs3_config
logger.info(f"macs3_config: {macs3_config}")
use rule macs3_callpeak from macs3 as PeakCalling_macs3_callpeak
