shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
outfiles = config.get("outfiles", [])
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

sample_data = config.get("sample_data", "")

rule all:
    input:
        outfiles

fumitools_config = {
    "indir": indir,
    "outdir": f"{outdir}/fastq/umi_fastq",
    "logdir": logdir,
    "paired_samples": config.get("paired_samples", []),
    "ROOT_DIR": ROOT_DIR,
    "Procedure": {
        "fumitools": config.get("Procedure", {}).get("fumitools")
    },
    "Params": {
        "umi_length": config.get("Params", {}).get("fumitools", {}).get("umi_length", 12),
        "tag_umi": config.get("Params", {}).get("fumitools", {}).get("tag_umi", False),
    }
}
logger.info(f"fumitools_config: {fumitools_config}")
module fumitools:
    snakefile: "../modules/fumitools/fumitools.smk"
    config: fumitools_config
use rule fumitools_copy_umi from fumitools as tRNAseq_fumitools_copy_umi
use rule fumitools_copy_umi_single from fumitools as tRNAseq_fumitools_copy_umi_single

cutadapt_config = {
            "indir": fumitools_config["outdir"],
            "outdir":  f"{outdir}/fastq/trimmed_fastq",
            "logdir": logdir,
            "mode": "UMI",
            "Procedure": {
                "trim_galore": config.get('Procedure',{}).get('trim_galore')
            }
}
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
logger.info(f"cutadapt_config: {cutadapt_config}")
use rule trimming_Paired from cutadapt as RNAseq_trimming_Paireds
use rule trimming_Single from cutadapt as RNAseq_trimming_Single

# Module config dict
mimseq_config = {
    "indir": indir,
    "outdir": f"{outdir}/mimseq",
    "logdir": logdir,
    "samples": samples,
    "sample_data": sample_data,
    "ROOT_DIR": ROOT_DIR,
    "Procedure": {
    },
    "Params": {
        "mimseq": config.get("Params", {}).get("mimseq", {})
    },
    "genome": {
        "trnas": config.get("genome", {}).get("trnas", ""),
        "trnaout": config.get("genome", {}).get("trnaout", ""),
        "mito_trnas": config.get("genome", {}).get("mito_trnas", ""),
        "plastid_trnas": config.get("genome", {}).get("plastid_trnas", ""),
    }
}
logger.info(f"mimseq_config: {mimseq_config}")

module mimseq:
    snakefile: "../modules/mimseq/mimseq.smk"
    config: mimseq_config

use rule mimseq_run from mimseq as tRNAseq_mimseq_run
use rule mimseq_result from mimseq as tRNAseq_mimseq_result
