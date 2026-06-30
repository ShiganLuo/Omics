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
