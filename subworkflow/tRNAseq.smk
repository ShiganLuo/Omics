shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "log")
outfiles = config.get("outfiles", [])
samples = config.get("sample_params", {}).get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

sample_data = config.get("sample_params", {}).get("sample_data", "")
trimmer = config.get('Params', {}).get('workflow', {}).get('trimmer') or "trim_galore"

genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(genome, {}) if genome else {}

rule all:
    input:
        outfiles

fumitools_config = {
    "resources": config.get("resources", {}),
    "indir": indir,
    "outdir": f"{outdir}/fastq/umi_fastq",
    "logdir": logdir,
    "paired_samples": config.get("sample_params", {}).get("paired_samples", []),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "Procedure": {
        "fumitools": config.get("Procedure", {}).get("fumitools")
    },
    "Params": {
        "fumitools": config.get("Params", {}).get("fumitools", {}),
    }
}
logger.info(f"fumitools_config: {fumitools_config}")
module fumitools:
    snakefile: "../modules/fumitools/fumitools.smk"
    config: fumitools_config
use rule fumitools_copy_umi from fumitools as tRNAseq_fumitools_copy_umi
use rule fumitools_copy_umi_single from fumitools as tRNAseq_fumitools_copy_umi_single

trimmed_dir = f"{outdir}/fastq/trimmed_fastq"
if trimmer == "cutadapt":
    cutadapt_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": fumitools_config["outdir"],
        "outdir": trimmed_dir,
        "logdir": logdir,
        "mode": "UMI",
        "Procedure": {
            "cutadapt": config.get('Procedure',{}).get('cutadapt')
        },
        "Params": {
            "cutadapt": config.get('Params',{}).get("cutadapt", {}),
        }
    }
    module cutadapt:
        snakefile: "../modules/cutadapt/cutadapt.smk"
        config: cutadapt_config
    logger.info(f"cutadapt_config: {cutadapt_config}")
    use rule trimming_Paired from cutadapt as tRNAseq_cutadapt_Paired
    use rule trimming_Single from cutadapt as tRNAseq_cutadapt_Single
elif trimmer == "trimmomatic":
    raise ValueError("trimmomatic 模块不支持 UMI mode,tRNAseq 的 Params.workflow.trimmer 只能是 'cutadapt' 或 'trim_galore'")
else:
    logger.info("use default trimmer: trim_galore")
    trim_galore_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": fumitools_config["outdir"],
        "outdir": trimmed_dir,
        "logdir": logdir,
        "mode": "UMI",
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": config.get('Params',{}).get("trim_galore", {}),
        }
    }
    module trim_galore:
        snakefile: "../modules/trim-galore/trim-galore.smk"
        config: trim_galore_config
    logger.info(f"trim_galore_config: {trim_galore_config}")
    use rule trimming_Paired from trim_galore as tRNAseq_trim_galore_Paired
    use rule trimming_Single from trim_galore as tRNAseq_trim_galore_Single

mimseq_config = {
    "resources": config.get("resources", {}),
    "indir": trimmed_dir,
    "outdir": f"{outdir}/mimseq",
    "logdir": logdir,
    "meta": config.get("Params", {}).get("workflow", {}).get("meta", ""),
    "samples": samples,
    "sample_data": sample_data,
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "data_dir": ((config.get("Params", {}).get("mimseq", {}) or {}).get("mimseq_tRNAtools", {}) or {}).get("data_dir") or "",
    "Procedure": {
    },
    "Params": {
        "mimseq": config.get("Params", {}).get("mimseq", {})
    },
    "genome": {
        "trnas": genome_ref.get("trnas", ""),
        "trnaout": genome_ref.get("trnaout", ""),
        "mito_trnas": genome_ref.get("mito_trnas", ""),
    }
}
logger.info(f"mimseq_config: {mimseq_config}")

module mimseq:
    snakefile: "../modules/mimseq/mimseq.smk"
    config: mimseq_config

use rule mimseq_all from mimseq as tRNAseq_mimseq_run
use rule mimseq_result from mimseq as tRNAseq_mimseq_result
use rule mimseq_prepare_sample_data from mimseq as tRNAseq_mimseq_prepare_sample_data
use rule mimseq_tRNAtools from mimseq as tRNAseq_mimseq_tRNAtools
use rule mimseq_align from mimseq as tRNAseq_mimseq_align
use rule mimseq_clusters from mimseq as tRNAseq_mimseq_clusters
use rule mimseq_mods from mimseq as tRNAseq_mimseq_mods
use rule mimseq_coverage from mimseq as tRNAseq_mimseq_coverage
use rule mimseq_deseq from mimseq as tRNAseq_mimseq_deseq
