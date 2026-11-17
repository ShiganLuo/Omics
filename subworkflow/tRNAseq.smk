shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "log")
outfiles = config.get("outfiles", [])
samples = config.get("Params", {}).get("workflow", {}).get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

sample_data = config.get("Params", {}).get("workflow", {}).get("sample_data", "")

# Genome reference resolution — nested default + references pattern
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
    "paired_samples": config.get("Params", {}).get("workflow", {}).get("paired_samples", []),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "Procedure": {
        "fumitools": config.get("Procedure", {}).get("fumitools")
    },
    "Params": {
        # 参数精确到 rule:整段透传(Params.fumitools.<rule>.<key>,rule 子段
        # 透明传递;模块内置默认同原转发默认)
        "fumitools": config.get("Params", {}).get("fumitools", {}),
    }
}
logger.info(f"fumitools_config: {fumitools_config}")
module fumitools:
    snakefile: "../modules/fumitools/fumitools.smk"
    config: fumitools_config
use rule fumitools_copy_umi from fumitools as tRNAseq_fumitools_copy_umi
use rule fumitools_copy_umi_single from fumitools as tRNAseq_fumitools_copy_umi_single

cutadapt_config = {

    "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
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
    "resources": config.get("resources", {}),
    "indir": cutadapt_config["outdir"],
    "outdir": f"{outdir}/mimseq",
    "logdir": logdir,
    "meta": config.get("Params", {}).get("workflow", {}).get("meta", ""),
    "samples": samples,
    "sample_data": sample_data,
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    # data_dir 精确到 rule:Params.mimseq.mimseq_tRNAtools.data_dir(config
    # 全量路径读取,无旧扁平键回退;顶层 data_dir 仍是接线键,透传给模块顶层)
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
