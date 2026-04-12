shell.prefix("set -x; set -e;")
import os
import logging
import sys
from typing import Dict
SNAKEFILE_FULL_PATH = workflow.snakefile
ROOT_DIR = os.path.dirname(SNAKEFILE_FULL_PATH)
sys.path.append(f"{ROOT_DIR}/src")
configfile: f"{ROOT_DIR}/main.json"
from common.MetaUtil import MetadataUtils
from common.LogUtil import setup_logger
fqdir = config.get('fqdir', 'data/fastq')
logdir = config.get('logdir', 'log')
outdir = config.get('outdir', 'output')
meta = config.get('meta', None)
logger = setup_logger("root", f"{logdir}/workflow.log")
logger.info(f"fastq input directory: {fqdir}; Output directory: {outdir}; Log file: {logdir}/workflow.log")
logger.info(f"Workflow Root directory: {ROOT_DIR}")


metadataUtil = MetadataUtils(
    meta=meta,
    outdir=outdir,
    fastq_dir=fqdir,
)

samples_info_dict, pairs, workdir = metadataUtil.run()


outfiles = []
paired_samples = []
single_samples = []

ip_samples = []
input_samples = []
treated_ip_samples = []
treated_input_samples = []

# function to get the output files for WES analysis
def get_MERIP_outfiles(samples_info_dict:Dict[str, any]):
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}_1.fq.gz")
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}_2.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/dedup/{sample_id}.dedup.bam")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}.fq.gz")
            # outfiles.append(f"{outdir}/hisat2/{sample_id}.bam")
            # outfiles.append(f"{outdir}/igv/{sample_id}.bigwig")
            outfiles.append(f"{outdir}/igv/dedup/{sample_id}.dedup.bam")
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
    outfiles.append(f"{outdir}/exomePeak/exomePeak_diff_peaks.tsv")
    outfiles.append(f"{outdir}/exomePeak/exomePeak_sig_siff_peaks.tsv")
    outfiles.append(f"{outdir}/exomePeak/exomePeak_con_sig_diff_peaks.tsv")
    include: "subworkflow/MERIP.smk"
    return outfiles

get_MERIP_outfiles(samples_info_dict)
logger.info(f"Paired samples: {paired_samples}")
logger.info(f"Single samples: {single_samples}")
logger.info(f"IP samples: {ip_samples}")
logger.info(f"Input samples: {input_samples}")
logger.info(f"Treated IP samples: {treated_ip_samples}")
logger.info(f"Treated Input samples: {treated_input_samples}")
logger.info(f"Final output files for MERIP subworkflow: {outfiles}")

rule all:
    input:
        outfiles