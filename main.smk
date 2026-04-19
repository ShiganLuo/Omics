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

#exomePeak
ip_samples = []
input_samples = []
treated_ip_samples = []
treated_input_samples = []

# TEtranscripts
single_sample_genome_pairs = []
paired_sample_genome_pairs = []
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
            # outfiles.append(f"{outdir}/cutadapt/{sample_id}.single.fq.gz")
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
    outfiles.append(f"{outdir}/exomePeak/sig_diff_peak.xls")
    include: "subworkflow/MERIP.smk"
    return outfiles

# get_MERIP_outfiles(samples_info_dict)

def get_CoCultrue_outfiles(samples_info_dict:Dict[str, any]):
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            paired_sample_genome_pairs.append((sample_id, sample_info.organism))
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/disambiguate/{sample_id}/{sample_id}.disambiguatedSpeciesA.bam")
            outfiles.append(f"{outdir}/disambiguate/{sample_id}/{sample_id}.disambiguatedSpeciesB.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            single_sample_genome_pairs.append((sample_id, sample_info.organism))
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/disambiguate/{sample_id}/{sample_id}.disambiguatedSpeciesA.bam")
            outfiles.append(f"{outdir}/disambiguate/{sample_id}/{sample_id}.disambiguatedSpeciesB.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    include: "subworkflow/CoCulture.smk"
    return outfiles
get_CoCultrue_outfiles(samples_info_dict)
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