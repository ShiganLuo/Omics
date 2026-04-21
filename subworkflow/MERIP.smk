shell.prefix("set -x; set -e;")
import os
import logging
import sys
from typing import Dict
SNAKEFILE_FULL_PATH = workflow.snakefile
ROOT_DIR = os.path.dirname(os.path.dirname(SNAKEFILE_FULL_PATH))
sys.path.append(f"{ROOT_DIR}/src")
from common.MetaUtil import MetadataUtils
from common.LogUtil import setup_logger
configfile: f"{ROOT_DIR}/config/MERIP.json"
fqdir = config.get('fqdir', 'data/fastq')
logdir = config.get('logdir', 'log')
outdir = config.get('outdir', 'output')
meta = config.get('meta', None)
logger = setup_logger("root", f"{logdir}/workflow.log")
logger.info(f"Workflow Root directory: {ROOT_DIR}")
logger.info(f"fastq input directory: {fqdir}; Output directory: {outdir}; Log file: {logdir}/workflow.log; meta file: {meta}")
metadataUtil = MetadataUtils(
    meta=meta,
    outdir=outdir,
    fastq_dir=fqdir,
)
samples_info_dict, pairs, raw_fastq_dir = metadataUtil.run()
# Define global variables
outfiles = []
paired_samples = []
single_samples = []
#exomePeak
ip_samples = []
input_samples = []
treated_ip_samples = []
treated_input_samples = []
def get_MERIP_outfiles(samples_info_dict:Dict[str,any]):
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
    return outfiles
get_MERIP_outfiles(samples_info_dict)
logger.info(f"Parser control Variables:")
logger.info(f"Paired samples: {paired_samples}")
logger.info(f"Single samples: {single_samples}")
logger.info(f"Final output files: {outfiles}")
logger.info(f"IP samples: {ip_samples}")
logger.info(f"Input samples: {input_samples}")
logger.info(f"Treated IP samples: {treated_ip_samples}")
logger.info(f"Treated Input samples: {treated_input_samples}")
rule all:
    input:
        outfiles
cutadapt_config = {
        "indir": raw_fastq_dir,
        "outdir":  f"{outdir}/cutadapt",
        "logdir": logdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        }
    }
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
logger.info(f"cutadapt_config: {cutadapt_config}")
use rule trimming_Paired from cutadapt as MERIP_trimming_Paired

hisat2_config = {
        "indir": cutadapt_config["outdir"],
        "outdir":  f"{outdir}/hisat2",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "hisat2": config.get('Procedure',{}).get('hisat2')
        },
        "genome": {
            "fasta": config.get('genome',{}).get('fasta'),
            "index_prefix": config.get('genome',{}).get('hisat2_index_prefix'),
        }
    }
module hisat2:
    snakefile: "../modules/hisat2/hisat2.smk"
    config: hisat2_config
logger.info(f"hisat2_config: {hisat2_config}")

use rule hisat2_align from hisat2 as MERIP_hisat2_align
use rule hisat2_index from hisat2 as MERIP_hisat2_index

igv_config = {
        "indir": hisat2_config["outdir"],
        "outdir":  f"{outdir}/igv",
        "logdir": logdir,
        "Procedure": {
            "samtools": config.get('Procedure',{}).get('samtools'),
            "deepTools": config.get('Procedure',{}).get('deepTools')
        }
    }
module igv:
    snakefile: "../modules/igv/igv.smk"
    config: igv_config
logger.info(f"igv_config: {igv_config}")
use rule dedup from igv as MERIP_dedup
use rule wig from igv as MERIP_wig

exomePeak_config = {
        "indir": igv_config["outdir"] + "/dedup",
        "outdir": f"{outdir}/exomePeak",
        "logdir": logdir,
        "gtf": config.get("gtf",""),
        "ip_samples": ip_samples,
        "input_samples": input_samples,
        "treated_ip_samples": treated_ip_samples,
        "treated_input_samples": treated_input_samples,
    }
module exomePeak:
    snakefile: "../modules/exomePeak/exomePeak.smk"
    config: exomePeak_config
logger.info(f"exomePeak_config: {exomePeak_config}")
use rule diff_exomePeak from exomePeak as MERIP_diff_exomePeak
use rule call_exomePeak from exomePeak as MERIP_call_exomePeak
