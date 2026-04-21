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
configfile: f"{ROOT_DIR}/config/CoCulture.json"
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
single_sample_genome_pairs = []
paired_sample_genome_pairs = []
def get_CoCultrue_outfiles(samples_info_dict:Dict[str,any]):
    for sample_id, sample_info in samples_info_dict.items():
        if sample_info.layout == "PE":
            paired_samples.append(sample_id)
            paired_sample_genome_pairs.append((sample_id, sample_info.organism))
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_1.fq.gz")
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}_2.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
        elif sample_info.layout == "SE":
            single_samples.append(sample_id)
            single_sample_genome_pairs.append((sample_id, sample_info.organism))
            outfiles.append(f"{outdir}/SOAPnuke/{sample_id}.single.fq.gz")
            outfiles.append(f"{outdir}/hisat2/GRCm39/{sample_id}.bam")
            outfiles.append(f"{outdir}/hisat2/GRCh38/{sample_id}.bam")
            # outfiles.append(f"{outdir}/disambiguate/{sample_id}/{sample_id}.disambiguatedSpecies_GRCm39.bam")
            # outfiles.append(f"{outdir}/disambiguate/{sample_id}/{sample_id}.disambiguatedSpecies_GRCh38.bam")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCm39/all_TEcount.tsv")
            outfiles.append(f"{outdir}/TEtranscripts/TEcount/GRCh38/all_TEcount.tsv")
        else:
            logger.error(f"Unknown layout type for sample {sample_id}: {sample_info.layout}")
    logger.info(f"Paired sample genome pairs: {paired_sample_genome_pairs}")
    logger.info(f"Single sample genome pairs: {single_sample_genome_pairs}")
    logger.info(f" raw_fastq_dir: {raw_fastq_dir}")
get_CoCultrue_outfiles(samples_info_dict)
logger.info(f"Parser control Variables:")
logger.info(f"Paired samples: {paired_samples}")
logger.info(f"Single samples: {single_samples}")
logger.info(f"Final output files: {outfiles}")
logger.info(f"Single sample genome pairs: {single_sample_genome_pairs}")
logger.info(f"Paired sample genome pairs: {paired_sample_genome_pairs}")
rule all:
    input:
        outfiles


aligner = config.get('Procedure',{}).get('aligner')
SOAPnuke_cofig = {
        "indir": raw_fastq_dir,
        "outdir":  f"{outdir}/SOAPnuke",
        "logdir": logdir
}
module SOAPnuke:
    snakefile: "../modules/SOAPnuke/SOAPnuke.smk"
    config: SOAPnuke_cofig
logger.info(f"SOAPnuke config: {SOAPnuke_cofig}")
use rule soapnuke_filter_paired from SOAPnuke as CoCulture_soapnuke_filter_paired
use rule soapnuke_filter_single from SOAPnuke as CoCulture_soapnuke_filter_single

    
hisat2_config = {
        "indir": SOAPnuke_cofig["outdir"],
        "outdir":  f"{outdir}/hisat2",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "hisat2": config.get('Procedure',{}).get('hisat2')
        },
        "genome": config.get("genome", {})
    }

module hisat2:
    snakefile: "../modules/hisat2/polygenomes/hisat2.smk"
    config: hisat2_config
logger.info(f"hisat2_config: {hisat2_config}")

use rule hisat2_align from hisat2 as CoCulture_hisat2_align
use rule hisat2_index from hisat2 as CoCulture_hisat2_index

disambiguate_config = {
        "indir": hisat2_config["outdir"],
        "outdir": f"{outdir}/disambiguate",
        "logdir": logdir,
        "genome_pairs": config.get("genome_pairs", []),
        "bam": {
            "aligner": config.get('Procedure',{}).get('aligner') or 'hisat2'
        },
        "Procedure": {
            "ngs_disambiguate": config.get('Procedure',{}).get('ngs_disambiguate') or 'ngs_disambiguate'
        }
    }
module disambiguate:
    snakefile: "../modules/disambiguate/disambiguate.smk"
    config: disambiguate_config
logger.info(f"disambiguate config: {disambiguate_config}")
use rule ngs_disambiguate from disambiguate as CoCulture_ngs_disambiguate
use rule disambiguate_sort_rename from disambiguate as CoCulture_disambiguate_sort_rename

TEtranscripts_config = {
        "indir": disambiguate_config["outdir"],
        "outdir": f"{outdir}/TEtranscripts",
        "logdir": logdir,
        "ROOT_DIR": ROOT_DIR,
        "genome_pairs": disambiguate_config["genome_pairs"],
        "single_samples": single_samples,
        "paired_samples": paired_samples,
        "genome": config.get("genome", {}),
        "Procedure": {
            "TEcount": config.get('Procedure',{}).get('TEcount') or 'TEcount',
            "TElocal": config.get('Procedure',{}).get('TElocal') or 'TElocal'
        }
    }
module TEtranscripts:
    snakefile: "../modules/TEtranscripts/disambiguate/TEtranscripts.smk"
    config: TEtranscripts_config
logger.info(f"TEtranscripts config: {TEtranscripts_config}")
use rule * from TEtranscripts as CoCulture_*