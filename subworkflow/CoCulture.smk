shell.prefix("set -x; set -e;")
from snakemake.logging import logger
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
genome_pairs = config.get("genome_pairs", [])
genome = config.get("genome", {})
outfiles = config.get("outfiles", [])
rule all:
    input:
        outfiles

aligner = config.get('Procedure',{}).get('aligner')
SOAPnuke_cofig = {
        "indir": indir,
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