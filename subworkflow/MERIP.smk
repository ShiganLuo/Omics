import logging
logger = logging.getLogger(__name__)

cutadapt_config = {
        "indir": workdir,
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
            "fasta": config.get('genome',{}).get('fasta')
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
            "bamCoverage": config.get('Procedure',{}).get('deepTools')
        }
    }
module igv:
    snakefile: "../modules/igv/igv.smk"
    config: igv_config
logger.info(f"igv_config: {igv_config}")
use rule dedup from igv as MERIP_dedup
use rule wig from igv as MERIP_wig
