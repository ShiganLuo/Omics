shell.prefix("set -x; set -e;")
from snakemake.logging import logger
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
aligner = config.get('aligner', 'hisat2')
outfiles = config.get("outfiles", [])
rule all:
    input:
        outfiles
cutadapt_config = {
        "indir": indir,
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
use rule trimming_Paired from cutadapt as RNAseq_trimming_Paired
use rule trimming_Single from cutadapt as RNAseq_trimming_Single

if aligner == 'hisat2':
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
        snakefile: "../modules/hisat2/TEtranscripts/hisat2.smk"
        config: hisat2_config
    logger.info(f"hisat2_config: {hisat2_config}")
    use rule hisat2_align from hisat2 as RNAseq_hisat2_align
    use rule hisat2_index from hisat2 as RNAseq_hisat2_index
elif aligner == 'star':
    star_config = {
            "indir": cutadapt_config["outdir"],
            "outdir":  f"{outdir}/star",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "star": config.get('Procedure',{}).get('star')
            },
            "genome": {
                "fasta": config.get('genome',{}).get('fasta'),
                "gtf": config.get('genome',{}).get('gtf'),
                "index_dir": config.get('genome',{}).get('star_index_dir')
            }
        }
    module star:
        snakefile: "../modules/star/TEtranscripts/star.smk"
        config: star_config
    logger.info(f"star_config: {star_config}")
    use rule star_align from star as RNAseq_star_align
    use rule star_index from star as RNAseq_star_index
else:
    raise ValueError(f"Unsupported aligner: {aligner}")


TEtranscripts_config = {
        "indir": f"{outdir}/star" if aligner == 'star' else f"{outdir}/hisat2",
        "outdir":  f"{outdir}/TEtranscripts",
        "logdir": logdir,
        "single_samples": single_samples,
        "paired_samples": paired_samples,
        "Procedure": {
            "TEcount": config.get('Procedure',{}).get('TEcount') or 'TEcount',
            "TElocal": config.get('Procedure',{}).get('TElocal') or 'TElocal'
        },
        "genome": {
            "gtf": config.get('genome',{}).get('gtf'),
            "TEind": config.get('genome',{}).get('TEind'),
            "TE_gtf": config.get('genome',{}).get('TE_gtf')
        }
    }
logger.info(f"TEtranscripts_config: {TEtranscripts_config}")
module TEtranscripts:
    snakefile: "../modules/TEtranscripts/TEtranscripts.smk"
    config: TEtranscripts_config

use rule * from TEtranscripts as RNAseq_*

DESeq2_config = {
        "indir": TEtranscripts_config["outdir"],
        "outdir":  f"{outdir}/DESeq2",
        "logdir": logdir,
        "control_samples": config.get("control_samples", []),
        "control_group_name": config.get("control_group_name", "control"),
        "treatment_samples": config.get("treatment_samples", []),
        "treatment_group_name": config.get("treatment_group_name", "treatment"),
        "genome": {
            "geneIDAnno": config.get('genome',{}).get('geneIDAnno')
        },
        "Procedure": {
            "DESeq2": config.get('Procedure',{}).get('DESeq2') or 'DESeq2'
        }
    }
module DESeq2:
    snakefile: "../modules/DESeq2/DESeq2.smk"
    config: DESeq2_config
logger.info(f"DESeq2_config: {DESeq2_config}")
use rule DESeq2_TEcount from DESeq2 as RNAseq_DESeq2_TEcount