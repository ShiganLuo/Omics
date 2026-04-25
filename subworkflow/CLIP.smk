shell.prefix("set -x; set -e;")
from snakemake.logging import logger
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
aligner = config.get('aligner', 'star')

rule all:
    input:
        outfiles
fastqc_raw_config = {
        "indir": indir,
        "outdir":  f"{outdir}/fastqc/raw",
        "logdir": logdir,
        "log_suffix": "raw.txt",
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "fastqc": config.get("Procedure", {}).get("fastqc") or "fastqc"
        }
    }
module fastqc_raw:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_raw_config
logger.info(f"fastqc_raw_config: {fastqc_raw_config}")
use rule fastqc from fastqc_raw as CLIP_fastqc

cutadapt_config = {
        "indir": indir,
        "outdir":  f"{outdir}/cutadapt",
        "logdir": logdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": {
                "quality": config.get('Params',{}).get("trim_galore", {}).get('quality') or 25
            }
        }
    }
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
logger.info(f"cutadapt_config: {cutadapt_config}")
use rule trimming_Paired from cutadapt as CLIP_trimming_Paired
use rule trimming_Single from cutadapt as CLIP_trimming_Single

fastqc_trimmed_config = {
        "indir": cutadapt_config["outdir"],
        "outdir":  f"{outdir}/fastqc/trimmed",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "log_suffix": "trimmed.txt",
        "Procedure": {
            "fastqc": config.get("Procedure", {}).get("fastqc")
        }
    }
module fastqc_trimmed:
    snakefile: "../modules/fastqc/fastqc.smk"
    config: fastqc_trimmed_config
logger.info(f"fastqc_trimmed_config: {fastqc_trimmed_config}")
use rule fastqc from fastqc_trimmed as CLIP_fastqc_trimmed


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
                "fasta": config.get('genome',{}).get('fasta'),
                "index_dir": config.get('genome',{}).get('hisat2_index_dir')
            }
        }
    module hisat2:
        snakefile: "../modules/hisat2/TEtranscripts/hisat2.smk"
        config: hisat2_config
    logger.info(f"hisat2_config: {hisat2_config}")
    use rule hisat2_align from hisat2 as CLIP_hisat2_align
    use rule hisat2_index from hisat2 as CLIP_hisat2_index
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
        snakefile: "../modules/star/star.smk"
        config: star_config
    logger.info(f"star_config: {star_config}")
    use rule star_align from star as CLIP_star_align
    use rule star_index from star as CLIP_star_index
else:
    raise ValueError(f"Unsupported aligner: {aligner}")

# igv_config = {
#         "indir": hisat2_config["outdir"],
#         "outdir":  f"{outdir}/igv",
#         "logdir": logdir,
#         "Procedure": {
#             "samtools": config.get('Procedure',{}).get('samtools'),
#             "deepTools": config.get('Procedure',{}).get('deepTools')
#         }
#     }
# module igv:
#     snakefile: "../modules/igv/igv.smk"
#     config: igv_config
# logger.info(f"igv_config: {igv_config}")
# use rule dedup from igv as CLIP_dedup
# use rule wig from igv as CLIP_wig

