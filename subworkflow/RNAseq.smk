shell.prefix("set -x; set -e;")
from snakemake.logging import logger
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
aligner = config.get('aligner', 'hisat2')
trimmer = config.get('trimmer', 'cutadapt')
outfiles = config.get("outfiles", [])
rule all:
    input:
        outfiles
if trimmer == "cutadapt":
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
    use rule trimming_Paired from cutadapt as RNAseq_trimming_Paireds
    use rule trimming_Single from cutadapt as RNAseq_trimming_Single
elif trimmer == "trimmomatic":
    trimmomatic_config = {
            "indir": indir,
            "outdir":  f"{outdir}/trimmomatic",
            "logdir": logdir,
            "Procedure": {
                "trimmomatic": config.get('Procedure',{}).get('trimmomatic')
            },
            "Params": {
                "trimmomatic": {
                    "adapter_pe": config.get('Params',{}).get("trimmomatic", {}).get('adapter_pe'),
                    "adapter_se": config.get('Params',{}).get("trimmomatic", {}).get('adapter_se')
                }
            }
        }
    module trimmomatic:
        snakefile: "../modules/trimmomatic/trimmomatic.smk"
        config: trimmomatic_config
    logger.info(f"trimmomatic_config: {trimmomatic_config}")
    use rule trimmomatic_Paired from trimmomatic as RNAseq_trimmomatic_Paireds
    use rule trimmomatic_Single from trimmomatic as RNAseq_trimmomatic_Singles
else:
    raise ValueError(f"Unsupported trimmer: {trimmer}")

if aligner == 'hisat2':
    hisat2_config_for_TEtranscripts = {
            "indir": cutadapt_config["outdir"] if trimmer == "cutadapt" else trimmomatic_config["outdir"],
            "outdir":  f"{outdir}/TEtranscripts/bam",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "hisat2": config.get('Procedure',{}).get('hisat2')
            },
            "Params": {
                "hisat2": {
                    "score_min": config.get('Params',{}).get('hisat2', {}).get('score_min') or "L,0,-0.2",
                    "flag_params": config.get('Params',{}).get('hisat2', {}).get('flag_params') or "--no-mixed --no-discordant",
                    "k": config.get('Params',{}).get('hisat2', {}).get('k') or 100
                }
            },
            "genome": {
                "fasta": config.get('genome',{}).get('fasta'),
                "index_prefix": config.get('genome',{}).get('hisat2_index_prefix')
            }
        }
    module hisat2_for_TEtranscripts:
        snakefile: "../modules/hisat2/hisat2.smk"
        config: hisat2_config_for_TEtranscripts
    logger.info(f"hisat2_config: {hisat2_config_for_TEtranscripts}")
    use rule hisat2_align from hisat2 as RNAseq_hisat2_align_for_TEtranscripts
    use rule hisat2_index from hisat2 as RNAseq_hisat2_index_for_TEtranscripts
elif aligner == 'star':
    star_config_for_TEtranscripts = {
            "indir": cutadapt_config["outdir"] if trimmer == "cutadapt" else trimmomatic_config["outdir"],
            "outdir":  f"{outdir}/TEtranscripts/bam",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "star": config.get('Procedure',{}).get('star')
            },
            "Params": {
                "STAR": {
                    "alignEndsType": config.get('Params',{}).get('STAR', {}).get('alignEndsType') or "Local",
                    "outFilterMismatchNoverReadLmax": config.get('Params',{}).get('STAR', {}).get('outFilterMismatchNoverReadLmax') or 1.0,
                    "outFilterMismatchNmax": config.get('Params',{}).get('STAR', {}).get('outFilterMismatchNmax') or 10,
                    "outFilterMultimapNmax": config.get('Params',{}).get('STAR', {}).get('outFilterMultimapNmax') or 100,
                    "winAnchorMultimapNmax": config.get('Params',{}).get('STAR', {}).get('winAnchorMultimapNmax') or 100
                }
            },
            "genome": {
                "fasta": config.get('genome',{}).get('fasta'),
                "gtf": config.get('genome',{}).get('gtf'),
                "index_dir": config.get('genome',{}).get('star_index_dir')
            }
        }
    module star_for_TEtranscripts:
        snakefile: "../modules/star/star.smk"
        config: star_config_for_TEtranscripts
    logger.info(f"star_config: {star_config_for_TEtranscripts}")
    use rule star_align from star as RNAseq_star_align_for_TEtranscripts
    use rule star_index from star as RNAseq_star_index_for_TEtranscripts
else:
    raise ValueError(f"Unsupported aligner: {aligner}")


TEtranscripts_config = {
        "indir": f"{outdir}/star" if aligner == 'star' else f"{outdir}/hisat2",
        "outdir":  f"{outdir}/TEtranscripts",
        "logdir": logdir,
        "samples": single_samples + paired_samples,
        "ROOT_DIR": ROOT_DIR,
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
        "ROOT_DIR": ROOT_DIR,
        "control_samples": config.get("control_samples", []),
        "control_group_name": config.get("control_group_name", "control"),
        "treatment_samples": config.get("treatment_samples", []),
        "experimental_group_name": config.get("experimental_group_name", "treatment"),
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



StringTie_config = {
        "indir": f"{outdir}/star" if aligner == 'star' else f"{outdir}/hisat2",
        "outdir":  f"{outdir}/stringTie",
        "logdir": logdir,
        "samples": single_samples + paired_samples,
        "genome": {
            "gtf": config.get('genome',{}).get('gtf')
        },
        "Procedure": {
            "stringtie": config.get("Procedure", {}).get("stringtie") or "stringtie"
        }
    }
logger.info(f"StringTie_config: {StringTie_config}")
module StringTie:
    snakefile: "../modules/StringTie/StringTie.smk"
    config: StringTie_config
use rule * from StringTie as RNAseq_*