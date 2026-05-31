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
UmiTools_extract_config = {
        "indir": indir,
        "outdir": f"{outdir}/umi_tools_extract",
        "logdir": logdir,
        "Procedure": {
            "umi_tools": config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
            "extract_method": config.get('Procedure',{}).get('extract_method') or 'string'
        },
        "Params": {
            "umi_tools": {
                "bc_pattern": config.get('Params',{}).get('umi_tools',{}).get('bc_pattern') or 'NNNXXXXNN',
                "bc_pattern2": config.get('Params',{}).get('umi_tools',{}).get('bc_pattern2') or 'NNNXXXXNN',
            }
        }
    }
module UmiTools_extract:
    snakefile: "../modules/UmiTools/extract/UmiTools.smk"
    config: UmiTools_extract_config
logger.info(f"UmiTools_extract_config: {UmiTools_extract_config}")
use rule UmiTools_extract_single from UmiTools_extract as CLIP_UmiTools_extract_single
use rule UmiTools_extract_paired from UmiTools_extract as CLIP_UmiTools_extract_paired
cutadapt_config = {
        "indir": UmiTools_extract_config["outdir"],
        "outdir":  f"{outdir}/cutadapt",
        "logdir": logdir,
        "mode": "UMI",
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": {
                "quality": config.get('Params',{}).get("trim_galore", {}).get('quality')
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
        snakefile: "../modules/hisat2/hisat2.smk"
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
                "STAR": config.get('Procedure',{}).get('STAR')
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
umi_tools_dedup_config = {
        "indir": star_config["outdir"] if aligner == 'star' else hisat2_config["outdir"],
        "outdir":  f"{outdir}/umi_tools_dedup",
        "logdir": logdir,
        "Procedure": {
            "umi_tools": config.get('Procedure',{}).get('umi_tools')
        },
        "Params": {
            "umi_tools": {
                "method": config.get('Params',{}).get('umi_tools',{}).get('method') or 'unique'
            }
        }
    }
module UmiTools_dedup:
    snakefile: "../modules/UmiTools/UmiTools.smk"
    config: umi_tools_dedup_config
logger.info(f"umi_tools_dedup_config: {umi_tools_dedup_config}")
use rule umi_tools_dedup_for_star from UmiTools_dedup as CLIP_umi_tools_dedup_for_star
use rule umi_tools_dedup_for_hisat2 from UmiTools_dedup as CLIP_umi_tools_dedup_for_hisat2

genome_config = {
        "genome": {
            "fasta": config.get('genome',{}).get('fasta')
        },
        "outdir": outdir,
        "logdir": logdir,
        "Procedure": {
            "samtools": config.get('Procedure',{}).get('samtools')
        }
    }
module genome:
    snakefile: "../modules/genome/genome.smk"
    config: genome_config
logger.info(f"genome_config: {genome_config}")
use rule chromosome_sizes from genome as CLIP_chromosome_sizes

bedtools_config = {
        "indir": umi_tools_dedup_config["outdir"],
        "outdir":  f"{outdir}/bedtools",
        "logdir": logdir,
        "Procedure": {
            "bedtools": config.get('Procedure',{}).get('bedtools')
        },
        "genome": {
            "chrom_sizes": config.get('genome',{}).get('chrom_sizes') or f"{outdir}/genome/chrom.sizes"
        }
    }
module bedtools:
    snakefile: "../modules/bedtools/bedtools.smk"
    config: bedtools_config
logger.info(f"bedtools_config: {bedtools_config}")
use rule iCLIP_bedtools from bedtools as CLIP_bedtools

PureCLIP_config = {
        "indir": umi_tools_dedup_config["outdir"],
        "outdir":  f"{outdir}/PureCLIP",
        "logdir": logdir,
        "Procedure": {
            "PureCLIP": config.get('Procedure',{}).get('PureCLIP')
        },
        "genome": {
            "fasta": config.get('genome',{}).get('fasta')
        }
    }
module PureCLIP:
    snakefile: "../modules/PureCLIP/PureCLIP.smk"
    config: PureCLIP_config
logger.info(f"PureCLIP_config: {PureCLIP_config}")
use rule pureclip from PureCLIP as CLIP_pureclip

track_config = {
        "indir": bedtools_config["outdir"],
        "outdir":  f"{outdir}/track",
        "logdir": logdir,
        "samples": single_samples + paired_samples,
        "igv": config.get('Params', {}).get('igv', {}),
    }
module track:
    snakefile: "../modules/track/track.smk"
    config: track_config
logger.info(f"track_config: {track_config}")
use rule igv_track_iclip from track as CLIP_igv_track_iclip
use rule ucsc_track_iclip from track as CLIP_ucsc_track_iclip
use rule ucsc_track_bedtools from track as CLIP_ucsc_track_bedtools
use rule igv_track_bedtools from track as CLIP_igv_track_bed



