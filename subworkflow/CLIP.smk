shell.prefix("set -x; set -e;")
from snakemake.logging import logger
indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
ROOT_DIR = config.get("ROOT_DIR","")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
aligner = config.get('aligner', 'star')

rule all:
    input:
        outfiles
fastqc_raw_config = {
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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
            "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
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
            "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
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
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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

# ============================================
# Sex check (coverage on sex chromosomes)
# ============================================
# Wired unconditionally over the alignment BAMs; node.py registers the
# sexcheck outputs (subworkflow.md QC wiring).
# aligner-conditional BAM source (config dicts are both defined above)
_sexcheck_bam_dir = star_config["outdir"] if aligner == "star" else hisat2_config["outdir"]
sexcheck_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/sexcheck",
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "sample_bams": {
        s: f"{_sexcheck_bam_dir}/{s}/{s}.bam"
        for s in paired_samples + single_samples
    },
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools")
    },
    "genome": {
        "sex_chroms": (config.get("genome", {}) or {}).get("references", {})
            .get(config.get("genome", {}).get("default"), {})
            .get("sex_chroms")
    },
    "Params": {
        "sexcheck": config.get("Params", {}).get("sexcheck", {})
    },
}
module sexcheck:
    snakefile: "../modules/sexcheck/sexcheck.smk"
    config: sexcheck_config
logger.info(f"sexcheck_config: {sexcheck_config}")
use rule sexcheck_idxstats from sexcheck as CLIP_sexcheck_idxstats
use rule sexcheck from sexcheck as CLIP_sexcheck


# ============================================
# QC Summary (cross-sample metric aggregation)
# ============================================
# Wired unconditionally; node.py registers the qc_summary outputs (see
# subworkflow.md QC wiring spec). Only stages with qc_summary parsers are
# declared; unregistered stages stay best-effort scans.
qc_summary_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "qc_dirs": {
        "fastqc_raw": fastqc_raw_config["outdir"],
        "fastqc_trimmed": fastqc_trimmed_config["outdir"],
        "cutadapt": cutadapt_config["outdir"],
        "sexcheck": sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "cutadapt": "{sample}/{sample}.cutadapt_report.txt",
        "sexcheck": "sexcheck.tsv",
    },
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as CLIP_qc_summary



genome_config = {
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "genome": {
            "fasta": config.get('genome',{}).get('fasta'),
            "gtf": config.get('genome',{}).get('gtf')
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
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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


# ============================================================
# Region-focused multi-sample track stacking plots (modules/regionplot)
#   Consumes this workflow's bigwigs; gated by Params.regionplot in
#   node.py (subworkflow wiring is unconditional per repo convention).
# ============================================================
regionplot_params = config.get("Params", {}).get("regionplot", {})
regionplot_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/regionplot",
    "logdir": f"{logdir}/group",
    "samples": single_samples + paired_samples,
    "tracks": [{"label": "plus", "dir": bedtools_config["outdir"], "suffix": ".plus.bw"},
              {"label": "minus", "dir": bedtools_config["outdir"], "suffix": ".minus.bw"}],
    "Params": {
        "regionplot": regionplot_params
    },
    "genome": {
        "gtf": config.get("genome", {}).get("gtf"),
    },
    "resources": config.get("resources", {}),
}
module regionplot:
    snakefile: "../modules/regionplot/regionplot.smk"
    config: regionplot_config
logger.debug(f"regionplot_config: {regionplot_config}")
use rule regionplot from regionplot as CLIP_regionplot
use rule regionplot_result from regionplot as CLIP_regionplot_result

PureCLIP_config = {
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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
        "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
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



