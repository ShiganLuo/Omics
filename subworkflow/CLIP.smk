shell.prefix("set -x; set -e;")
from snakemake.logging import logger
indir = config.get("Params", {}).get("workflow", {}).get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow", {}).get("outdir", "output")
logdir = config.get("Params", {}).get("workflow", {}).get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR","")
outfiles = config.get("outfiles", [])
paired_samples = config.get("sample_params", {}).get("paired_samples", [])
single_samples = config.get("sample_params", {}).get("single_samples", [])
aligner = config.get("Params", {}).get("workflow", {}).get("aligner") or "star"
trimmer = config.get('Params', {}).get('workflow', {}).get('trimmer') or "trim_galore"

genome = config.get("genome", {}).get("default")
genome_ref = config.get("genome", {}).get("references", {}).get(genome, {}) if genome else {}

rule all:
    input:
        outfiles
fastqc_raw_config = {
        "resources": config.get("resources", {}),
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
        "resources": config.get("resources", {}),
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
            "umi_tools": config.get('Params',{}).get('umi_tools', {}),
            "UmiTools": {
                "keep_intermediates": config.get('Params',{}).get('UmiTools',{}).get('keep_intermediates', True)
            },
        }
    }
module UmiTools_extract:
    snakefile: "../modules/UmiTools/extract/UmiTools.smk"
    config: UmiTools_extract_config
logger.info(f"UmiTools_extract_config: {UmiTools_extract_config}")
use rule UmiTools_extract_single from UmiTools_extract as CLIP_UmiTools_extract_single
use rule UmiTools_extract_paired from UmiTools_extract as CLIP_UmiTools_extract_paired
trimmed_dir = f"{outdir}/cutadapt"
if trimmer == "cutadapt":
    cutadapt_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": UmiTools_extract_config["outdir"],
        "outdir": trimmed_dir,
        "logdir": logdir,
        "mode": "UMI",
        "Procedure": {
            "cutadapt": config.get('Procedure',{}).get('cutadapt')
        },
        "Params": {
            "cutadapt": config.get('Params',{}).get("cutadapt", {}),
        }
    }
    module cutadapt:
        snakefile: "../modules/cutadapt/cutadapt.smk"
        config: cutadapt_config
    logger.info(f"cutadapt_config: {cutadapt_config}")
    use rule trimming_Paired from cutadapt as CLIP_cutadapt_Paired
    use rule trimming_Single from cutadapt as CLIP_cutadapt_Single
elif trimmer == "trimmomatic":
    raise ValueError("trimmomatic 模块不支持 UMI mode,CLIP 的 Params.workflow.trimmer 只能是 'cutadapt' 或 'trim_galore'")
else:
    logger.info("use default trimmer: trim_galore")
    trim_galore_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": UmiTools_extract_config["outdir"],
        "outdir": trimmed_dir,
        "logdir": logdir,
        "mode": "UMI",
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": config.get('Params',{}).get("trim_galore", {}),
        }
    }
    module trim_galore:
        snakefile: "../modules/trim-galore/trim-galore.smk"
        config: trim_galore_config
    logger.info(f"trim_galore_config: {trim_galore_config}")
    use rule trimming_Paired from trim_galore as CLIP_trim_galore_Paired
    use rule trimming_Single from trim_galore as CLIP_trim_galore_Single

fastqc_trimmed_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": trimmed_dir,
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
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "indir": trimmed_dir,
            "outdir":  f"{outdir}/hisat2",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "hisat2": config.get('Procedure',{}).get('hisat2')
            },
            "Params": {
                "hisat2": config.get('Params',{}).get('hisat2', {})
            },
            "genome": {
                "fasta": genome_ref.get("fasta"),
                "index_dir": genome_ref.get("hisat2_index_dir")
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
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "indir": trimmed_dir,
            "outdir":  f"{outdir}/star",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Procedure": {
                "STAR": config.get('Procedure',{}).get('STAR')
            },
            "Params": {
                "star": config.get('Params',{}).get('star', {})
            },
            "genome": {
                "fasta": genome_ref.get("fasta"),
                "gtf": genome_ref.get("gtf"),
                "index_dir": genome_ref.get("star_index_dir")
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
        "resources": config.get("resources", {}),
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
            "umi_tools": config.get('Params',{}).get('umi_tools', {}),
            "UmiTools": {
                "keep_intermediates": config.get('Params',{}).get('UmiTools',{}).get('keep_intermediates', True)
            },
        }
    }
module UmiTools_dedup:
    snakefile: "../modules/UmiTools/UmiTools.smk"
    config: umi_tools_dedup_config
logger.info(f"umi_tools_dedup_config: {umi_tools_dedup_config}")
use rule umi_tools_dedup_for_star from UmiTools_dedup as CLIP_umi_tools_dedup_for_star
use rule umi_tools_dedup_for_hisat2 from UmiTools_dedup as CLIP_umi_tools_dedup_for_hisat2

_sexcheck_bam_dir = star_config["outdir"] if aligner == "star" else hisat2_config["outdir"]
sexcheck_config = {
    "resources": config.get("resources", {}),
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
            .get("sex_chroms"),
        "gtf": (config.get("genome", {}) or {}).get("references", {})
            .get(config.get("genome", {}).get("default"), {})
            .get("gtf")
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

qc_summary_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "qc_dirs": {
        "fastqc_raw": fastqc_raw_config["outdir"],
        "fastqc_trimmed": fastqc_trimmed_config["outdir"],
        "sexcheck": sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "sexcheck": "sexcheck.tsv",
    },
    "qc_anchor_samples": {},
}
if trimmer == "cutadapt":
    qc_summary_config["qc_dirs"]["cutadapt"] = trimmed_dir
    if paired_samples:
        qc_summary_config["qc_anchors"]["cutadapt"] = "{sample}/{sample}.cutadapt_report.txt"
        qc_summary_config["qc_anchor_samples"]["cutadapt"] = paired_samples
    if single_samples:
        qc_summary_config["qc_dirs"]["cutadapt_single"] = trimmed_dir
        qc_summary_config["qc_anchors"]["cutadapt_single"] = "{sample}/{sample}.single.cutadapt_report.txt"
        qc_summary_config["qc_anchor_samples"]["cutadapt_single"] = single_samples
else:
    qc_summary_config["qc_dirs"]["trim_galore"] = trimmed_dir
    if paired_samples:
        qc_summary_config["qc_anchors"]["trim_galore"] = "{sample}/trimming_statistics_1.txt"
        qc_summary_config["qc_anchor_samples"]["trim_galore"] = paired_samples
    if single_samples:
        qc_summary_config["qc_dirs"]["trim_galore_single"] = trimmed_dir
        qc_summary_config["qc_anchors"]["trim_galore_single"] = "{sample}/trimming_statistics.txt"
        qc_summary_config["qc_anchor_samples"]["trim_galore_single"] = single_samples
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
            "fasta": genome_ref.get("fasta"),
            "gtf": genome_ref.get("gtf")
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
            "chrom_sizes": genome_ref.get("chrom_sizes") or f"{outdir}/genome/chrom.sizes"
        }
    }
module bedtools:
    snakefile: "../modules/bedtools/bedtools.smk"
    config: bedtools_config
logger.info(f"bedtools_config: {bedtools_config}")
use rule iCLIP_bedtools from bedtools as CLIP_bedtools

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
        "gtf": genome_ref.get("gtf"),
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
        "resources": config.get("resources", {}),
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
            "fasta": genome_ref.get("fasta")
        }
    }
module PureCLIP:
    snakefile: "../modules/PureCLIP/PureCLIP.smk"
    config: PureCLIP_config
logger.info(f"PureCLIP_config: {PureCLIP_config}")
use rule pureclip from PureCLIP as CLIP_pureclip

track_config = {
        "resources": config.get("resources", {}),
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

