shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
samples = config.get("samples", [])
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})

rule all:
    input:
        outfiles

# ==============================================================================
# Step 1: Raw FastQC
# ==============================================================================
fastqc_raw_config = {
        "indir": indir,
        "outdir":  f"{outdir}/QC/1_raw_fastqc",
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
use rule fastqc from fastqc_raw as PeakCalling_fastqc_raw

# ==============================================================================
# Step 2: Trim Galore
# ==============================================================================
trim_galore_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": indir,
        "outdir": f"{outdir}/common/2_trimmed_fastq",
        "logdir": logdir,
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": {
                "quality": config.get('Params',{}).get("trim_galore", {}).get('quality')
            }
        },
    }
module trim_galore:
    snakefile: "../modules/trim-galore/trim-galore.smk"
    config: trim_galore_config
logger.info(f"TrimGalore parameters: {trim_galore_config}")
use rule trimming_Paired from trim_galore as PeakCalling_trimming_Paired
use rule trimming_Single from trim_galore as PeakCalling_trimming_Single

# ==============================================================================
# Step 3: Trimmed FastQC
# ==============================================================================
fastqc_trimmed_config = {
        "indir": trim_galore_config["outdir"],
        "outdir":  f"{outdir}/QC/2_trimmed_fastqc",
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
use rule fastqc from fastqc_trimmed as PeakCalling_fastqc_trimmed

# ==============================================================================
# Step 4: Bowtie2 Alignment
# ==============================================================================
bowtie2_config = {
    "indir": trim_galore_config["outdir"],
    "outdir": f"{outdir}/common/3_raw_bam",
    "logdir": logdir,
    "Procedure": {
        "bowtie2-build": config.get("Procedure", {}).get("bowtie2-build"),
        "bowtie2": config.get("Procedure", {}).get("bowtie2")
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta"),
        "index_prefix": config.get("genome", {}).get("bowtie2_index_prefix")
    }
}
module bowtie2:
    snakefile: "../modules/bowtie2/bowtie2.smk"
    config: bowtie2_config
logger.info(f"bowtie2_config: {bowtie2_config}")
use rule bowtie2_index from bowtie2 as PeakCalling_bowtie2_index
use rule bowtie2_align_paired from bowtie2 as PeakCalling_bowtie2_align_paired
use rule bowtie2_align_single from bowtie2 as PeakCalling_bowtie2_align_single

# =============================================================================
# Step 5: Add Read Groups + Mark Duplicates (GATK4)
# Consistent with nf-core/chipseq: AddOrReplaceReadGroups + MarkDuplicates.
# Duplicates are flagged but not removed; downstream tools (MACS3, FRiP)
# handle duplicate filtering via --keep-dup.
# =============================================================================
gatk_prepare_config = {
    "indir": bowtie2_config["outdir"],
    "outdir": f"{outdir}/common/4_markdup_bam",
    "logdir": logdir,
    "input_bam_substring": "",
    "Procedure": {
        "gatk": config.get("Procedure", {}).get("gatk") or "gatk",
        "samtools": config.get("Procedure", {}).get("samtools") or "samtools"
    },
    "Params": {
        "gatk": config.get("Params", {}).get("gatk", {})
    },
    "addReadsGroup": config.get("addReadsGroup", {}),
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module gatk_prepare:
    snakefile: "../modules/gatk/gatk_prepare.smk"
    config: gatk_prepare_config
logger.info(f"gatk_prepare_config: {gatk_prepare_config}")
use rule addReadsGroup from gatk_prepare as PeakCalling_addReadsGroup
use rule MarkDuplicates from gatk_prepare as PeakCalling_MarkDuplicates

# ==============================================================================
# Step 6: BigWig Track Generation (bamCoverage)
# Generates normalized coverage tracks for visualization in genome browsers.
# Uses the existing igv module's dedup + wig rules.
# ==============================================================================
igv_config = {
    "indir": bowtie2_config["outdir"],
    "outdir": f"{outdir}/tracks",
    "logdir": logdir,
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools") or "samtools",
        "bamCoverage": config.get("Procedure", {}).get("bamCoverage") or "bamCoverage"
    },
    "Params": {
        "bamCoverage": {
            "binSize": config.get("Params", {}).get("bamCoverage", {}).get("binSize") or 50,
            "normalizeUsing": config.get("Params", {}).get("bamCoverage", {}).get("normalizeUsing") or "CPM",
            "offset": config.get("Params", {}).get("bamCoverage", {}).get("offset"),
            "extendReads": config.get("Params", {}).get("bamCoverage", {}).get("extendReads") or False
        }
    }
}
module igv:
    snakefile: "../modules/igv/igv.smk"
    config: igv_config
logger.info(f"igv_config: {igv_config}")
use rule samtools_dedup from igv as PeakCalling_dedup
use rule wig from igv as PeakCalling_bigwig

track_config = {
        "indir": igv_config["outdir"],
        "outdir":  igv_config["outdir"],
        "logdir": logdir,
        "samples": single_samples + paired_samples,
        "igv": config.get('Params', {}).get('igv', {}),
    }

module track:
    snakefile: "../modules/track/track.smk"
    config: track_config
logger.info(f"track_config: {track_config}")
use rule * from track as PeakCalling_*
# =============================================================================
# Step 7: MACS3 Peak Calling
# =============================================================================
macs3_config = {
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/peaks",
    "logdir": logdir,
    "samples": ip_samples,
    "ip_samples": ip_samples,
    "input_samples": input_samples,
    "sample_ip_input_map": sample_ip_input_map,
    "Procedure": {
        "macs3": config.get("Procedure", {}).get("macs3")
    },
    "Params": {
        "macs3": config.get("Params", {}).get("macs3", {})
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module macs3:
    snakefile: "../modules/macs3/macs3.smk"
    config: macs3_config
logger.info(f"macs3_config: {macs3_config}")
use rule macs3_callpeak from macs3 as PeakCalling_macs3_callpeak

# =============================================================================
# Step 8: FRIP Score (Fraction of Reads in Peaks)
# Key ChIP-seq QC metric: measures enrichment quality.
# FRiP = reads_in_peaks / total_mapped_reads (target >= 0.2 for good data)
# =============================================================================
frip_score_config = {
    "indir": gatk_prepare_config["outdir"],
    "outdir": f"{outdir}/QC/3_frip_score",
    "logdir": logdir,
    "peaks_indir": macs3_config["outdir"],
    "samples": ip_samples,
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools") or "samtools",
        "bedtools": config.get("Procedure", {}).get("bedtools") or "bedtools"
    }
}
module frip_score:
    snakefile: "../modules/frip_score/frip_score.smk"
    config: frip_score_config
logger.info(f"frip_score_config: {frip_score_config}")
use rule frip_score from frip_score as PeakCalling_frip_score

# =============================================================================
# Step 9: HOMER Peak Annotation
# Annotates peaks with genomic features (promoter, intron, intergenic, etc.)
# and nearest gene information.
# ==============================================================================
homer_config = {
    "indir": macs3_config["outdir"],
    "outdir": f"{outdir}/annotation",
    "logdir": logdir,
    "samples": ip_samples,
    "Procedure": {
        "annotatePeaks": config.get("Procedure", {}).get("annotatePeaks") or "annotatePeaks.pl"
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta"),
        "gtf": config.get("genome", {}).get("gtf")
    }
}
module homer:
    snakefile: "../modules/homer/homer.smk"
    config: homer_config
logger.info(f"homer_config: {homer_config}")
use rule homer_annotatepeaks from homer as PeakCalling_homer_annotatepeaks
