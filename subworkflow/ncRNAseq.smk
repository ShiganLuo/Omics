shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
all_samples = config.get("samples", [])
outfiles = config.get("outfiles", [])
aligner = config.get("Procedure", {}).get("aligner") or "star"

rule all:
    input:
        outfiles
fastqc_raw_config = {
        "ROOT_DIR": ROOT_DIR,
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

# ── 0. Demultiplex: 3' adapter removal + PCR duplicate removal ───────────────
demultiplexer_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": indir,
        "outdir": f"{outdir}/common/2_trimmed_dedup_fastq/jla-demultiplexer",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Params": {
            "demultiplexer": config.get("Params", {}).get("demultiplexer", {})
        },
    }
module demultiplexer:
    snakefile: "../modules/demultiplexer/demultiplexer.smk"
    config: demultiplexer_config
logger.info(f"demultiplexer_config: {demultiplexer_config}")
use rule demultiplex_trim_dedup from demultiplexer as ncRNAseq_demultiplex_trim_dedup

trim_galore_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": demultiplexer_config["outdir"],
        "outdir": f"{outdir}/common/2_trimmed_dedup_fastq",
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
use rule trimming_Paired from trim_galore as ncRNAseq_trimming_Paired
use rule trimming_Single from trim_galore as ncRNAseq_trimming_Single

# ── 0.5 Subsample: seqtk subsample for abundant small RNAs ───────────────────
subsample_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": trim_galore_config["outdir"],
        "outdir": f"{outdir}/common/trimmed_subsampled_fastq",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Params": {
            "subsample": {
                "abundant_rnas": config.get("Params", {}).get("ncRNAseq", {}).get("abund_small_rnas", []),
                "n": config.get("Params", {}).get("ncRNAseq", {}).get("subsample_n", 100000),
                "seed": config.get("Params", {}).get("ncRNAseq", {}).get("subsample_seed", 42),
            }
        },
    }
module subsample:
    snakefile: "../modules/subsample/subsample.smk"
    config: subsample_config
logger.info(f"subsample_config: {subsample_config}")
use rule subsample_fastq from subsample as ncRNAseq_subsample_fastq

fastqc_trimmed_config = {
        "ROOT_DIR": ROOT_DIR,
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

# ── 2. Align ─────────────────────────────────────────────────────────────────
STAR = config.get("Procedure", {}).get("STAR") or "STAR"
SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
genome_fasta = config.get("genome", {}).get("fasta")
star_index_dir = config.get("genome", {}).get("star_index_dir")
smallrna_fasta = config.get("genome", {}).get("smallrna_fasta")
smallrna_bed = config.get("genome", {}).get("smallrna_bed")
smallrna_star_index = config.get("genome", {}).get("smallrna_star_index")

if aligner == "hisat2":
    hisat2_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": subsample_config["outdir"],
        "outdir": f"{outdir}/common/3_raw_bam",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "hisat2": config.get("Procedure", {}).get("hisat2"),
            "hisat2-build": config.get("Procedure", {}).get("hisat2-build")
        },
        "genome": {
            "fasta": genome_fasta,
            "hisat2_index_prefix": config.get("genome", {}).get("hisat2_index_prefix")
        }
    }
    logger.info(f"hisat2_config: {hisat2_config}")
    module hisat2:
        snakefile: "../modules/hisat2/ncRNAseq/hisat2.smk"
        config: hisat2_config
    use rule hisat2_align_ncRNAseq_single from hisat2 as ncRNAseq_hisat2_align

elif aligner == "star":
    # ── STAR index for genome (auto-build when star_index_dir is null) ──
    if not star_index_dir:
        star_genome_idx_config = {
            "ROOT_DIR": ROOT_DIR,
            "outdir": f"{outdir}/common/3_raw_bam",
            "logdir": logdir,
            "Procedure": {"STAR": STAR},
            "Params": {"STAR": {}},
            "genome": {
                "fasta": genome_fasta,
                "gtf": config.get("genome", {}).get("gtf"),
            }
        }
        logger.info(f"star_genome_idx_config: {star_genome_idx_config}")

        module star_genome_idx:
            snakefile: "../modules/star/star.smk"
            config: star_genome_idx_config

        use rule star_index from star_genome_idx as ncRNAseq_star_index_genome

        star_index_dir = f"{outdir}/common/3_raw_bam/index"

    star_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": subsample_config["outdir"],
        "outdir": f"{outdir}/common/3_raw_bam",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "STAR": STAR
        },
        "Params": {
            "STAR": {
                "genomeLoad": config.get("Params", {}).get("STAR", {}).get("genomeLoad") or "LoadAndRemove",
                "limitBAMsortRAM": config.get("Params", {}).get("STAR", {}).get("limitBAMsortRAM") or 20000000000,
                "outReadsUnmapped": config.get("Params", {}).get("STAR", {}).get("outReadsUnmapped") or "Fastx",
                "outFilterMultimapNmax": config.get("Params", {}).get("STAR", {}).get("outFilterMultimapNmax") or 99999,
                "outFilterMismatchNoverLmax": config.get("Params", {}).get("STAR", {}).get("outFilterMismatchNoverLmax") or 0.1,
                "outFilterMatchNminOverLread": config.get("Params", {}).get("STAR", {}).get("outFilterMatchNminOverLread") or 0.66,
                "alignSJoverhangMin": config.get("Params", {}).get("STAR", {}).get("alignSJoverhangMin") or 999,
                "alignSJDBoverhangMin": config.get("Params", {}).get("STAR", {}).get("alignSJDBoverhangMin") or 999
            }
        },
        "genome": {
            "fasta": genome_fasta,
            "gtf": config.get("genome", {}).get("gtf"),
            "index_dir": star_index_dir
        }
    }
    logger.info(f"star_config: {star_config}")
    module star:
        snakefile: "../modules/star/star.smk"
        config: star_config
    use rule star_align from star as ncRNAseq_star_align

elif aligner == "star_3pass":
    # ================================================================
    # Three-pass STAR alignment for canonical small RNA quantification
    #
    # Pass 1:  genome alignment (relaxed, multimapping allowed)
    #          → extract reads overlapping small RNA genes
    # Pass 2:  align extracted reads to canonical small RNA FASTA
    #          (EndToEnd, clipped, strict read mismatch)
    # Pass 3a: re-align canonically-mapped reads to genome (strict)
    #          → extract those still overlapping small RNA genes
    # Pass 3b: re-align unmapped reads from pass 2 to genome (strict)
    # Merge:   combine pass3a canonical + pass3b reads
    # ================================================================

    # ── Resolve smallRNA derived paths BEFORE any config dict uses them ──
    smallrna_bed = f"{outdir}/genome/smallrna/smallrna_genes.bed"
    smallrna_fasta = f"{outdir}/genome/smallrna/smallrna_genes_flank.fa"
    smallrna_star_index = f"{outdir}/genome/smallrna/index"

    # ── STAR index for genome (auto-build when star_index_dir is null) ──
    if not star_index_dir:
        star_genome_idx_config = {
            "ROOT_DIR": ROOT_DIR,
            "outdir": f"{outdir}/genome",
            "logdir": logdir,
            "Procedure": {"STAR": STAR},
            "Params": {"STAR": {}},
            "genome": {
                "fasta": genome_fasta,
                "gtf": config.get("genome", {}).get("gtf"),
            }
        }
        logger.info(f"star_genome_idx_config: {star_genome_idx_config}")

        module star_genome_idx:
            snakefile: "../modules/star/star.smk"
            config: star_genome_idx_config

        use rule star_index from star_genome_idx as ncRNAseq_star_index_genome

        star_index_dir = f"{outdir}/genome/index"

    # ── Read three-pass params from config (with defaults) ──────────────
    p3p = config.get("Params", {}).get("star_3pass", {})
    pass1_params = p3p.get("pass1", {})
    pass2_params = p3p.get("pass2", {})
    pass3_params = p3p.get("pass3", {})

    # ── Pass 1 config: relaxed genome alignment ─────────────────────────
    star_pass1_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": subsample_config["outdir"],
        "outdir": f"{outdir}/common/3_raw_bam/pass1",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {"STAR": STAR},
        "Params": {"STAR": {
            "outFilterMultimapNmax": pass1_params.get("outFilterMultimapNmax", 1000),
            "alignIntronMin": pass1_params.get("alignIntronMin", 9999999),
            "outFilterMultimapScoreRange": pass1_params.get("outFilterMultimapScoreRange", 0),
            "outFilterMismatchNoverLmax": pass1_params.get("outFilterMismatchNoverLmax", 0.2),
        }},
        "genome": {"fasta": genome_fasta, "index_dir": star_index_dir}
    }
    logger.info(f"star_pass1_config: {star_pass1_config}")

    # ── Pass 2 config: canonical small RNA FASTA alignment ──────────────
    star_pass2_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": f"{outdir}/common/3_raw_bam/pass1_extract",
        "outdir": f"{outdir}/common/3_raw_bam/pass2",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {"STAR": STAR},
        "Params": {"STAR": {
            "outFilterMultimapNmax": pass2_params.get("outFilterMultimapNmax", 1000),
            "outFilterMultimapScoreRange": pass2_params.get("outFilterMultimapScoreRange", 0),
            "outFilterMismatchNoverLmax": pass2_params.get("outFilterMismatchNoverLmax", 0.2),
            "outFilterMismatchNoverReadLmax": pass2_params.get("outFilterMismatchNoverReadLmax", 0.05),
            "clip5pNbases": pass2_params.get("clip5pNbases", "20 0"),
            "clip3pNbases": pass2_params.get("clip3pNbases", "0 20"),
            "alignIntronMin": pass2_params.get("alignIntronMin", 9999999),
            "alignMatesGapMax": pass2_params.get("alignMatesGapMax", 500),
            "alignEndsType": pass2_params.get("alignEndsType", "EndToEnd"),
            "outReadsUnmapped": pass2_params.get("outReadsUnmapped", "Fastx"),
        }},
        "genome": {"fasta": smallrna_fasta, "index_dir": smallrna_star_index}
    }
    logger.info(f"star_pass2_config: {star_pass2_config}")

    # ── Pass 3a config: strict genome re-alignment (canonical reads) ────
    star_pass3a_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": f"{outdir}/common/3_raw_bam/pass2_mapped_fq",
        "outdir": f"{outdir}/common/3_raw_bam/pass3a",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {"STAR": STAR},
        "Params": {"STAR": {
            "outFilterMultimapNmax": pass3_params.get("outFilterMultimapNmax", 1000),
            "outFilterMultimapScoreRange": pass3_params.get("outFilterMultimapScoreRange", 0),
            "outFilterMismatchNoverLmax": pass3_params.get("outFilterMismatchNoverLmax", 0.025),
            "alignIntronMin": pass3_params.get("alignIntronMin", 9999999),
            "alignMatesGapMax": pass3_params.get("alignMatesGapMax", 500),
            "alignEndsType": pass3_params.get("alignEndsType", "Local"),
        }},
        "genome": {"fasta": genome_fasta, "index_dir": star_index_dir}
    }
    logger.info(f"star_pass3a_config: {star_pass3a_config}")

    # ── Pass 3b config: strict genome alignment (unmapped from pass 2) ──
    star_pass3b_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": f"{outdir}/common/3_raw_bam/pass2_unmapped_fq",
        "outdir": f"{outdir}/common/3_raw_bam/pass3b",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {"STAR": STAR},
        "Params": {"STAR": {
            "outFilterMultimapNmax": pass3_params.get("outFilterMultimapNmax", 1000),
            "outFilterMultimapScoreRange": pass3_params.get("outFilterMultimapScoreRange", 0),
            "outFilterMismatchNoverLmax": pass3_params.get("outFilterMismatchNoverLmax", 0.025),
            "alignIntronMin": pass3_params.get("alignIntronMin", 9999999),
            "alignMatesGapMax": pass3_params.get("alignMatesGapMax", 500),
            "alignEndsType": pass3_params.get("alignEndsType", "Local"),
        }},
        "genome": {"fasta": genome_fasta, "index_dir": star_index_dir}
    }
    logger.info(f"star_pass3b_config: {star_pass3b_config}")

    # ── Import genome module (extract smallRNA BED/FASTA) ──────────────
    genome_sm_config = {
        "ROOT_DIR": ROOT_DIR,
        "outdir": outdir,
        "logdir": logdir,
        "Procedure": {
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "Params": {
            "smallrna_types": config.get("Params", {}).get("ncRNAseq", {}).get(
                "smallrna_types", ["miRNA", "snRNA", "snoRNA", "rRNA", "misc_RNA", "scRNA", "scaRNA", "vaultRNA"]),
            "smallrna_flank": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_flank", 50)
        },
        "genome": {
            "fasta": genome_fasta,
            "gtf": config.get("genome", {}).get("gtf")
        }
    }
    logger.info(f"genome_sm_config: {genome_sm_config}")

    module genome_sm:
        snakefile: "../modules/genome/genome.smk"
        config: genome_sm_config

    use rule chromosome_sizes from genome_sm as ncRNAseq_chromosome_sizes
    use rule extract_smallrna from genome_sm as ncRNAseq_extract_smallrna

    # ── STAR index for smallRNA FASTA (reuses star module) ─────────────
    star_smallrna_idx_config = {
        "ROOT_DIR": ROOT_DIR,
        "indir": f"{outdir}/genome/smallrna",
        "outdir": f"{outdir}/genome/smallrna",
        "logdir": logdir,
        "Procedure": {"STAR": STAR},
        "Params": {"STAR": {}},
        "genome": {
            "fasta": smallrna_fasta,
            "gtf": None,
        }
    }
    logger.info(f"star_smallrna_idx_config: {star_smallrna_idx_config}")

    module star_smallrna_idx:
        snakefile: "../modules/star/star.smk"
        config: star_smallrna_idx_config

    use rule star_index from star_smallrna_idx as ncRNAseq_star_index_smallrna

    smallrna_star_index = f"{outdir}/genome/smallrna/index"

    # ── Read three-pass params from config (with defaults) ──────────────
    module star_pass1:
        snakefile: "../modules/star/star.smk"
        config: star_pass1_config
    module star_pass2:
        snakefile: "../modules/star/star.smk"
        config: star_pass2_config
    module star_pass3a:
        snakefile: "../modules/star/star.smk"
        config: star_pass3a_config
    module star_pass3b:
        snakefile: "../modules/star/star.smk"
        config: star_pass3b_config

    use rule star_align from star_pass1 as ncRNAseq_star3p_pass1
    use rule star_align from star_pass2 as ncRNAseq_star3p_pass2
    use rule star_align from star_pass3a as ncRNAseq_star3p_pass3a
    use rule star_align from star_pass3b as ncRNAseq_star3p_pass3b

    # ── Import auxiliary rules (extract, merge) from star_3pass ──────────
    star_3pass_config = {
        "ROOT_DIR": ROOT_DIR,
        "outdir": f"{outdir}/common/3_raw_bam",
        "pass2_outdir": f"{outdir}/common/3_raw_bam/pass2",
        "logdir": logdir,
        "Procedure": {
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "genome": {
            "smallrna_bed": smallrna_bed
        }
    }
    logger.info(f"star_3pass_config: {star_3pass_config}")

    module star_3pass:
        snakefile: "../modules/star/star_3pass/star_3pass.smk"
        config: star_3pass_config

    use rule star_3p_extract_smallrna from star_3pass as ncRNAseq_star3p_extract_smallrna
    use rule star_3p_pass2_to_fq from star_3pass as ncRNAseq_star3p_pass2_to_fq
    use rule star_3p_pass2_unmapped from star_3pass as ncRNAseq_star3p_pass2_unmapped
    use rule star_3p_pass3a_extract from star_3pass as ncRNAseq_star3p_pass3a_extract
    use rule star_3p_merge from star_3pass as ncRNAseq_star3p_merge

else:
    raise ValueError(f"Unsupported aligner: {aligner}. Please choose 'hisat2', 'star', or 'star_3pass'.")

# ── 3. Quantify (featureCounts) ──────────────────────────────────────────────
featureCounts_config = {
    "ROOT_DIR": ROOT_DIR,
    "indir": f"{outdir}/ncRNAseq/bam",
    "outdir": f"{outdir}/ncRNAseq/counts",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "featureCounts": config.get("Procedure", {}).get("featureCounts")
    },
    "genome": {
        "gtf": config.get("genome", {}).get("gtf")
    }
}
logger.info(f"featureCounts_config: {featureCounts_config}")
module featureCounts:
    snakefile: "../modules/featureCounts/featureCounts.smk"
    config: featureCounts_config
use rule featureCounts_single_noMultiple from featureCounts as ncRNAseq_featureCounts_single
use rule featureCounts_paired_noMultiple from featureCounts as ncRNAseq_featureCounts_paired
use rule featureCounts_result from featureCounts as ncRNAseq_featureCounts_result

# ── 4. Tailer (3' end analysis) ──────────────────────────────────────────────
tailer_config = {
    "ROOT_DIR": ROOT_DIR,
    "indir": f"{outdir}/common/3_raw_bam",
    "outdir": f"{outdir}/results/tailer",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "tailer": "Tailer"
    },
    "Params": {
        "tailer": config.get("Params", {}).get("tailer", {})
    },
    "genome": {
        "gtf": config.get("genome", {}).get("gtf")
    }
}
logger.info(f"tailer_config: {tailer_config}")
module tailer:
    snakefile: "../modules/tailer/tailer.smk"
    config: tailer_config
use rule tailer_global from tailer as ncRNAseq_tailer_global
