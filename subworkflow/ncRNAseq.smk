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
        "env": config.get("env", {}),
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
        "env": config.get("env", {}),
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
        "env": config.get("env", {}),
        "indir": demultiplexer_config["outdir"],
        "outdir": f"{outdir}/common/2_trimmed_dedup_fastq/final_trimmed_fastq",
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
if not config.get("Params", {}).get("ncRNAseq", {}).get("subsample", False):
    logger.info("Subsampling is disabled. Skipping subsample step.")
    final_fastq_outdir = trim_galore_config["outdir"]
else:
    logger.info("Subsampling is enabled. Proceeding with subsample step.")
    subsample_config = {
            "ROOT_DIR": ROOT_DIR,
            "env": config.get("env", {}),
            "indir": trim_galore_config["outdir"],
            "outdir": f"{outdir}/common/2_trimmed_dedup_fastq/trimmed_subsampled_fastq",
            "logdir": logdir,
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Params": {
                "subsample": {
                    "abund_small_rnas": config.get("Params", {}).get("ncRNAseq", {}).get("abund_small_rnas", []),
                    "subsample_n": config.get("Params", {}).get("ncRNAseq", {}).get("subsample_n", 100000),
                    "subsample_seed": config.get("Params", {}).get("ncRNAseq", {}).get("subsample_seed", 42),
                    "hard_clip_direction": config.get("Params", {}).get("ncRNAseq", {}).get("hard_clip_direction", "5prime"),
                    "hard_clip_length": config.get("Params", {}).get("ncRNAseq", {}).get("hard_clip_length", 0)
                }
            },
        }
    module subsample:
        snakefile: "../modules/subsample/subsample.smk"
        config: subsample_config
    logger.info(f"subsample_config: {subsample_config}")
    use rule subsample_fastq from subsample as ncRNAseq_subsample_fastq
    final_fastq_outdir = subsample_config["outdir"]


fastqc_trimmed_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
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


if aligner == "hisat2":
    hisat2_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
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
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/whole_genome",
            "logdir": logdir,
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star", {}).get("index", {})}},
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

        star_index_dir = f"{outdir}/genome/whole_genome/index"

    star_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
        "outdir": f"{outdir}/common/3_raw_bam",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "STAR": STAR
        },
        "Params": {
            "star": {
                "genomeLoad": config.get("Params", {}).get("star", {}).get("genomeLoad") or "LoadAndRemove",
                "limitBAMsortRAM": config.get("Params", {}).get("star", {}).get("limitBAMsortRAM") or 20000000000,
                "outReadsUnmapped": config.get("Params", {}).get("star", {}).get("outReadsUnmapped") or "Fastx",
                "outFilterMultimapNmax": config.get("Params", {}).get("star", {}).get("outFilterMultimapNmax") or 99999,
                "outFilterMismatchNoverLmax": config.get("Params", {}).get("star", {}).get("outFilterMismatchNoverLmax") or 0.1,
                "outFilterMatchNminOverLread": config.get("Params", {}).get("star", {}).get("outFilterMatchNminOverLread") or 0.66,
                "alignSJoverhangMin": config.get("Params", {}).get("star", {}).get("alignSJoverhangMin") or 999,
                "alignSJDBoverhangMin": config.get("Params", {}).get("star", {}).get("alignSJDBoverhangMin") or 999
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
    smallrna_star_index_dir = f"{outdir}/genome/smallrna/index"

    # ── STAR index for genome (auto-build when star_index_dir is null) ──
    if not star_index_dir or not os.path.exists(star_index_dir):
        star_genome_idx_config = {
            "ROOT_DIR": ROOT_DIR,
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/whole_genome",
            "logdir": f"{logdir}/genome/whole_genome",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("genome", {})}},
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

        star_index_dir = f"{outdir}/genome/whole_genome/index"

    # ── Import genome module (extract smallRNA BED/FASTA) ──────────────
    genome_sm_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "outdir": outdir,
        "logdir": f"{logdir}/genome",
        "Procedure": {
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "Params": {
            "smallrna_types": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_types") or ["miRNA", "snRNA", "snoRNA", "rRNA", "misc_RNA", "scRNA", "scaRNA", "vaultRNA"],
            "smallrna_flank": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_flank") or 50
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
    smallrna_star_index_dir = config.get("genome", {}).get("smallrna_star_index_dir")
    if not smallrna_star_index_dir or not os.path.exists(smallrna_star_index_dir):
        star_smallrna_idx_config = {
            "ROOT_DIR": ROOT_DIR,
            "env": config.get("env", {}),
            "indir": f"{outdir}/genome/smallrna",
            "outdir": f"{outdir}/genome/smallrna",
            "logdir": f"{logdir}/genome/smallrna",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("smallrna", {})}},
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

        smallrna_star_index_dir = f"{outdir}/genome/smallrna/index"

    # ── Canonical three-pass alignment in one execution rule ─────────────
    star_3pass_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
        "outdir": f"{outdir}/common/3_raw_bam/final_bam",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "STAR": STAR,
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "genome": {
            "genome_index": star_index_dir,
            "smallrna_index": smallrna_star_index_dir,
            "smallrna_bed": smallrna_bed,
        }
    }
    logger.info(f"star_3pass_config: {star_3pass_config}")
    module star_3pass:
        snakefile: "../modules/star/star_3pass/star_3pass.smk"
        config: star_3pass_config
    use rule star_3p_align from star_3pass as ncRNAseq_star3p_align

elif aligner == "star_3pass_gene":
    # ================================================================
    # Three-pass STAR alignment with per-gene re-alignment
    # (Ma et al, 2024 modification)
    #
    # Pass 1:  genome end-to-end alignment with 5' hard clip (10 nt)
    #          -> facilitates alignment of reads with post-transcriptional
    #             nucleotide modifications
    # Extract: group aligned reads by gene, reconvert to per-gene FASTQ
    #          using bedtools + samtools
    # Pass 2:  per-gene local alignment to single-gene genomic sequences
    # Merge:   combine per-gene BAMs into sample BAM
    # Tailer:  gene-specific 3'-end info in global alignment mode
    # ================================================================

    # ── Resolve smallRNA derived paths ──────────────────────────────
    smallrna_bed = f"{outdir}/genome/smallrna/smallrna_genes.bed"
    smallrna_fasta = f"{outdir}/genome/smallrna/smallrna_genes_flank.fa"
    smallrna_star_index_dir = f"{outdir}/genome/smallrna/index"

    # ── STAR index for genome (auto-build when star_index_dir is null) ──
    if not star_index_dir or not os.path.exists(star_index_dir):
        star_genome_idx_config = {
            "ROOT_DIR": ROOT_DIR,
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/whole_genome",
            "logdir": f"{logdir}/genome/whole_genome",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("genome", {})}},
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

        star_index_dir = f"{outdir}/genome/whole_genome/index"

    # ── Import genome module (extract smallRNA BED/FASTA) ──────────────
    genome_sm_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "outdir": outdir,
        "logdir": f"{logdir}/genome",
        "Procedure": {
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "Params": {
            "smallrna_types": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_types") or ["snRNA", "misc_RNA", "rRNA", "rRNA_pseudogene", "snoRNA", "scaRNA", "ribozyme", "TERC"],
            "smallrna_flank": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_flank") or 50
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

    # The upstream three-pass method still requires the canonical small-RNA
    # reference index. Build it from the extracted branch-specific FASTA when
    # an existing index was not supplied.
    configured_smallrna_index = config.get("genome", {}).get("smallrna_star_index_dir")
    if configured_smallrna_index and os.path.exists(configured_smallrna_index):
        smallrna_star_index_dir = configured_smallrna_index
    elif not os.path.exists(smallrna_star_index_dir):
        star_smallrna_idx_config = {
            "ROOT_DIR": ROOT_DIR,
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/smallrna",
            "logdir": f"{logdir}/genome/smallrna",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("smallrna", {})}},
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

    # The branch reuses the canonical three-pass implementation, with the
    # paper's end-to-end adaptation made explicit.
    p3g = config.get("Params", {}).get("star_3pass_gene", {})
    star_3pass_gene_upstream_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
        "outdir": f"{outdir}/common/3_raw_bam",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "STAR": STAR,
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "Params": {
            "star_3pass": config.get("Params", {}).get("star_3pass", {}),
        },
        "genome": {
            "genome_index": star_index_dir,
            "smallrna_index": smallrna_star_index_dir,
            "smallrna_bed": smallrna_bed,
        }
    }
    module star_3pass_gene_upstream:
        snakefile: "../modules/star/star_3pass/star_3pass.smk"
        config: star_3pass_gene_upstream_config
    use rule star_3p_align from star_3pass_gene_upstream as ncRNAseq_star3pg_upstream_align

    # ── Strict post-final-BAM gene-specific local alignment ──────────────
    star_3pass_gene_config = {
        "ROOT_DIR": ROOT_DIR,
        "env": config.get("env", {}),
        "outdir": f"{outdir}/common/4_per_gene_bam",
        "final_bam_dir": star_3pass_gene_upstream_config["outdir"],
        "logdir": logdir,
        "Procedure": {
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
            "STAR": STAR,
            "tailer": config.get("Procedure", {}).get("tailer") or config.get("Procedure", {}).get("Tailer") or "Tailer",
        },
        "Params": {
            "star_3pass_gene": p3g,
            "tailer": config.get("Params", {}).get("tailer", {}),
        },
        "genome": {
            "smallrna_bed": smallrna_bed,
            "genome_fasta": genome_fasta,
        }
    }
    logger.info(f"star_3pass_gene_config: {star_3pass_gene_config}")

    module star_3pass_gene:
        snakefile: "../modules/star/star_3pass/star_3pass_gene.smk"
        config: star_3pass_gene_config

    use rule star_3pg_gene_specific from star_3pass_gene as ncRNAseq_star3pg_gene_specific

else:
    raise ValueError(f"Unsupported aligner: {aligner}. Please choose 'hisat2', 'star', 'star_3pass', or 'star_3pass_gene'.")

# ── 3. Quantify (featureCounts) ──────────────────────────────────────────────
# Resolve the alignment output directory based on the chosen aligner.
if aligner in ("star_3pass", "star_3pass_gene"):
    # Gene-local BAM coordinates are incompatible with the whole-genome GTF.
    # Keep featureCounts attached to the canonical final BAM; gene-specific
    # BAMs are consumed only with their matching local annotations by Tailer.
    align_bam_dir = f"{outdir}/common/3_raw_bam/final_bam"
else:
    align_bam_dir = f"{outdir}/common/3_raw_bam"

featureCounts_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": align_bam_dir,
    "outdir": f"{outdir}/counts",
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
if aligner == "star_3pass_gene":
    align_bam_dir = f"{outdir}/common/3_raw_bam/per_gene"
    tailer_gtf = config.get("genome", {}).get("gtf")
else:
    tailer_gtf = config.get("genome", {}).get("gtf")
tailer_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": align_bam_dir,
    "outdir": f"{outdir}/results/tailer",
    "logdir": logdir,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Procedure": {
        "tailer": config.get("Procedure", {}).get("tailer") or config.get("Procedure", {}).get("Tailer") or "Tailer"
    },
    "Params": {
        "tailer": config.get("Params", {}).get("tailer", {})
    },
    "genome": {
        "gtf": tailer_gtf
    }
}
logger.info(f"tailer_config: {tailer_config}")
module tailer:
    snakefile: "../modules/tailer/tailer.smk"
    config: tailer_config
if aligner != "star_3pass_gene":
    use rule tailer_global from tailer as ncRNAseq_tailer_global

# ── 5. Report ────────────────────────────────────────────────────────────────
ncRNAseq_report_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "outdir": outdir,
    "logdir": logdir,
    "samples": all_samples,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "Params": {
        "report": config.get("Params", {}).get("report", {}),
    },
}
module ncRNAseq_report:
    snakefile: "../modules/ncRNAseq_report/ncRNAseq_report.smk"
    config: ncRNAseq_report_config
logger.info(f"ncRNAseq_report_config: {ncRNAseq_report_config}")
use rule generate_report from ncRNAseq_report as ncRNAseq_generate_report
use rule report_result from ncRNAseq_report as ncRNAseq_report_result
