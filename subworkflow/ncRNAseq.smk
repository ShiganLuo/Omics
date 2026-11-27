shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os

ROOT_DIR = config.get("ROOT_DIR", ".")
outfiles = config.get("outfiles", [])
indir = config.get("Params", {}).get("workflow").get("indir", "data/fastq")
outdir = config.get("Params", {}).get("workflow").get("outdir", "output")
logdir = config.get("Params", {}).get("workflow").get("logdir", "logs")
paired_samples = config.get("sample_params", {}).get("paired_samples", [])
single_samples = config.get("sample_params", {}).get("single_samples", [])
all_samples = config.get("sample_params", {}).get("samples", [])
aligner = config.get("Params", {}).get("workflow").get("aligner") or "star"

ip_samples = config.get("sample_params", {}).get("ip_samples", [])
input_samples = config.get("sample_params", {}).get("input_samples", [])
sample_ip_input_map = config.get("sample_params", {}).get("sample_ip_input_map", {})
genome = config.get("genome", {}).get("default")

rule all:
    input:
        outfiles
fastqc_raw_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": indir,
        "outdir":  f"{outdir}/QC/1_raw_fastqc",
        "logdir": os.path.join(logdir,"sample"),
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
use rule fastqc from fastqc_raw as ncRNAseq_fastqc_raw

if aligner in ["stat_3pass", "star_3pass_gene"]:
    demultiplexer_config = {
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "indir": indir,
            "outdir": f"{outdir}/common/2_trimmed_dedup_fastq/jla-demultiplexer",
            "logdir": os.path.join(logdir,"sample"),
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
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": demultiplexer_config["outdir"] if aligner in ["stat_3pass", "star_3pass_gene"] else indir,
        "outdir": f"{outdir}/common/2_trimmed_dedup_fastq/final_trimmed_fastq",
        "logdir": os.path.join(logdir,"sample"),
        "Procedure": {
            "trim_galore": config.get('Procedure',{}).get('trim_galore')
        },
        "Params": {
            "trim_galore": config.get('Params',{}).get("trim_galore", {})
        },
    }
module trim_galore:
    snakefile: "../modules/trim-galore/trim-galore.smk"
    config: trim_galore_config
logger.info(f"TrimGalore parameters: {trim_galore_config}")
use rule trimming_Paired from trim_galore as ncRNAseq_trimming_Paired
use rule trimming_Single from trim_galore as ncRNAseq_trimming_Single

if not config.get("Params", {}).get("ncRNAseq", {}).get("subsample", False):
    logger.info("Subsampling is disabled. Skipping subsample step.")
    final_fastq_outdir = trim_galore_config["outdir"]
else:
    logger.info("Subsampling is enabled. Proceeding with subsample step.")
    subsample_config = {
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "indir": trim_galore_config["outdir"],
            "outdir": f"{outdir}/common/2_trimmed_dedup_fastq/trimmed_subsampled_fastq",
            "logdir": os.path.join(logdir,"sample"),
            "paired_samples": paired_samples,
            "single_samples": single_samples,
            "Params": {
                "subsample": {
                    "keep_intermediates": config.get('Params',{}).get('subsample',{}).get('keep_intermediates', True),
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
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
        "outdir":  f"{outdir}/QC/2_trimmed_fastqc",
        "logdir": os.path.join(logdir,"sample"),
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
use rule fastqc from fastqc_trimmed as ncRNAseq_fastqc_trimmed

STAR = config.get("Procedure", {}).get("STAR") or "STAR"
SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
star_index_dir = config.get("genome", {}).get("references", {}).get(genome, {}).get("star_index_dir")

if aligner == "hisat2":
    hisat2_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
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
            "fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("fasta"),
            "hisat2_index_prefix": config.get("genome", {}).get("references", {}).get(genome, {}).get("hisat2_index_prefix")
        },
        "Params": {
            "hisat2": {
                "score_min": config.get("Params", {}).get("hisat2", {}).get("score_min") or "L,0,-0.6",
                "no-spliced-alignment": config.get("Params", {}).get("hisat2", {}).get("no-spliced-alignment") or True,
                "k": config.get("Params", {}).get("hisat2", {}).get("k") or 99999,
            }
        }
    }
    logger.info(f"hisat2_config: {hisat2_config}")
    module hisat2:
        snakefile: "../modules/hisat2/hisat2.smk"
        config: hisat2_config
    use rule hisat2_align_ncRNAseq_single from hisat2 as ncRNAseq_hisat2_align
    gatk_prepare_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": hisat2_config["outdir"],
        "outdir": f"{outdir}/common/4_markdup_bam",
        "logdir": os.path.join(logdir,"sample"),
        "input_bam_substring": "",
        "Procedure": {
            "gatk": config.get("Procedure", {}).get("gatk") or "gatk",
            "samtools": config.get("Procedure", {}).get("samtools") or "samtools"
        },
        "Params": {
            "gatk": config.get("Params", {}).get("gatk", {}),
            "gatk_prepare": {
                "keep_intermediates": config.get("Params", {}).get("gatk_prepare", {}).get("keep_intermediates", True)
            }
        },
        "genome": {
            "fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("fasta")
        }
    }
    module gatk_prepare:
        snakefile: "../modules/gatk/gatk_prepare.smk"
        config: gatk_prepare_config
    logger.info(f"gatk_prepare_config: {gatk_prepare_config}")
    use rule addReadsGroup from gatk_prepare as ncRNAseq_addReadsGroup
    use rule MarkDuplicates from gatk_prepare as ncRNAseq_MarkDuplicates
elif aligner == "star":
    if not star_index_dir:
        star_genome_idx_config = {
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/whole_genome",
            "logdir": os.path.join(logdir,"group"),
            "logdir_index": os.path.join(logdir,"group"),
            "log_index_substring": "star_index_genome",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star", {}).get("index", {})}},
            "genome": {
                "fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("fasta"),
                "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf"),
            }
        }
        logger.info(f"star_genome_idx_config: {star_genome_idx_config}")

        module star_genome_idx:
            snakefile: "../modules/star/star.smk"
            config: star_genome_idx_config

        use rule star_index from star_genome_idx as ncRNAseq_star_index_genome

        star_index_dir = f"{outdir}/genome/whole_genome/index"

    star_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
        "outdir": f"{outdir}/common/3_raw_bam",
        "logdir": os.path.join(logdir,"sample"),
        "logdir_index": os.path.join(logdir,"group"),
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "STAR": STAR
        },
        "Params": {
            "star": {
                "star_align": {
                    "twopassMode": config.get('Params', {}).get('star', {}).get('star_align', {}).get('twopassMode'),
                    "genomeLoad": config.get('Params', {}).get('star', {}).get('star_align', {}).get('genomeLoad') or "LoadAndRemove",
                    "limitBAMsortRAM": config.get('Params', {}).get('star', {}).get('star_align', {}).get('limitBAMsortRAM') or 20000000000,
                    "outReadsUnmapped": config.get('Params', {}).get('star', {}).get('star_align', {}).get('outReadsUnmapped') or "Fastx",
                    "outFilterMultimapNmax": config.get('Params', {}).get('star', {}).get('star_align', {}).get('outFilterMultimapNmax') or 99999,
                    "outFilterMismatchNoverLmax": config.get('Params', {}).get('star', {}).get('star_align', {}).get('outFilterMismatchNoverLmax') or 0.1,
                    "outFilterMatchNminOverLread": config.get('Params', {}).get('star', {}).get('star_align', {}).get('outFilterMatchNminOverLread') or 0.66,
                    "alignSJoverhangMin": config.get('Params', {}).get('star', {}).get('star_align', {}).get('alignSJoverhangMin') or 999,
                    "alignSJDBoverhangMin": config.get('Params', {}).get('star', {}).get('star_align', {}).get('alignSJDBoverhangMin') or 999,
                }
            }
        },
        "genome": {
            "fasta": config.get("genome",{}).get("references", {}).get(genome, {}).get("fasta"),
            "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf"),
            "index_dir": star_index_dir
        }
    }
    logger.info(f"star_config: {star_config}")
    module star:
        snakefile: "../modules/star/star.smk"
        config: star_config
    use rule star_align from star as ncRNAseq_star_align
    gatk_prepare_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": star_config["outdir"],
        "outdir": f"{outdir}/common/4_markdup_bam",
        "logdir": os.path.join(logdir,"sample"),
        "input_bam_substring": "",
        "Procedure": {
            "gatk": config.get("Procedure", {}).get("gatk") or "gatk",
            "samtools": config.get("Procedure", {}).get("samtools") or "samtools"
        },
        "Params": {
            "gatk": config.get("Params", {}).get("gatk", {}),
            "gatk_prepare": {
                "keep_intermediates": config.get("Params", {}).get("gatk_prepare", {}).get("keep_intermediates", True)
            }
        },
        "genome": {
            "fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("fasta")
        }
    }
    module gatk_prepare:
        snakefile: "../modules/gatk/gatk_prepare.smk"
        config: gatk_prepare_config
    logger.info(f"gatk_prepare_config: {gatk_prepare_config}")
    use rule addReadsGroup from gatk_prepare as ncRNAseq_addReadsGroup
    use rule MarkDuplicates from gatk_prepare as ncRNAseq_MarkDuplicates
elif aligner == "star_3pass":

    smallrna_bed = f"{outdir}/genome/smallrna/smallrna_genes.bed"
    smallrna_fasta = f"{outdir}/genome/smallrna/smallrna_genes_flank.fa"
    smallrna_star_index_dir = f"{outdir}/genome/smallrna/index"

    if not star_index_dir or not os.path.exists(star_index_dir):
        star_genome_idx_config = {
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/whole_genome",
            "logdir": os.path.join(logdir,"sample"),
            "logdir_index": os.path.join(logdir,"group"),
            "log_index_substring": "star_index_genome",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("genome", {})}},
            "genome": {
                "fasta": config.get("genome",{}).get("references", {}).get(genome, {}).get("fasta"),
                "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf"),
            }
        }
        logger.info(f"star_genome_idx_config: {star_genome_idx_config}")

        module star_genome_idx:
            snakefile: "../modules/star/star.smk"
            config: star_genome_idx_config

        use rule star_index from star_genome_idx as ncRNAseq_star_index_genome

        star_index_dir = f"{outdir}/genome/whole_genome/index"

    genome_sm_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "outdir": outdir,
        "logdir": os.path.join(logdir,"group"),
        "Procedure": {
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "Params": {
            "smallrna_types": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_types") or ["miRNA", "snRNA", "snoRNA", "rRNA", "misc_RNA", "scRNA", "scaRNA", "vaultRNA"],
            "smallrna_flank": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_flank") or 50
        },
        "genome": {
            "fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("fasta"),
            "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf")
        }
    }
    logger.info(f"genome_sm_config: {genome_sm_config}")

    module genome_sm:
        snakefile: "../modules/genome/genome.smk"
        config: genome_sm_config

    use rule chromosome_sizes from genome_sm as ncRNAseq_chromosome_sizes
    use rule extract_smallrna from genome_sm as ncRNAseq_extract_smallrna

    smallrna_star_index_dir = config.get("genome", {}).get("references", {}).get(genome, {}).get("smallrna_star_index_dir")
    if not smallrna_star_index_dir or not os.path.exists(smallrna_star_index_dir):
        star_smallrna_idx_config = {
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "indir": f"{outdir}/genome/smallrna",
            "outdir": f"{outdir}/genome/smallrna",
            "logdir": os.path.join(logdir,"sample"),
            "logdir_index": os.path.join(logdir,"group"),
            "log_index_substring": "star_index_smallrna",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("smallrna", {})}},
            "genome": {
                "fasta": smallrna_fasta,
                "gtf": None,
                "is_fasta_intermediate": True
            }
        }
        logger.info(f"star_smallrna_idx_config: {star_smallrna_idx_config}")
        module star_smallrna_idx:
            snakefile: "../modules/star/star.smk"
            config: star_smallrna_idx_config

        use rule star_index from star_smallrna_idx as ncRNAseq_star_index_smallrna

        smallrna_star_index_dir = f"{outdir}/genome/smallrna/index"

    star_3pass_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
        "outdir": f"{outdir}/common/3_raw_bam/final_bam",
        "logdir": os.path.join(logdir,"sample"),
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
    tailer_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": star_3pass_config["outdir"],
        "outdir": f"{outdir}/results/tailer",
        "logdir": os.path.join(logdir,"sample"),
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "tailer": config.get("Procedure", {}).get("tailer") or config.get("Procedure", {}).get("Tailer") or "Tailer"
        },
        "Params": {
            "tailer": config.get("Params", {}).get("tailer", {})
        },
        "genome": {
            "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf")
        }
    }
    logger.info(f"tailer_config: {tailer_config}")
    module tailer:
        snakefile: "../modules/tailer/tailer.smk"
        config: tailer_config
    use rule tailer_global from tailer as ncRNAseq_tailer_global
elif aligner == "star_3pass_gene":

    smallrna_bed = f"{outdir}/genome/smallrna/smallrna_genes.bed"
    smallrna_fasta = f"{outdir}/genome/smallrna/smallrna_genes_flank.fa"
    smallrna_star_index_dir = f"{outdir}/genome/smallrna/index"

    if not star_index_dir or not os.path.exists(star_index_dir):
        star_genome_idx_config = {
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/whole_genome",
            "logdir": os.path.join(logdir,"sample"),
            "logdir_index": os.path.join(logdir,"group"),
            "log_index_substring": "star_index_genome",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("genome", {})}},
            "genome": {
                "fasta": config.get("genome",{}).get("references", {}).get(genome, {}).get("fasta"),
                "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf"),
            }
        }
        logger.info(f"star_genome_idx_config: {star_genome_idx_config}")

        module star_genome_idx:
            snakefile: "../modules/star/star.smk"
            config: star_genome_idx_config

        use rule star_index from star_genome_idx as ncRNAseq_star_index_genome

        star_index_dir = f"{outdir}/genome/whole_genome/index"

    genome_sm_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "outdir": outdir,
        "logdir": os.path.join(logdir,"group"),
        "Procedure": {
            "samtools": SAMTOOLS,
            "bedtools": BEDTOOLS,
        },
        "Params": {
            "smallrna_types": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_types") or ["snRNA", "misc_RNA", "rRNA", "rRNA_pseudogene", "snoRNA", "scaRNA", "ribozyme", "TERC"],
            "smallrna_flank": config.get("Params", {}).get("ncRNAseq", {}).get("smallrna_flank") or 50
        },
        "genome": {
            "fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("fasta"),
            "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf")
        }
    }
    logger.info(f"genome_sm_config: {genome_sm_config}")

    module genome_sm:
        snakefile: "../modules/genome/genome.smk"
        config: genome_sm_config

    use rule chromosome_sizes from genome_sm as ncRNAseq_chromosome_sizes
    use rule extract_smallrna from genome_sm as ncRNAseq_extract_smallrna

    configured_smallrna_index = config.get("genome", {}).get("references", {}).get(genome, {}).get("smallrna_star_index_dir")
    if configured_smallrna_index and os.path.exists(configured_smallrna_index):
        smallrna_star_index_dir = configured_smallrna_index
    elif not os.path.exists(smallrna_star_index_dir):
        star_smallrna_idx_config = {
            "resources": config.get("resources", {}),
            "ROOT_DIR": ROOT_DIR,
            "OUTDIR_ROOT": outdir,
            "env": config.get("env", {}),
            "outdir": f"{outdir}/genome/smallrna",
            "logdir": os.path.join(logdir,"sample"),
            "logdir_index": os.path.join(logdir,"group"),
            "log_index_substring": "star_index_smallrna",
            "Procedure": {"STAR": STAR},
            "Params": {"star": {"index": config.get("Params", {}).get("star_3pass", {}).get("index", {}).get("smallrna", {})}},
            "genome": {
                "fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("smallrna_fasta"),
                "gtf": None,
                "is_fasta_intermediate": True
            }
        }
        logger.info(f"star_smallrna_idx_config: {star_smallrna_idx_config}")
        module star_smallrna_idx:
            snakefile: "../modules/star/star.smk"
            config: star_smallrna_idx_config
        use rule star_index from star_smallrna_idx as ncRNAseq_star_index_smallrna

    p3g = config.get("Params", {}).get("star_3pass_gene", {})
    star_3pass_gene_upstream_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "indir": final_fastq_outdir,
        "outdir": f"{outdir}/common/3_raw_bam",
        "logdir": os.path.join(logdir,"sample"),
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

    star_3pass_gene_config = {
        "resources": config.get("resources", {}),
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "outdir": f"{outdir}/common/4_per_gene_bam",
        "final_bam_dir": star_3pass_gene_upstream_config["outdir"],
        "logdir": os.path.join(logdir,"sample"),
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
            "genome_fasta": config.get("genome", {}).get("references", {}).get(genome, {}).get("fasta"),
        }
    }
    logger.info(f"star_3pass_gene_config: {star_3pass_gene_config}")

    module star_3pass_gene:
        snakefile: "../modules/star/star_3pass/star_3pass_gene.smk"
        config: star_3pass_gene_config

    use rule star_3pg_gene_specific from star_3pass_gene as ncRNAseq_star3pg_gene_specific
else:
    raise ValueError(f"Unsupported aligner: {aligner}. Please choose 'hisat2', 'star', 'star_3pass', or 'star_3pass_gene'.")

if aligner in ("star_3pass", "star_3pass_gene"):
    align_bam_dir = f"{outdir}/common/3_raw_bam"
else:
    align_bam_dir = f"{outdir}/common/4_markdup_bam"

igv_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": align_bam_dir,
    "outdir": f"{outdir}/tracks",
    "logdir": os.path.join(logdir,"sample"),
    "bam_substring": "sorted_markdup",
    "Procedure": {
        "samtools": config.get("Procedure", {}).get("samtools") or "samtools",
        "bamCoverage": config.get("Procedure", {}).get("bamCoverage") or "bamCoverage"
    },
    "Params": {
        "bamCoverage": config.get("Params", {}).get("bamCoverage", {}),
    }
}
module igv:
    snakefile: "../modules/igv/igv.smk"
    config: igv_config
logger.info(f"igv_config: {igv_config}")
use rule samtools_dedup from igv as ncRNAseq_dedup
use rule wig from igv as ncRNAseq_bigwig

regionplot_params = config.get("Params", {}).get("regionplot", {})
regionplot_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/regionplot",
    "logdir": f"{logdir}/group",
    "samples": all_samples,
    "tracks": [{"label": "bigwig", "dir": igv_config["outdir"], "suffix": ".bigwig"}],
    "Params": {
        "regionplot": regionplot_params
    },
    "genome": {
        "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf"),
    },
    "resources": config.get("resources", {}),
}
module regionplot:
    snakefile: "../modules/regionplot/regionplot.smk"
    config: regionplot_config
logger.debug(f"regionplot_config: {regionplot_config}")
use rule regionplot from regionplot as ncRNAseq_regionplot
use rule regionplot_result from regionplot as ncRNAseq_regionplot_result

track_config = {
        "resources": config.get("resources", {}),
        "indir": igv_config["outdir"],
        "outdir":  igv_config["outdir"],
        "logdir": os.path.join(logdir,"group"),
        "samples": single_samples + paired_samples,
        "ROOT_DIR": ROOT_DIR,
        "OUTDIR_ROOT": outdir,
        "env": config.get("env", {}),
        "igv": config.get('Params', {}).get('track', {}).get('igv', {}).get('configs', {}).get(genome, {}),
    }

module track:
    snakefile: "../modules/track/track.smk"
    config: track_config
logger.info(f"track_config: {track_config}")
use rule * from track as ncRNAseq_*

featureCounts_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "indir": align_bam_dir,
    "outdir": f"{outdir}/counts",
    "logdir": os.path.join(logdir,"sample"),
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "bam_substring": "sorted_markdup",
    "Procedure": {
        "featureCounts": config.get("Procedure", {}).get("featureCounts")
    },
    "genome": {
        "gtf": config.get("genome", {}).get("references", {}).get(genome, {}).get("gtf")
    }
}
logger.info(f"featureCounts_config: {featureCounts_config}")
module featureCounts:
    snakefile: "../modules/featureCounts/featureCounts.smk"
    config: featureCounts_config
use rule featureCounts_single_noMultiple from featureCounts as ncRNAseq_featureCounts_single
use rule featureCounts_paired_noMultiple from featureCounts as ncRNAseq_featureCounts_paired
use rule featureCounts_result from featureCounts as ncRNAseq_featureCounts_result

ncRNAseq_report_config = {
    "resources": config.get("resources", {}),
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": outdir,
    "logdir": os.path.join(logdir, "sample"),
    "samples": all_samples,
    "paired_samples": paired_samples,
    "single_samples": single_samples,
    "aligner": aligner,
    "sample_groups": config.get("sample_groups", {}),
    "short_names": config.get("short_names", {}),
    "Params": {
        "report": config.get("Params", {}).get("report", {}),
    },
}
module ncRNAseq_report:
    snakefile: "../modules/ncRNAseq_report/ncRNAseq_report.smk"
    config: ncRNAseq_report_config
logger.info(f"ncRNAseq_report_config: {ncRNAseq_report_config}")
use rule generate_report from ncRNAseq_report as ncRNAseq_generate_report

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
        s: f"{gatk_prepare_config['outdir']}/{s}/{s}.sorted_markdup.bam"
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
use rule sexcheck_idxstats from sexcheck as ncRNAseq_sexcheck_idxstats
use rule sexcheck from sexcheck as ncRNAseq_sexcheck

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
        "markdup": gatk_prepare_config["outdir"],
        "sexcheck": sexcheck_config["outdir"],
    },
    "qc_anchors": {
        "markdup": "{sample}/{sample}.Markdup-metrics.txt",
        "sexcheck": "sexcheck.tsv",
    },
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
logger.info(f"qc_summary_config: {qc_summary_config}")
use rule qc_summary from qc_summary as ncRNAseq_qc_summary

use rule report_result from ncRNAseq_report as ncRNAseq_report_result

