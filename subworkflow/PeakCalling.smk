shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outfiles = config.get("outfiles", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
samples = config.get("samples", [])
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})
genomes = config.get("genomes", ["mm"])

rule all:
    input:
        outfiles

# ============================================
# Step 1: Trimming with cutadapt/trim_galore
# ============================================
cutadapt_config = {
    "indir": indir,
    "outdir": f"{outdir}/cutadapt",
    "logdir": logdir,
    "Procedure": {
        "trim_galore": config.get("Procedure", {}).get("trim_galore")
    },
    "Params": {
        "trim_galore": {
            "quality": config.get("Params", {}).get("trim_galore", {}).get("quality") or 30
        }
    }
}
module cutadapt:
    snakefile: "../modules/cutadapt/cutadapt.smk"
    config: cutadapt_config
logger.info(f"cutadapt_config: {cutadapt_config}")
use rule trimming_Paired from cutadapt as PeakCalling_trimming_Paired
use rule trimming_Single from cutadapt as PeakCalling_trimming_Single

# ============================================
# Step 2: Bowtie2 Index
# ============================================
for genome in genomes:
    bowtie2_index_config = {
        "indir": cutadapt_config["outdir"],
        "outdir": f"{outdir}/bowtie2/{genome}",
        "logdir": logdir,
        "Procedure": {
            "bowtie2-build": config.get("Procedure", {}).get("bowtie2-build"),
            "bowtie2": config.get("Procedure", {}).get("bowtie2")
        },
        "genome": {
            "fasta": config.get("genome", {}).get(genome, {}).get("fasta") or config.get("genome", {}).get("fasta")
        }
    }
    module bowtie2_index:
        snakefile: "../modules/bowtie2/bowtie2.smk"
        config: bowtie2_index_config
    logger.info(f"bowtie2_index_config_{genome}: {bowtie2_index_config}")
    use rule bowtie2_index from bowtie2_index as PeakCalling_bowtie2_index

# ============================================
# Step 3: Bowtie2 Align
# ============================================
for genome in genomes:
    bowtie2_align_config = {
        "indir": cutadapt_config["outdir"],
        "outdir": f"{outdir}/bowtie2/{genome}",
        "logdir": logdir,
        "paired_samples": paired_samples,
        "single_samples": single_samples,
        "Procedure": {
            "bowtie2": config.get("Procedure", {}).get("bowtie2"),
            "samtools": config.get("Procedure", {}).get("samtools")
        },
        "Params": {
            "bowtie2": config.get("Params", {}).get("bowtie2", {})
        },
        "genome": {
            "fasta": config.get("genome", {}).get(genome, {}).get("fasta") or config.get("genome", {}).get("fasta"),
            "index_prefix": f"{outdir}/bowtie2/{genome}/index/genome"
        }
    }
    module bowtie2_align:
        snakefile: "../modules/bowtie2/bowtie2.smk"
        config: bowtie2_align_config
    logger.info(f"bowtie2_align_config_{genome}: {bowtie2_align_config}")
    use rule bowtie2_align_paired from bowtie2_align as PeakCalling_bowtie2_align_paired
    use rule bowtie2_align_single from bowtie2_align as PeakCalling_bowtie2_align_single

# ============================================
# Step 4: MACS3 Peak Calling
# ============================================
for genome in genomes:
    macs3_config = {
        "indir": f"{outdir}/bowtie2/{genome}",
        "outdir": f"{outdir}/macs3/{genome}",
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
            "fasta": config.get("genome", {}).get(genome, {}).get("fasta") or config.get("genome", {}).get("fasta")
        }
    }
    module macs3:
        snakefile: "../modules/macs3/macs3.smk"
        config: macs3_config
    logger.info(f"macs3_config_{genome}: {macs3_config}")
    use rule macs3_callpeak from macs3 as PeakCalling_macs3_callpeak
