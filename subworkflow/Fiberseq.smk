"""Fiber-seq subworkflow for single-molecule chromatin accessibility analysis.

Workflow:
  Step 1: fiberseq/1_prepare — predict-m6a (auto-skip if SPRQ) + add-nuc + fire
  Step 2: common/2_aligned  — pbmm2 align → aligned fiberseq BAM
  Step 3: fiberseq/3_extract — extract BED12 + call-peaks + qc

Reference:
  - Stergachis et al., 2020, Science (DOI: 10.1126/science.aaz1646)
  - Jha, Bohaczuk et al., 2024, Genome Research (DOI: 10.1101/gr.279095.124)
  - https://fiberseq.github.io/
"""
shell.prefix("set -x; set -e;")
from snakemake.logging import logger

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
outfiles = config.get("outfiles") or []
samples = config.get("samples") or []

# ---- Resolve genome reference ----
genome_cfg = config.get("genome", {})
default_genome = genome_cfg.get("default", "")
genome_ref = genome_cfg.get("references", {}).get(default_genome, {})
fasta = genome_ref.get("fasta") or genome_cfg.get("fasta") or ""
logger.info(f"Fiber-seq genome: {default_genome}, fasta={fasta}")

rule all:
    input:
        outfiles


# ============================================================
# Step 1: Prepare fiberseq BAM → common/2_fiberseq_bam
#   predict-m6a (auto-skip if SPRQ) → add-nucleosomes → fire
# ============================================================
prepare_outdir = f"{outdir}/common/2_fiberseq_bam"

fibertools_prepare_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": indir,
    "outdir": prepare_outdir,
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "Params": {
        "fibertools": config.get("Params", {}).get("fibertools", {})
    },
    "genome": {
        "fasta": fasta
    }
}
module fibertools_prepare:
    snakefile: "../modules/fibertools/fibertools_prepare.smk"
    config: fibertools_prepare_config
logger.debug(f"Fiber-seq prepare config: {fibertools_prepare_config}")
use rule ft_predict_m6a from fibertools_prepare as Fiberseq_ft_predict_m6a
use rule ft_add_nucleosomes from fibertools_prepare as Fiberseq_ft_add_nucleosomes
use rule ft_fire from fibertools_prepare as Fiberseq_ft_fire


# ============================================================
# Step 2: Align to genome → common/3_align_bam
# ============================================================
aligned_dir = f"{outdir}/common/3_align_bam"

pbmm2_config = {
    "indir": prepare_outdir,
    "outdir": aligned_dir,
    "bam_suffix": ".fiberseq.fire.bam",
    "output_bam_suffix": ".fiberseq.fire.sorted",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "Procedure": {
        "pbmm2": config.get("Procedure", {}).get("pbmm2"),
        "samtools": config.get("Procedure", {}).get("samtools"),
    },
    "genome": {
        "fasta": fasta
    }
}
module pbmm2:
    snakefile: "../modules/pbmm2/pbmm2.smk"
    config: pbmm2_config
logger.debug(f"Fiber-seq pbmm2 config: {pbmm2_config}")
use rule pbmm2_align from pbmm2 as Fiberseq_pbmm2_align


# ============================================================
# Step 3: Extract BED + peaks + QC → results
# ============================================================
analysis_outdir = f"{outdir}/results"

fibertools_analysis_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": analysis_outdir,
    "bam_suffix": ".fiberseq.fire.sorted.bam",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "Params": {
        "fibertools": config.get("Params", {}).get("fibertools", {})
    },
    "genome": {
        "fasta": fasta
    }
}
module fibertools_analysis:
    snakefile: "../modules/fibertools/fibertools_analysis.smk"
    config: fibertools_analysis_config
logger.debug(f"Fiber-seq analysis config: {fibertools_analysis_config}")
use rule ft_extract from fibertools_analysis as Fiberseq_ft_extract
use rule ft_call_peaks from fibertools_analysis as Fiberseq_ft_call_peaks
use rule ft_qc from fibertools_analysis as Fiberseq_ft_qc