"""Fiber-seq subworkflow for single-molecule chromatin accessibility analysis.

Workflow:
  common/1_raw_data  — raw CCS BAM (set by MetaUtil)
  common/2_aligned   — pbmm2 aligned BAM (shared preprocessing)
  fiberseq/1_nucleosomes — ft add-nucleosomes
  fiberseq/2_fire        — ft fire
  fiberseq/3_extract     — ft extract --all
  fiberseq/4_peaks       — ft call-peaks
  fiberseq/5_qc          — ft qc

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
# Step 1: Align CCS reads (pbmm2 module) → common/2_aligned
# ============================================================
aligned_dir = f"{outdir}/common/2_aligned"

pbmm2_config = {
    "indir": indir,
    "outdir": aligned_dir,
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
# Step 2: Add nucleosome calls → fiberseq/1_nucleosomes
# ============================================================
nuc_outdir = f"{outdir}/fiberseq/1_nucleosomes"

fibertools_nuc_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": aligned_dir,
    "outdir": nuc_outdir,
    "upstream_outdir": aligned_dir,
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "bam_suffix": ".sorted.bam",
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "genome": {
        "fasta": fasta
    }
}
module fibertools_nuc:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_nuc_config
logger.debug(f"Fiber-seq nucleosome config: {fibertools_nuc_config}")
use rule ft_add_nucleosomes from fibertools_nuc as Fiberseq_ft_add_nucleosomes


# ============================================================
# Step 3: FIRE calling → fiberseq/2_fire
# ============================================================
fire_outdir = f"{outdir}/fiberseq/2_fire"

fibertools_fire_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": nuc_outdir,
    "outdir": fire_outdir,
    "upstream_outdir": nuc_outdir,
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
module fibertools_fire:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_fire_config
logger.debug(f"Fiber-seq FIRE config: {fibertools_fire_config}")
use rule ft_fire from fibertools_fire as Fiberseq_ft_fire


# ============================================================
# Step 4: Extract data to BED → fiberseq/3_extract
# ============================================================
extract_outdir = f"{outdir}/fiberseq/3_extract"

fibertools_extract_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fire_outdir,
    "outdir": extract_outdir,
    "upstream_outdir": fire_outdir,
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "genome": {
        "fasta": fasta
    }
}
module fibertools_extract:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_extract_config
logger.debug(f"Fiber-seq extract config: {fibertools_extract_config}")
use rule ft_extract from fibertools_extract as Fiberseq_ft_extract


# ============================================================
# Step 5: FIRE peak calling → fiberseq/4_peaks
# ============================================================
peaks_outdir = f"{outdir}/fiberseq/4_peaks"

fibertools_peaks_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fire_outdir,
    "outdir": peaks_outdir,
    "upstream_outdir": fire_outdir,
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
module fibertools_peaks:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_peaks_config
logger.debug(f"Fiber-seq peaks config: {fibertools_peaks_config}")
use rule ft_call_peaks from fibertools_peaks as Fiberseq_ft_call_peaks


# ============================================================
# Step 6: QC metrics → fiberseq/5_qc
# ============================================================
qc_outdir = f"{outdir}/fiberseq/5_qc"

fibertools_qc_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fire_outdir,
    "outdir": qc_outdir,
    "upstream_outdir": fire_outdir,
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
module fibertools_qc:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_qc_config
logger.debug(f"Fiber-seq QC config: {fibertools_qc_config}")
use rule ft_qc from fibertools_qc as Fiberseq_ft_qc
