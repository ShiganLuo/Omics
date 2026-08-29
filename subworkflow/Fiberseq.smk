"""Fiber-seq subworkflow for single-molecule chromatin accessibility analysis.

Workflow:
  1. Predict m6A modifications (ft predict-m6a)
  2. Add nucleosome calls (ft add-nucleosomes)
  3. Call FIREs (ft fire)
  4. Extract data to BED format (ft extract)

Optional:
  - Align with pbmm2 (if input is unaligned CCS BAM)
  - Phase with HiPhase (for allele-specific analysis)

Reference:
  - Stergachis et al., 2020, Science (Fiber-seq original paper, DOI: 10.1126/science.aaz1646)
  - Jha, Bohaczuk et al., 2024, Genome Research (fibertools-rs, DOI: 10.1101/gr.279095.124)
  - https://fiberseq.github.io/
"""
from snakemake.logging import logger

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
outfiles = config.get("outfiles") or []
samples = config.get("samples") or []
skip_align = config.get("Params", {}).get("skip_align", False)
skip_phase = config.get("Params", {}).get("skip_phase", False)

rule all:
    input:
        outfiles


# ============================================================
# Step 1: m6A prediction (if input is PacBio CCS BAM with kinetics)
# ============================================================
fibertools_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": indir,
    "outdir": f"{outdir}/fiberseq/1_m6a",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module fibertools_m6a:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_config
logger.info(f"Fiber-seq m6a config: {fibertools_config}")
use rule ft_predict_m6a from fibertools_m6a as Fiberseq_ft_predict_m6a


# ============================================================
# Step 2: Add nucleosome calls
# ============================================================
fibertools_nuc_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fibertools_config["outdir"],
    "outdir": f"{outdir}/fiberseq/2_nucleosomes",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module fibertools_nuc:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_nuc_config
logger.info(f"Fiber-seq nucleosome config: {fibertools_nuc_config}")
use rule ft_add_nucleosomes from fibertools_nuc as Fiberseq_ft_add_nucleosomes


# ============================================================
# Step 3: FIRE calling
# ============================================================
fibertools_fire_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fibertools_nuc_config["outdir"],
    "outdir": f"{outdir}/fiberseq/3_fire",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "Params": {
        "fibertools": config.get("Params", {}).get("fibertools", {})
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module fibertools_fire:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_fire_config
logger.info(f"Fiber-seq FIRE config: {fibertools_fire_config}")
use rule ft_fire from fibertools_fire as Fiberseq_ft_fire


# ============================================================
# Step 4: Extract data to BED format
# ============================================================
fibertools_extract_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fibertools_fire_config["outdir"],
    "outdir": f"{outdir}/fiberseq/4_extract",
    "logdir": f"{logdir}/sample",
    "samples": samples,
    "Procedure": {
        "fibertools": config.get("Procedure", {}).get("fibertools") or "ft"
    },
    "genome": {
        "fasta": config.get("genome", {}).get("fasta")
    }
}
module fibertools_extract:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_extract_config
logger.info(f"Fiber-seq extract config: {fibertools_extract_config}")
use rule ft_extract from fibertools_extract as Fiberseq_ft_extract
