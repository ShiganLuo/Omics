"""Fiber-seq subworkflow for single-molecule chromatin accessibility analysis.

Workflow (Revio SPRQ — m6A calls already in CCS BAM from jasmine):
  1. Add nucleosome calls (ft add-nucleosomes)
  2. Call FIREs (ft fire)
  3. Extract data to BED format (ft extract --all)
  4. Call FIRE peaks (ft call-peaks, FDR-based)
  5. Collect QC metrics (ft qc)

For older instruments (Sequel II/IIe) without jasmine m6A calls,
a separate ft predict-m6a step is needed before step 1.

Reference:
  - Stergachis et al., 2020, Science (Fiber-seq original paper, DOI: 10.1126/science.aaz1646)
  - Jha, Bohaczuk et al., 2024, Genome Research (fibertools-rs, DOI: 10.1101/gr.279095.124)
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

rule all:
    input:
        outfiles


# ============================================================
# Step 1: Add nucleosome calls (input: raw CCS BAM with m6A from jasmine)
# ============================================================
fibertools_nuc_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": indir,
    "outdir": f"{outdir}/fiberseq/1_nucleosomes",
    "upstream_outdir": indir,
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
# Step 2: FIRE calling
# ============================================================
fibertools_fire_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fibertools_nuc_config["outdir"],
    "outdir": f"{outdir}/fiberseq/2_fire",
    "upstream_outdir": fibertools_nuc_config["outdir"],
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
# Step 3: Extract data to BED format
# ============================================================
fibertools_extract_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fibertools_fire_config["outdir"],
    "outdir": f"{outdir}/fiberseq/3_extract",
    "upstream_outdir": fibertools_fire_config["outdir"],
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


# ============================================================
# Step 4: FIRE peak calling
# ============================================================
fibertools_peaks_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fibertools_fire_config["outdir"],
    "outdir": f"{outdir}/fiberseq/4_peaks",
    "upstream_outdir": fibertools_fire_config["outdir"],
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
module fibertools_peaks:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_peaks_config
logger.info(f"Fiber-seq peaks config: {fibertools_peaks_config}")
use rule ft_call_peaks from fibertools_peaks as Fiberseq_ft_call_peaks


# ============================================================
# Step 5: QC metrics
# ============================================================
fibertools_qc_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "indir": fibertools_fire_config["outdir"],
    "outdir": f"{outdir}/fiberseq/5_qc",
    "upstream_outdir": fibertools_fire_config["outdir"],
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
module fibertools_qc:
    snakefile: "../modules/fibertools/fibertools.smk"
    config: fibertools_qc_config
logger.info(f"Fiber-seq QC config: {fibertools_qc_config}")
use rule ft_qc from fibertools_qc as Fiberseq_ft_qc
