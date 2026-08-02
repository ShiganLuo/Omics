"""Orchestrate atomic Scanpy modules for standardized scRNA-seq analysis."""

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
input_h5ad = config.get("input_h5ad", "")
sample_h5ad = config.get("sample_h5ad", {})
analysis = config.get("Params", {}).get("scanpy", {})
advanced = analysis.get("advanced", {})

if not input_h5ad and not sample_h5ad:
    raise ValueError("scRNAseq requires input_h5ad or sample_h5ad")

qc_h5ad = f"{outdir}/scanpy/qc/filtered.h5ad"
cluster_h5ad = f"{outdir}/scanpy/cluster/clustered.h5ad"
advanced_enabled = advanced.get("trajectory", True) or advanced.get("velocity", False) or advanced.get("communication", False) or advanced.get("cnv", False)
advanced_h5ad = f"{outdir}/scanpy/advanced/advanced.h5ad"
de_h5ad = f"{outdir}/scanpy/de/differential_expression.h5ad"
outfiles = config.get("outfiles") or [
    qc_h5ad,
    f"{outdir}/scanpy/qc/qc_metrics.tsv",
    cluster_h5ad,
    f"{outdir}/scanpy/cluster/markers.tsv",
    *( [advanced_h5ad] if advanced_enabled else [] ),
    de_h5ad,
    f"{outdir}/scanpy/de/markers.tsv",
]

scanpy_config = dict(config)
scanpy_config.update({
    "ROOT_DIR": ROOT_DIR,
    "indir": indir,
    "outdir": outdir,
    "logdir": logdir,
    "input_h5ad": input_h5ad,
    "sample_h5ad": sample_h5ad,
    "outfiles": outfiles,
})

module scanpy:
    snakefile: "../modules/scanpy/scRNAseq_scanpy.smk"
    config: scanpy_config

use rule scanpy_qc from scanpy as scRNAseq_scanpy_qc
use rule scanpy_cluster from scanpy as scRNAseq_scanpy_cluster
if advanced_enabled:
    use rule scanpy_advanced from scanpy as scRNAseq_scanpy_advanced
use rule scanpy_differential_expression from scanpy as scRNAseq_scanpy_de

rule all:
    input:
        outfiles
