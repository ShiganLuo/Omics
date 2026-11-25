
"""Orchestrate atomic spatial transcriptomics modules."""
ROOT_DIR = config.get("ROOT_DIR", ".")
_wf = config.get("Params", {}).get("workflow", {})
outdir = _wf.get("outdir") or "output"
logdir = _wf.get("logdir") or "logs"
input_h5ad = _wf.get("input_h5ad") or ""
visium_h5 = _wf.get("visium_h5") or ""
spatial_dir = _wf.get("spatial_dir") or ""
advanced = config.get("Params", {}).get("spatial", {}).get("spatial_scanpy_advanced", {})
if not input_h5ad and not visium_h5:
    raise ValueError("spatial transcriptomics requires input_h5ad or visium_h5")
advanced_enabled = advanced.get("spatial_autocorrelation", True)
qc = f"{outdir}/spatial/qc/filtered.h5ad"
cluster = f"{outdir}/spatial/cluster/clustered.h5ad"
advanced_output = f"{outdir}/spatial/advanced/spatial_autocorrelation.h5ad"
outfiles = config.get("outfiles") or [qc, f"{outdir}/spatial/qc/qc_metrics.tsv", cluster, f"{outdir}/spatial/cluster/spatial_markers.tsv"]
if advanced_enabled:
    outfiles.append(advanced_output)
spatial_config = dict(config)
spatial_config.update({
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "outdir": outdir,
    "logdir": logdir,
    "input_h5ad": input_h5ad,
    "visium_h5": visium_h5,
    "spatial_dir": spatial_dir,
    "outfiles": outfiles,
})
module spatial_scanpy:
    snakefile: "../modules/spatial_scanpy/spatial_scanpy.smk"
    config: spatial_config
use rule spatial_scanpy_qc from spatial_scanpy as spatial_transcriptomics_qc
use rule spatial_scanpy_cluster from spatial_scanpy as spatial_transcriptomics_cluster
if advanced_enabled:
    use rule spatial_scanpy_advanced from spatial_scanpy as spatial_transcriptomics_advanced
rule all:
    input:
        outfiles
