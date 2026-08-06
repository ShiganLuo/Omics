
include: "../common/common.smk"
import os
import shlex
import time

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
input_h5ad = config.get("input_h5ad", "")
visium_h5 = config.get("visium_h5", "")
spatial_dir = config.get("spatial_dir", "")
params = config.get("Params", {}).get("spatial", {})
advanced = params.get("advanced", {})
python = config.get("Procedure", {}).get("python") or "python"
script = os.path.join(ROOT_DIR, "modules", "spatial_scanpy", "bin", "spatial_transcriptomics.py")
qc_h5ad = outdir + "/spatial/qc/filtered.h5ad"
cluster_h5ad = outdir + "/spatial/cluster/clustered.h5ad"
advanced_h5ad = outdir + "/spatial/advanced/spatial_autocorrelation.h5ad"

if not input_h5ad and not visium_h5:
    raise ValueError("spatial transcriptomics requires input_h5ad or visium_h5")

def _common_command(mode, output, extra):
    cmd = [python, script, "--mode", mode]
    if input_h5ad:
        cmd += ["--input-h5ad", input_h5ad]
    else:
        cmd += ["--visium-h5", visium_h5]
        if spatial_dir:
            cmd += ["--spatial-dir", spatial_dir]
    cmd += ["--output", output]
    cmd += extra
    return cmd

rule spatial_scanpy_qc:
    input:
        h5ad=input_h5ad if input_h5ad else visium_h5,
        spatial=lambda wc: spatial_dir if spatial_dir else []
    output:
        h5ad=qc_h5ad,
        metrics=outdir + "/spatial/qc/qc_metrics.tsv"
    log:
        logdir + "/spatial/spatial_scanpy_qc.log"
    threads: 4
    conda:
        "spatial_scanpy.yaml"
    container:
        sif("spatial_scanpy.yaml")
    params:
        command=lambda wc: _common_command("qc", qc_h5ad, ["--metrics", outdir + "/spatial/qc/qc_metrics.tsv", "--min-genes", str(params.get("min_genes", 100)), "--max-pct-mt", str(params.get("max_pct_mt", 25)), "--n-top-genes", str(params.get("n_top_genes", 3000))])
    run:
        log_path=str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True); os.makedirs(os.path.dirname(output.h5ad), exist_ok=True)
            stamp=time.strftime("%Y%m%d_%H%M%S", time.localtime()); command_script=os.path.join(os.path.dirname(output.h5ad), f"spatial_qc_{stamp}.sh")
            with open(command_script,"w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(shlex.quote(str(x)) for x in params.command)+"\n")
            shell(f"bash {shlex.quote(command_script)} >> {shlex.quote(log_path)} 2>&1")
        except Exception as exc:
            with open(log_path,"a") as handle: handle.write(f"spatial QC failed: {exc}\n")
            raise

rule spatial_scanpy_cluster:
    input: h5ad=qc_h5ad
    output: h5ad=cluster_h5ad, genes=outdir + "/spatial/cluster/spatial_markers.tsv"
    log: logdir + "/spatial/spatial_scanpy_cluster.log"
    threads: 4
    conda: "spatial_scanpy.yaml"
    container:
        sif("spatial_scanpy.yaml")
    params:
        command=lambda wc: _common_command("cluster", cluster_h5ad, ["--genes", outdir + "/spatial/cluster/spatial_markers.tsv", "--n-pcs", str(params.get("n_pcs", 50)), "--n-neighbors", str(params.get("n_neighbors", 15)), "--resolution", str(params.get("resolution", 0.8))])
    run:
        log_path=str(log)
        try:
            os.makedirs(os.path.dirname(log_path),exist_ok=True);os.makedirs(os.path.dirname(output.h5ad),exist_ok=True);stamp=time.strftime("%Y%m%d_%H%M%S",time.localtime());command_script=os.path.join(os.path.dirname(output.h5ad),f"spatial_cluster_{stamp}.sh")
            with open(command_script,"w") as handle: handle.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(shlex.quote(str(x)) for x in params.command)+"\n")
            shell(f"bash {shlex.quote(command_script)} >> {shlex.quote(log_path)} 2>&1")
        except Exception as exc:
            with open(log_path,"a") as handle: handle.write(f"spatial clustering failed: {exc}\n")
            raise

if advanced.get("spatial_autocorrelation", True):
    rule spatial_scanpy_advanced:
        input: h5ad=cluster_h5ad
        output: h5ad=advanced_h5ad
        log: logdir + "/spatial/spatial_scanpy_advanced.log"
        threads: 4
        conda: "spatial_scanpy.yaml"
        container:
            sif("spatial_scanpy.yaml")
        params: command=lambda wc: _common_command("advanced", advanced_h5ad, [])
        run:
            log_path=str(log)
            try:
                os.makedirs(os.path.dirname(log_path),exist_ok=True);os.makedirs(os.path.dirname(output.h5ad),exist_ok=True);stamp=time.strftime("%Y%m%d_%H%M%S",time.localtime());command_script=os.path.join(os.path.dirname(output.h5ad),f"spatial_advanced_{stamp}.sh")
                with open(command_script,"w") as handle: handle.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(shlex.quote(str(x)) for x in params.command)+"\n")
                shell(f"bash {shlex.quote(command_script)} >> {shlex.quote(log_path)} 2>&1")
            except Exception as exc:
                with open(log_path,"a") as handle: handle.write(f"spatial advanced analysis failed: {exc}\n")
                raise
