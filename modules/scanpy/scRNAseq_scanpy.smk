
include: "../common/common.smk"
import os
import time
import shlex
from pathlib import Path
from snakemake.logging import logger

ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
input_h5ad = config.get("input_h5ad", "")
sample_h5ad = config.get("sample_h5ad", {})
inputs = list(sample_h5ad.values()) if sample_h5ad else [input_h5ad]
if not inputs or not inputs[0]:
    raise ValueError("scRNAseq requires input_h5ad or sample_h5ad")
procedure = config.get("Procedure", {})
params = config.get("Params", {}).get("scanpy", {})
advanced = params.get("advanced", {})
python = procedure.get("python") or "python"
script = os.path.join(ROOT_DIR, "modules", "scanpy", "bin", "scRNAseq.py")

qc_h5ad = outdir + "/scanpy/qc/filtered.h5ad"
cluster_h5ad = outdir + "/scanpy/cluster/clustered.h5ad"
advanced_h5ad = outdir + "/scanpy/advanced/advanced.h5ad"
de_h5ad = outdir + "/scanpy/de/differential_expression.h5ad"
de_input = advanced_h5ad if advanced.get("trajectory", True) or advanced.get("velocity", False) else cluster_h5ad

rule scanpy_qc:
    input:
        h5ad=inputs
    output:
        h5ad=qc_h5ad,
        metrics=outdir + "/scanpy/qc/qc_metrics.tsv"
    log:
        logdir + "/scanpy/scanpy_qc.log"
    threads: 4
    conda:
        "scRNAseq_scanpy.yaml"
    container:
        sif("scRNAseq_scanpy.yaml")
    params:
        python=python,
        script=script,
        min_genes=params.get("min_genes", 200),
        max_genes=params.get("max_genes", 6000),
        max_pct_mt=params.get("max_pct_mt", 20),
        n_top_genes=params.get("n_top_genes", 3000),
        batch_key=params.get("batch_key", "")
    run:
        log_path = str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(output.h5ad), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            command_script = os.path.join(os.path.dirname(output.h5ad), f"scanpy_qc_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "qc", "--input"] + list(input.h5ad) + ["--output", output.h5ad, "--metrics", output.metrics, "--min-genes", str(params.min_genes), "--max-genes", str(params.max_genes), "--max-pct-mt", str(params.max_pct_mt), "--n-top-genes", str(params.n_top_genes)]
            if params.batch_key:
                cmd += ["--batch-key", params.batch_key]
            with open(command_script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
            shell(f"bash {shlex.quote(command_script)} >> {shlex.quote(log_path)} 2>&1")
        except Exception as exc:
            with open(log_path, "a") as handle:
                handle.write(f"scanpy QC failed: {exc}\n")
            raise

rule scanpy_cluster:
    input:
        h5ad=qc_h5ad
    output:
        h5ad=cluster_h5ad,
        markers=outdir + "/scanpy/cluster/markers.tsv"
    log:
        logdir + "/scanpy/scanpy_cluster.log"
    threads: 4
    conda:
        "scRNAseq_scanpy.yaml"
    container:
        sif("scRNAseq_scanpy.yaml")
    params:
        python=python,
        script=script,
        n_pcs=params.get("n_pcs", 50),
        n_neighbors=params.get("n_neighbors", 15),
        resolution=params.get("resolution", 0.8)
    run:
        log_path = str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(output.h5ad), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            command_script = os.path.join(os.path.dirname(output.h5ad), f"scanpy_cluster_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "cluster", "--input", input.h5ad, "--output", output.h5ad, "--markers", output.markers, "--n-pcs", str(params.n_pcs), "--n-neighbors", str(params.n_neighbors), "--resolution", str(params.resolution)]
            with open(command_script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
            shell(f"bash {shlex.quote(command_script)} >> {shlex.quote(log_path)} 2>&1")
        except Exception as exc:
            with open(log_path, "a") as handle:
                handle.write(f"scanpy clustering failed: {exc}\n")
            raise

if advanced.get("trajectory", True) or advanced.get("velocity", False):
    rule scanpy_advanced:
        input:
            h5ad=cluster_h5ad
        output:
            h5ad=advanced_h5ad
        log:
            logdir + "/scanpy/scanpy_advanced.log"
        threads: 4
        conda:
            "scRNAseq_scanpy.yaml"
        container:
            sif("scRNAseq_scanpy.yaml")
        params:
            python=python,
            script=script,
            n_pcs=params.get("n_pcs", 50),
            n_neighbors=params.get("n_neighbors", 15),
            trajectory=advanced.get("trajectory", True),
            velocity=advanced.get("velocity", False),
            communication=advanced.get("communication", False),
            cnv=advanced.get("cnv", False)
        run:
            log_path = str(log)
            try:
                os.makedirs(os.path.dirname(log_path), exist_ok=True)
                os.makedirs(os.path.dirname(output.h5ad), exist_ok=True)
                current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
                command_script = os.path.join(os.path.dirname(output.h5ad), f"scanpy_advanced_{current_time}.sh")
                cmd = [params.python, params.script, "--mode", "advanced", "--input", input.h5ad, "--output", output.h5ad, "--n-pcs", str(params.n_pcs), "--n-neighbors", str(params.n_neighbors)]
                if params.trajectory:
                    cmd.append("--trajectory")
                if params.velocity:
                    cmd.append("--velocity")
                if params.communication:
                    cmd.append("--communication")
                if params.cnv:
                    cmd.append("--cnv")
                with open(command_script, "w") as handle:
                    handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                    handle.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                shell(f"bash {shlex.quote(command_script)} >> {shlex.quote(log_path)} 2>&1")
            except Exception as exc:
                with open(log_path, "a") as handle:
                    handle.write(f"scanpy advanced analysis failed: {exc}\n")
                raise

    rule scanpy_differential_expression:
        input:
            h5ad=de_input
        output:
            h5ad=de_h5ad,
            table=outdir + "/scanpy/de/markers.tsv"
        log:
            logdir + "/scanpy/scanpy_differential_expression.log"
        threads: 2
        conda:
            "scRNAseq_scanpy.yaml"
        container:
            sif("scRNAseq_scanpy.yaml")
        params:
            python=python,
            script=script
        run:
            log_path = str(log)
            try:
                os.makedirs(os.path.dirname(log_path), exist_ok=True)
                os.makedirs(os.path.dirname(output.h5ad), exist_ok=True)
                current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
                command_script = os.path.join(os.path.dirname(output.h5ad), f"scanpy_de_{current_time}.sh")
                cmd = [params.python, params.script, "--mode", "de", "--input", input.h5ad, "--output", output.h5ad, "--deg", output.table]
                with open(command_script, "w") as handle:
                    handle.write("#!/usr/bin/env bash\\nset -euo pipefail\\n")
                    handle.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\\n")
                shell(f"bash {shlex.quote(command_script)} >> {shlex.quote(log_path)} 2>&1")
            except Exception as exc:
                with open(log_path, "a") as handle:
                    handle.write(f"scanpy differential expression failed: {exc}\\n")
                raise
