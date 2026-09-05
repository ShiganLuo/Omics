
include: "../common/common.smk"
import os
import time
import shlex
from snakemake.logging import logger

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outdir_combine = config.get("outdir_combine", outdir)
logdir_combine = config.get("logdir_combine", logdir)


procedure = config.get("Procedure", {})
params = config.get("Params", {}).get("scanpy", {})
# tissue_samples: {tissue: [sample_id, ...]} — from node.py
tissue_samples = params.get("tissue_samples", {})
tissues = sorted(tissue_samples.keys())
python = procedure.get("python") or "python"
script = os.path.join(ROOT_DIR, "modules", "scanpy", "bin", "scRNAseq.py")
# ---------------------------------------------------------------------------
# Per-sample QC
# ---------------------------------------------------------------------------
rule scanpy_qc:
    input:
        h5ad = indir + "/{sample_id}/{sample_id}_{counter}.h5ad"
    output:
        h5ad = outdir + "/{sample_id}/{sample_id}_{counter}_qc.h5ad",
        metrics = outdir + "/{sample_id}/{sample_id}_{counter}_qc_metrics.tsv",
        plot_dir = directory(outdir + "/{sample_id}/plots/{counter}")
    log:
        logdir + "/{sample_id}/scanpy_qc_{counter}.log"
    threads: 4
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script,
        min_genes = lambda wildcards: params.get(wildcards.counter, {}).get("qc",{}).get("min_genes", 200),
        max_genes=lambda wildcards: params.get(wildcards.counter, {}).get("qc",{}).get("max_genes", 6000),
        max_pct_mt=lambda wildcards: params.get(wildcards.counter, {}).get("qc",{}).get("max_pct_mt", 20),
        n_top_genes=lambda wildcards: params.get(wildcards.counter, {}).get("qc",{}).get("n_top_genes", 3000),
        use_mad=lambda wildcards: params.get(wildcards.counter, {}).get("qc",{}).get("use_mad", False),
        scrublet=lambda wildcards: params.get(wildcards.counter, {}).get("qc",{}).get("scrublet", True),
        doublet_rate=lambda wildcards: params.get(wildcards.counter, {}).get("qc",{}).get("doublet_rate", 0.06)
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("scanpy_qc", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy QC for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_qc_{wildcards.counter}_{wildcards.sample_id}_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "qc",
                   "--input", input.h5ad,
                   "--output", output.h5ad,
                   "--metrics", output.metrics,
                   "--min-genes", str(params.min_genes),
                   "--max-genes", str(params.max_genes),
                   "--max-pct-mt", str(params.max_pct_mt),
                   "--n-top-genes", str(params.n_top_genes)]
            cmd += ["--plot-dir", str(output.plot_dir)]
            if params.use_mad:
                cmd += ["--use-mad"]
            if params.scrublet:
                cmd += ["--scrublet", "--doublet-rate", str(params.doublet_rate)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy QC for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during scanpy QC for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during scanpy QC for sample {wildcards.sample_id}: {e}")
            raise e

# ---------------------------------------------------------------------------
# Merge by tissue
# ---------------------------------------------------------------------------
def get_tissue_qc_files(wildcards):
    """Get QC'd h5ad paths for all samples in a tissue group."""
    sample_ids = tissue_samples[wildcards.tissue]
    return [outdir + f"/{sample_id}/{sample_id}_{wildcards.counter}_qc.h5ad" for sample_id in sample_ids]


rule scanpy_merge:
    input:
        h5ad = get_tissue_qc_files
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_merged.h5ad",
        plot_dir = directory(outdir_combine + "/{tissue}/plots/{counter}/merge")
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_merge_{counter}.log"
    threads: 2
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("scanpy_merge", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy merge for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_merge_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "merge",
                   "--input"] + list(input.h5ad) + ["--output", output.h5ad]
            cmd += ["--plot-dir", str(output.plot_dir)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy merge for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during scanpy merge for tissue {wildcards.tissue}: {e}\n")
            logger.error(f"Error occurred during scanpy merge for tissue {wildcards.tissue}: {e}")
            raise e



rule scanpy_cluster:
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_merged.h5ad"
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_clustered.h5ad",
        markers = outdir_combine + "/{tissue}/{tissue}_{counter}_markers.tsv",
        plot_dir = directory(outdir_combine + "/{tissue}/plots/{counter}/cluster")
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_cluster_{counter}.log"
    threads: 4
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script,
        n_pcs=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("n_pcs", 50),
        n_neighbors=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("n_neighbors", 50),
        resolution=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("resolution", 0.8),
        n_top_genes=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("n_top_genes", 3000),
        batch_method=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("batch_method", "harmony"),
        batch_key=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("batch_key", "sample_id"),
        auto_n_pcs=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("auto_n_pcs", False)
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("scanpy_cluster", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy cluster for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_cluster_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "cluster",
                   "--input", input.h5ad,
                   "--output", output.h5ad,
                   "--markers", output.markers,
                   "--n-pcs", str(params.n_pcs),
                   "--n-neighbors", str(params.n_neighbors),
                   "--resolution", str(params.resolution),
                   "--n-top-genes", str(params.n_top_genes),
                   "--batch-method", params.batch_method,
                   "--batch-key", params.batch_key]
            if params.auto_n_pcs:
                cmd.append("--auto-n-pcs")
            cmd += ["--plot-dir", str(output.plot_dir)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy cluster for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during scanpy clustering for tissue {wildcards.tissue}: {e}\n")
            logger.error(f"Error occurred during scanpy clustering for tissue {wildcards.tissue}: {e}")
            raise e


rule scanpy_annotate:
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_clustered.h5ad"
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_annotated.h5ad",
        plot_dir = directory(outdir_combine + "/{tissue}/plots/{counter}/annotate")
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_annotate_{counter}.log"
    threads: 4
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script,
        marker_file=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("marker_file", ""),
        celltypist_model=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("celltypist_model", ""),
        llm_method=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("llm_method", ""),
        llm_model=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("llm_model", ""),
        llm_api_key=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("llm_api_key", ""),
        llm_base_url=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("llm_base_url", ""),
        llm_top_genes=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("llm_top_genes", 30),
        annotate_group=lambda wildcards: params.get(wildcards.counter, {}).get("annotate", {}).get("annotate_group", ""),
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("scanpy_annotate", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy annotation for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_annotate_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "annotate",
                    "--input", input.h5ad,
                    "--output", output.h5ad,
                    "--tissue", wildcards.tissue]
            if params.marker_file:
                cmd += ["--marker-file", params.marker_file]
            if params.celltypist_model:
                cmd += ["--celltypist-model", params.celltypist_model]
            if params.llm_method:
                cmd += ["--llm-method", params.llm_method]
                if params.llm_model:
                    cmd += ["--llm-model", params.llm_model]
                if params.llm_api_key:
                    cmd += ["--llm-api-key", params.llm_api_key]
                if params.llm_base_url:
                    cmd += ["--llm-base-url", params.llm_base_url]
                cmd += ["--llm-top-genes", str(params.llm_top_genes)]
            if params.annotate_group:
                cmd += ["--annotate-group", params.annotate_group]
            cmd += ["--plot-dir", str(output.plot_dir)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy annotation for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during scanpy annotation for tissue {wildcards.tissue}: {e}\n")
            logger.error(f"Error occurred during scanpy annotation for tissue {wildcards.tissue}: {e}")
            raise e

rule scanpy_advanced:
    input:
        h5ad =  outdir_combine + "/{tissue}/{tissue}_{counter}_annotated.h5ad"
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_advanced.h5ad",
        plot_dir = directory(outdir_combine + "/{tissue}/plots/{counter}/advanced")
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_advanced_{counter}.log"
    threads: 4
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script,
        n_pcs=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("n_pcs", 50),
        n_neighbors=lambda wildcards: params.get(wildcards.counter,{}).get("cluster",{}).get("n_neighbors", 15),
        trajectory=lambda wildcards: params.get(wildcards.counter,{}).get("advanced",{}).get("trajectory", True),
        velocity=lambda wildcards: params.get(wildcards.counter,{}).get("advanced",{}).get("velocity", False),
        communication=lambda wildcards: params.get(wildcards.counter,{}).get("advanced",{}).get("communication", False),
        cnv=lambda wildcards: params.get(wildcards.counter,{}).get("advanced",{}).get("cnv", False),
        gtf=lambda wildcards: params.get(wildcards.counter,{}).get("advanced",{}).get("gtf", ""),
        cnv_reference=lambda wildcards: params.get(wildcards.counter,{}).get("advanced",{}).get("cnv_reference", "")
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("scanpy_advanced", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy advanced analysis for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_advanced_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "advanced",
                    "--input", input.h5ad,
                    "--output", output.h5ad,
                    "--n-pcs", str(params.n_pcs),
                    "--n-neighbors", str(params.n_neighbors)]
            if params.trajectory:
                cmd.append("--trajectory")
            if params.velocity:
                cmd.append("--velocity")
            if params.communication:
                cmd.append("--communication")
            if params.cnv:
                cmd.append("--cnv")
            if params.gtf:
                cmd += ["--gtf", params.gtf]
            if params.cnv_reference:
                cmd += ["--cnv-reference", params.cnv_reference]
            cmd += ["--plot-dir", str(output.plot_dir)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy advanced analysis for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during scanpy advanced analysis for tissue {wildcards.tissue}: {e}\n")
            logger.error(f"Error occurred during scanpy advanced analysis for tissue {wildcards.tissue}: {e}")
            raise e


rule scanpy_differential_expression:
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_advanced.h5ad"
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_de.h5ad",
        table = outdir_combine + "/{tissue}/{tissue}_{counter}_markers.tsv",
        plot_dir = directory(outdir_combine + "/{tissue}/plots/{counter}/de")
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_de_{counter}.log"
    threads: 2
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("scanpy_differential_expression", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy differential expression for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_de_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "--mode", "de",
                   "--input", input.h5ad,
                   "--output", output.h5ad,
                   "--deg", output.table]
            cmd += ["--plot-dir", str(output.plot_dir)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy differential expression for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during scanpy differential expression for tissue {wildcards.tissue}: {e}\n")
            logger.error(f"Error occurred during scanpy differential expression for tissue {wildcards.tissue}: {e}")
            raise e

# ---------------------------------------------------------------------------
# Result aggregation (for subworkflow use rule)
# ---------------------------------------------------------------------------
rule scanpy_result:
    input:
        h5ad = outdir_combine + "/{t}/{t}_{counter}_de.h5ad",
        table = outdir_combine + "/{t}/{t}_{counter}_markers.tsv"