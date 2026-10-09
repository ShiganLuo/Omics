
include: "../common/common.smk"


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

# Gene type annotation references
te_bed = config.get("te_bed", "")
gene_tsv = config.get("gene_tsv", "")
# Species for tissue-specific annotation
species = config.get("species", "")
# ---------------------------------------------------------------------------
# Per-sample QC
# ---------------------------------------------------------------------------
rule scanpy_qc:
    """Perform QC on single-cell data using Scanpy."""
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
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_qc", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy QC for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_qc_{wildcards.counter}_{wildcards.sample_id}_{current_time}.sh")
            cmd = [params.python, params.script, "qc",
                   "-i", input.h5ad,
                   "-o", output.h5ad,
                   "-m", output.metrics,
                   "-n", str(params.min_genes),
                   "-N", str(params.max_genes),
                   "-t", str(params.max_pct_mt)]
            cmd += ["-p", str(output.plot_dir)]
            if params.use_mad:
                cmd += ["-M"]
            if params.scrublet:
                cmd += ["-s", "-d", str(params.doublet_rate)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy QC for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy QC for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during scanpy QC for sample {wildcards.sample_id}: {e}")

# ---------------------------------------------------------------------------
# Merge by tissue
# ---------------------------------------------------------------------------
def get_tissue_qc_files(wildcards):
    """Get QC'd h5ad paths for all samples in a tissue group."""
    sample_ids = tissue_samples[wildcards.tissue]
    return [outdir + f"/{sample_id}/{sample_id}_{wildcards.counter}_qc.h5ad" for sample_id in sample_ids]


rule scanpy_merge:
    """Merge QC'd single-cell data for a tissue group using Scanpy."""
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
        script=script,
        te_bed=te_bed,
        gene_tsv=gene_tsv
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_merge", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy merge for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_merge_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "merge",
                   "-i"] + list(input.h5ad) + ["-o", output.h5ad]
            cmd += ["-p", str(output.plot_dir)]
            if params.te_bed and params.gene_tsv:
                cmd += ["-b", params.te_bed, "-G", params.gene_tsv]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy merge for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy merge for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy merge for tissue {wildcards.tissue}: {e}")


rule scanpy_auto:
    """Automated Scanpy analysis for merged tissue data, including clustering and annotation.
    """
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_merged.h5ad"
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad",
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_auto_{counter}.log"
    threads: 4
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script,
        tissue=lambda wildcards: wildcards.tissue,
        llm_method=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("llm_method", ""),
        llm_model=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("llm_model", ""),
        llm_api_key=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("llm_api_key", ""),
        llm_base_url=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("llm_base_url", ""),
        resolution=lambda wildcards: params.get(wildcards.counter, {}).get("cluster", {}).get("resolution", 0.8),
        max_iterations=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("max_iterations", 3),
        min_genes=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("min_genes", 800),
        min_counts=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("min_counts", 3000),
        max_pct_mt=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("max_pct_mt", 20),
        n_pcs=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("n_pcs", 50),
        n_neighbors=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("n_neighbors", 50),
        n_top_genes=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("n_top_genes", 3000),
        batch_method=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("batch_method", "harmony"),
        batch_key=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("batch_key", "sample_id"),
        auto_n_pcs=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("auto_n_pcs", True),
        skip_te=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("skip_te", True),
        save_iteration_h5ad=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("save_iteration_h5ad", False),
        species=species
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_auto", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy auto for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script_path = os.path.join(sample_outdir, f"scanpy_auto_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "auto",
                   "-i", input.h5ad,
                   "-o", output.h5ad,
                   "-T", params.tissue,
                   "-r", str(params.resolution),
                   "-I", str(params.max_iterations),
                   "-n", str(params.min_genes),
                   "-U", str(params.min_counts),
                   "-t", str(params.max_pct_mt),
                   "-c", str(params.n_pcs),
                   "-k", str(params.n_neighbors),
                   "-g", str(params.n_top_genes),
                   "-B", params.batch_method,
                   "-K", params.batch_key]
            if params.species:
                cmd += ["-S", params.species]
            if params.llm_method:
                cmd += ["-L", params.llm_method]
            if params.llm_model:
                cmd += ["-l", params.llm_model]
            if params.llm_api_key:
                cmd += ["-Q", params.llm_api_key]
            if params.llm_base_url:
                cmd += ["-u", params.llm_base_url]
            if params.auto_n_pcs:
                cmd.append("-A")
            if params.skip_te:
                cmd.append("-x")
            if params.save_iteration_h5ad:
                cmd.append("-H")
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy auto for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy auto for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy auto for tissue {wildcards.tissue}: {e}")



rule scanpy_advanced:
    """Advanced Scanpy analysis for merged tissue data, including trajectory inference, RNA velocity, cell-cell communication, and CNV analysis.
    """
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad"
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
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_advanced", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy advanced analysis for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_advanced_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "advanced",
                    "-i", input.h5ad,
                    "-o", output.h5ad,
                    "-c", str(params.n_pcs),
                    "-k", str(params.n_neighbors)]
            if params.trajectory:
                cmd.append("-T")
            if params.velocity:
                cmd.append("-V")
            if params.communication:
                cmd.append("-L")
            if params.cnv:
                cmd.append("-C")
            if params.gtf:
                cmd += ["-G", params.gtf]
            if params.cnv_reference:
                cmd += ["-R", params.cnv_reference]
            cmd += ["-p", str(output.plot_dir)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy advanced analysis for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy advanced analysis for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy advanced analysis for tissue {wildcards.tissue}: {e}")



rule scanpy_differential_expression:
    """Perform differential expression analysis on merged tissue data using Scanpy.
    """
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
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_differential_expression", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy differential expression for tissue {wildcards.tissue} at {current_time}")
            sample_outdir = os.path.dirname(str(output.h5ad))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(str(output.plot_dir), exist_ok=True)
            script = os.path.join(sample_outdir, f"scanpy_de_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            cmd = [params.python, params.script, "de",
                   "-i", input.h5ad,
                   "-o", output.h5ad,
                   "-D", output.table]
            cmd += ["-p", str(output.plot_dir)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy differential expression for {wildcards.tissue} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy differential expression for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy differential expression for tissue {wildcards.tissue}: {e}")


