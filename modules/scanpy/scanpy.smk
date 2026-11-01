
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
    threads: rule_res("scanpy_qc").get("threads", 4)
    resources:
        mem_mb=rule_res("scanpy_qc").get("mem_mb", 16000),
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
            script = os.path.join(RECEIPT_DIR, f"scanpy_qc_{wildcards.counter}_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"scanpy_qc_{wildcards.counter}_{wildcards.sample_id}.sh")
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
                # Cap numeric-library threads to the rule's thread budget
                # (config resources -> rule_res -> threads) so the process
                # honors what snakemake allocated instead of grabbing every
                # core; the value is injected at script generation, never
                # hardcoded here.
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy QC for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
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
    threads: rule_res("scanpy_merge").get("threads", 2)
    resources:
        mem_mb=rule_res("scanpy_merge").get("mem_mb", 32000),
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
            script = os.path.join(RECEIPT_DIR, f"scanpy_merge_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"scanpy_merge_{wildcards.counter}_{wildcards.tissue}.sh")
            cmd = [params.python, params.script, "merge",
                   "-i"] + list(input.h5ad) + ["-o", output.h5ad]
            cmd += ["-p", str(output.plot_dir)]
            if params.te_bed and params.gene_tsv:
                cmd += ["-b", params.te_bed, "-G", params.gene_tsv]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                # Cap numeric-library threads to the rule's thread budget
                # (config resources -> rule_res -> threads) so the process
                # honors what snakemake allocated instead of grabbing every
                # core; the value is injected at script generation, never
                # hardcoded here.
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy merge for {wildcards.tissue} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
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
    threads: rule_res("scanpy_auto").get("threads", 4)
    resources:
        mem_mb=rule_res("scanpy_auto").get("mem_mb", 64000),
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
        regress_cell_cycle=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("regress_cell_cycle", True),
        save_iteration_h5ad=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("save_iteration_h5ad", False),
        marker_supplement=lambda wildcards: (params.get(wildcards.counter, {}).get("auto", {}).get("marker_supplement", {}) or {}).get(wildcards.tissue, ""),
        marker_supplement_llm=lambda wildcards: bool(params.get(wildcards.counter, {}).get("auto", {}).get("marker_supplement_llm", True)),
        exclude_genes=lambda wildcards: params.get(wildcards.counter, {}).get("auto", {}).get("exclude_genes", "") or "",
        exclude_from_deg=lambda wildcards: bool(params.get(wildcards.counter, {}).get("auto", {}).get("exclude_from_deg", False)),
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
            script_path = os.path.join(RECEIPT_DIR, f"scanpy_auto_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"scanpy_auto_{wildcards.counter}_{wildcards.tissue}.sh")
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
            if not params.regress_cell_cycle:
                cmd.append("-C")
            if params.marker_supplement:
                cmd += ["-X", str(params.marker_supplement)]
            if not params.marker_supplement_llm:
                cmd.append("-M")
            if params.exclude_genes:
                cmd += ["-E", str(params.exclude_genes)]
            if params.exclude_from_deg:
                cmd.append("-D")
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                # Cap numeric-library threads to the rule's thread budget
                # (config resources -> rule_res -> threads) so the process
                # honors what snakemake allocated instead of grabbing every
                # core; the value is injected at script generation, never
                # hardcoded here.
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy auto for {wildcards.tissue} at {current_time} completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy auto for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy auto for tissue {wildcards.tissue}: {e}")

def _reannotate_override(wildcards):
    """Override JSON path for this tissue/counter, or [] when disabled."""
    return (params.get(wildcards.counter, {}).get("auto", {})
            .get("cluster_annotations", {}) or {}
            ).get(wildcards.tissue, "") or []


def _reannotate_marker(wildcards):
    """Marker file input enforcing auto -> reannotate -> consumers ordering.

    The reannotate step overwrites the auto h5ad in place; downstream rules
    (advanced, pseudobulk) depend on this marker so they never read a
    pre-override auto.h5ad. Empty when no overrides are configured."""
    return (outdir_combine + f"/{wildcards.tissue}/{wildcards.tissue}_{wildcards.counter}_auto.h5ad.reannotated"
            if _reannotate_override(wildcards) else [])


rule scanpy_reannotate:
    """Apply manual cluster annotation overrides to a finished auto h5ad.

    IN-PLACE: backs the auto h5ad up to <auto>.bak, overwrites it with the
    corrected labels, and regenerates annotation_report.tsv /
    annotation_audit.json / marker plots inside the auto run's own report
    dir — no parallel _reannotated artifact set. The rule's declared
    output is a .reannotated marker file (DAG ordering only). No
    re-clustering, no LLM. Only enters the DAG for tissues with an
    override JSON configured.
    """
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad",
        overrides = _reannotate_override
    output:
        marker = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad.reannotated",
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_reannotate_{counter}.log"
    threads: rule_res("scanpy_reannotate").get("threads", 2)
    resources:
        mem_mb=rule_res("scanpy_reannotate").get("mem_mb", 16000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        python=python,
        script=script,
        marker_supplement=lambda wildcards: (params.get(wildcards.counter, {}).get("auto", {}).get("marker_supplement", {}) or {}).get(wildcards.tissue, ""),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_reannotate", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy reannotate (in-place) for tissue {wildcards.tissue}")
            sample_outdir = os.path.dirname(str(input.h5ad))
            script = os.path.join(RECEIPT_DIR, f"scanpy_reannotate_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"scanpy_reannotate_{wildcards.counter}_{wildcards.tissue}.sh")
            cmd = [params.python, params.script, "reannotate",
                   "-i", input.h5ad,
                   "-o", input.h5ad,
                   "-Y", input.overrides,
                   "-T", wildcards.tissue]
            if params.marker_supplement:
                cmd += ["-X", str(params.marker_supplement)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy reannotate for {wildcards.tissue} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            with open(str(output.marker), "w") as f:
                f.write(f"reannotated at {current_time}\n")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy reannotate for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy reannotate for tissue {wildcards.tissue}: {e}")

rule scanpy_advanced:
    """Advanced Scanpy analysis for merged tissue data, including trajectory inference, RNA velocity, cell-cell communication, and CNV analysis.
    """
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad",
        # ordering guard: when overrides are configured for this tissue,
        # reannotate rewrites the auto h5ad in place; this marker ensures
        # advanced reads the post-override labels.
        reannotated = _reannotate_marker
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_advanced.h5ad",
        plot_dir = directory(outdir_combine + "/{tissue}/plots/{counter}/advanced")
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_advanced_{counter}.log"
    threads: rule_res("scanpy_advanced").get("threads", 4)
    resources:
        mem_mb=rule_res("scanpy_advanced").get("mem_mb", 64000),
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
            script = os.path.join(RECEIPT_DIR, f"scanpy_advanced_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"scanpy_advanced_{wildcards.counter}_{wildcards.tissue}.sh")
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
                # Cap numeric-library threads to the rule's thread budget
                # (config resources -> rule_res -> threads) so the process
                # honors what snakemake allocated instead of grabbing every
                # core; the value is injected at script generation, never
                # hardcoded here.
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy advanced analysis for {wildcards.tissue} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy advanced analysis for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy advanced analysis for tissue {wildcards.tissue}: {e}")

        # ---------------------------------------------------------------------------
        # Pseudo-bulk DEG (PyDESeq2) — modules/scanpy/bin/pseudobulk.py
        #   Two-file design: --h5ad (scanpy_auto output) supplies the cell_type
        #   annotations; --counts-h5ad (scanpy_merge output) supplies the
        #   full-gene integer layers["counts"] (the auto file's counts layer is
        #   HVG-subset and its raw slot is log1p-normalised — neither is valid
        #   pseudobulk input). Runs PyDESeq2 per
        #   Params.scanpy.pseudobulk.compare entry, plots volcano/heatmap/GO.
        #   Enabled per counter via
        #   Params.scanpy.<counter>.pseudobulk.{enabled,compare,...}
        # ---------------------------------------------------------------------------
rule scanpy_pseudobulk:
    """Pseudo-bulk DEG via PyDESeq2: annotations from scanpy_auto, full-gene
    raw counts from the merged h5ad. When overrides are configured for the
    tissue, the reannotate marker orders pb after the in-place label fix."""
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad",
        reannotated = _reannotate_marker,
        counts = outdir_combine + "/{tissue}/{tissue}_{counter}_merged.h5ad"
    output:
        report = outdir_combine + "/{tissue}/pseudobulk/{counter}/pseudobulk_report.html",
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_pseudobulk_{counter}.log"
    threads: rule_res("scanpy_pseudobulk").get("threads", 4)
    resources:
        mem_mb=rule_res("scanpy_pseudobulk").get("mem_mb", 32000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml")
    params:
        pb_cfg = lambda wildcards: params.get(wildcards.counter, {}).get("pseudobulk", {}),
        python = python,
        script = os.path.join(ROOT_DIR, "modules", "scanpy", "bin", "pseudobulk.py"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_pseudobulk", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy pseudobulk for tissue {wildcards.tissue} at {current_time}")
            pb = params.pb_cfg
            out_root = os.path.dirname(str(output.report))
            os.makedirs(out_root, exist_ok=True)
            # dual-write receipts: timestamped copy in log/cmd, plain copy in outdir
            script = os.path.join(RECEIPT_DIR, f"scanpy_pseudobulk_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            plain_script = os.path.join(out_root, f"scanpy_pseudobulk_{wildcards.counter}_{wildcards.tissue}.sh")
            cmd = [params.python, params.script,
                   "--h5ad", f"{wildcards.tissue}={input.h5ad}",
                   "--counts-h5ad", f"{wildcards.tissue}={input.counts}",
                   "--out-dir", out_root,
                   "--min-cells", str(pb.get("min_cells", 30)),
                   "--min-counts", str(pb.get("min_counts", 1000)),
                   "--n-pseudo-reps", str(pb.get("n_pseudo_reps", 3))]
            # entries: tissue:ref_sample:treat_sample[=label] — every entry
            # is scoped to ONE tissue (pseudobulk.py rejects a compare whose
            # tissue has no --h5ad definition), so this per-tissue job only
            # takes the entries prefixed with its own wildcard.
            tissue_compares = [
                str(c) for c in (pb.get("compare") or [])
                if str(c).split(":", 1)[0] == wildcards.tissue
            ]
            if not tissue_compares:
                raise RuntimeError(
                    f"Params.scanpy.{wildcards.counter}.pseudobulk.compare "
                    f"has no entry for tissue '{wildcards.tissue}' — expected "
                    f"entries like '{wildcards.tissue}:ref_sample:"
                    f"treat_sample[=label]'")
            for cmp_ in tissue_compares:
                cmd += ["--compare", cmp_]
            for ct in pb.get("top_cell_types", []) or []:
                cmd += ["--top-cell-types", str(ct)]
            for g in pb.get("highlight_genes", []) or []:
                cmd += ["--highlight-genes", str(g)]
            # --deg-mode is nargs='+': ONE flag carrying all values —
            # repeated flags overwrite each other in argparse and only the
            # last mode would run.
            cmd += ["--deg-mode"] + [str(m) for m in (pb.get("deg_mode") or ["all"])]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scanpy pseudobulk for {wildcards.tissue} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scanpy pseudobulk for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scanpy pseudobulk for tissue {wildcards.tissue}: {e}")
