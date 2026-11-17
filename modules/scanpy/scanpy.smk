
include: "../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
outdir_combine = config.get("outdir_combine", outdir)
# pseudobulk deliverables root (scRNAseq: <outdir>/diff/pseudobulk)
outdir_pseudobulk = config.get("outdir_pseudobulk", outdir_combine)
logdir_combine = config.get("logdir_combine", logdir)

procedure = config.get("Procedure", {})
params = config.get("Params", {}).get("scanpy", {})
# tissue_samples: {tissue: [sample_id, ...]} — from node.py
tissue_samples = params.get("tissue_samples", {})
tissues = sorted(tissue_samples.keys())
python = procedure.get("python") or "python"
script = os.path.join(ROOT_DIR, "modules", "scanpy", "bin", "scRNAseq.py")

# 参数精确到 rule:Params.scanpy.<counter>.<rule_name>.<key>(rule_param 见
# common.smk;无 counter 级共享回退)。原 qc/auto/advanced/cluster/pseudobulk
# 参数段逐 key 挂到消费它的 rule 名下(scanpy_qc / scanpy_merge / scanpy_auto /
# scanpy_reannotate / scanpy_advanced / scanpy_pseudobulk)各 rule 独立。


def _image_formats(wildcards, rule_name: str) -> list:
    """Resolve image formats for one rule's figures (modules.md #9).

    参数精确到 rule:Params.scanpy.<counter>.<rule_name>.formats -> ["png"]
    (today's fixed output names). No counter-level shared ``formats``
    fallback (rule_param, common.smk) — each rule carries its own key.
    Plot outputs are directory() (pre-existing shape, unchanged); the files
    written inside carry the requested extensions via the script's
    repeatable ``--format`` flag.
    """
    fmts = rule_param(params.get(wildcards.counter), rule_name, "formats", ["png"])
    return [fmts] if isinstance(fmts, str) else [str(f) for f in fmts]


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
    benchmark: logdir + "/{sample_id}/scanpy_qc_{counter}.benchmark.tsv",
    resources:
        mem_mb=rule_res("scanpy_qc").get("mem_mb", 16000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml", "scanpy_qc")
    params:
        python=python,
        script=script,
        # 参数精确到 rule:Params.scanpy.<counter>.scanpy_qc.<key>
        min_genes=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_qc", "min_genes", 200),
        max_genes=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_qc", "max_genes", 6000),
        max_pct_mt=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_qc", "max_pct_mt", 20),
        n_top_genes=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_qc", "n_top_genes", 3000),
        use_mad=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_qc", "use_mad", False),
        scrublet=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_qc", "scrublet", True),
        doublet_rate=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_qc", "doublet_rate", 0.06),
        formats=lambda wildcards: _image_formats(wildcards, "scanpy_qc")
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
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
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
    benchmark: logdir_combine + "/scanpy/{tissue}/scanpy_merge_{counter}.benchmark.tsv",
    resources:
        mem_mb=rule_res("scanpy_merge").get("mem_mb", 32000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml", "scanpy_merge")
    params:
        python=python,
        script=script,
        te_bed=te_bed,
        gene_tsv=gene_tsv,
        # 参数精确到 rule:Params.scanpy.<counter>.scanpy_merge.formats
        formats=lambda wildcards: _image_formats(wildcards, "scanpy_merge")
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
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
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
    benchmark: logdir_combine + "/scanpy/{tissue}/scanpy_auto_{counter}.benchmark.tsv",
    resources:
        mem_mb=rule_res("scanpy_auto").get("mem_mb", 64000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml", "scanpy_auto")
    params:
        python=python,
        script=script,
        tissue=lambda wildcards: wildcards.tissue,
        # 参数精确到 rule:Params.scanpy.<counter>.scanpy_auto.<key>
        # (原 auto/cluster 参数段逐 key 挂到 scanpy_auto 名下)
        llm_method=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "llm_method", ""),
        llm_model=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "llm_model", ""),
        llm_api_key=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "llm_api_key", ""),
        llm_base_url=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "llm_base_url", ""),
        resolution=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "resolution", 0.8),
        max_iterations=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "max_iterations", 3),
        min_genes=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "min_genes", 800),
        min_counts=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "min_counts", 3000),
        max_pct_mt=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "max_pct_mt", 20),
        n_pcs=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "n_pcs", 50),
        n_neighbors=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "n_neighbors", 50),
        n_top_genes=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "n_top_genes", 3000),
        batch_method=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "batch_method", "harmony"),
        batch_key=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "batch_key", "sample_id"),
        auto_n_pcs=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "auto_n_pcs", True),
        skip_te=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "skip_te", True),
        regress_cell_cycle=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "regress_cell_cycle", True),
        save_iteration_h5ad=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "save_iteration_h5ad", False),
        # marker_supplement 原为 scanpy_auto/scanpy_reannotate 共享 key;
        # 现按 rule 独立(Params.scanpy.<counter>.scanpy_auto.marker_supplement)
        marker_supplement=lambda wildcards: (rule_param(params.get(wildcards.counter), "scanpy_auto", "marker_supplement", {}) or {}).get(wildcards.tissue, ""),
        marker_supplement_llm=lambda wildcards: bool(rule_param(params.get(wildcards.counter), "scanpy_auto", "marker_supplement_llm", True)),
        exclude_genes=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_auto", "exclude_genes", "") or "",
        exclude_from_deg=lambda wildcards: bool(rule_param(params.get(wildcards.counter), "scanpy_auto", "exclude_from_deg", False)),
        formats=lambda wildcards: _image_formats(wildcards, "scanpy_auto"),
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
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
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

def _reannotate_override(rule_name):
    """Override JSON path for this tissue/counter under the rule's own key,
    or [] when disabled.

    参数精确到 rule:Params.scanpy.<counter>.<rule_name>.cluster_annotations
    (原 auto.cluster_annotations 共享 key 已按 rule 拆分;需要 override 或
    排序的 rule 各自配置自己的子段,与 scanpy_reannotate 保持同值以保证
    DAG 一致 — marker 生产端(scanpy_reannotate)与消费端
    (scanpy_advanced / scanpy_pseudobulk)按各自 key 独立解析)。
    """
    def _get(wildcards):
        ca = rule_param(params.get(wildcards.counter), rule_name, "cluster_annotations", {}) or {}
        return ca.get(wildcards.tissue, "") or []
    return _get


def _reannotate_marker(rule_name):
    """Marker file input enforcing auto -> reannotate -> consumers ordering.

    The reannotate step overwrites the auto h5ad in place; downstream rules
    (advanced, pseudobulk) depend on this marker so they never read a
    pre-override auto.h5ad. Empty when the consuming rule's own
    cluster_annotations key has no entry for the tissue. The condition is
    read from ``rule_name``'s sub-section (参数精确到 rule) — configure the
    same overrides under the consumer rule name as under scanpy_reannotate
    so the marker is requested whenever reannotate runs.
    """
    def _get(wildcards):
        return (outdir_combine + f"/{wildcards.tissue}/{wildcards.tissue}_{wildcards.counter}_auto.h5ad.reannotated"
                if _reannotate_override(rule_name)(wildcards) else [])
    return _get


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
        overrides = _reannotate_override("scanpy_reannotate")
    output:
        marker = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad.reannotated",
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_reannotate_{counter}.log"
    threads: rule_res("scanpy_reannotate").get("threads", 2)
    benchmark: logdir_combine + "/scanpy/{tissue}/scanpy_reannotate_{counter}.benchmark.tsv",
    resources:
        mem_mb=rule_res("scanpy_reannotate").get("mem_mb", 16000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml", "scanpy_reannotate")
    params:
        python=python,
        script=script,
        # 参数精确到 rule:Params.scanpy.<counter>.scanpy_reannotate.<key>
        # (marker_supplement 原与 scanpy_auto 共享,现各自独立 key)
        marker_supplement=lambda wildcards: (rule_param(params.get(wildcards.counter), "scanpy_reannotate", "marker_supplement", {}) or {}).get(wildcards.tissue, ""),
        formats=lambda wildcards: _image_formats(wildcards, "scanpy_reannotate"),
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
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
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
        reannotated = _reannotate_marker("scanpy_advanced")
    output:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_advanced.h5ad",
        plot_dir = directory(outdir_combine + "/{tissue}/plots/{counter}/advanced")
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_advanced_{counter}.log"
    threads: rule_res("scanpy_advanced").get("threads", 4)
    benchmark: logdir_combine + "/scanpy/{tissue}/scanpy_advanced_{counter}.benchmark.tsv",
    resources:
        mem_mb=rule_res("scanpy_advanced").get("mem_mb", 64000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml", "scanpy_advanced")
    params:
        python=python,
        script=script,
        # 参数精确到 rule:Params.scanpy.<counter>.scanpy_advanced.<key>
        # (原 cluster/advanced 参数段逐 key 挂到 scanpy_advanced 名下)
        n_pcs=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "n_pcs", 50),
        n_neighbors=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "n_neighbors", 15),
        trajectory=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "trajectory", True),
        velocity=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "velocity", False),
        communication=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "communication", False),
        cnv=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "cnv", False),
        gtf=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "gtf", ""),
        cnv_reference=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_advanced", "cnv_reference", ""),
        formats=lambda wildcards: _image_formats(wildcards, "scanpy_advanced")
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
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
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
        #   Pseudo-bulk DEG (PyDESeq2) — modules/scanpy/bin/pseudobulk.py
        #   Two-file design: --h5ad (scanpy_auto output) supplies the cell_type
        #   annotations; --counts-h5ad (scanpy_merge output) supplies the
        #   full-gene integer layers["counts"] (the auto file's counts layer is
        #   HVG-subset and its raw slot is log1p-normalised — neither is valid
        #   pseudobulk input). Runs PyDESeq2 per compare entry, plots
        #   volcano/heatmap/GO. 参数精确到 rule:compare 等逐 key 配置于
        #   Params.scanpy.<counter>.scanpy_pseudobulk.<key>(rule_param 见
        #   common.smk;原 pseudobulk 整段 dict 转发已展开逐 key,无共享回退)。
        # ---------------------------------------------------------------------------
rule scanpy_pseudobulk:
    """Pseudo-bulk DEG via PyDESeq2: annotations from scanpy_auto, full-gene
    raw counts from the merged h5ad. When overrides are configured for the
    tissue, the reannotate marker orders pb after the in-place label fix."""
    input:
        h5ad = outdir_combine + "/{tissue}/{tissue}_{counter}_auto.h5ad",
        reannotated = _reannotate_marker("scanpy_pseudobulk"),
        counts = outdir_combine + "/{tissue}/{tissue}_{counter}_merged.h5ad"
    output:
        report = outdir_pseudobulk + "/{tissue}/{counter}/pseudobulk_report.html",
    log:
        logdir_combine + "/scanpy/{tissue}/scanpy_pseudobulk_{counter}.log"
    threads: rule_res("scanpy_pseudobulk").get("threads", 4)
    benchmark: logdir_combine + "/scanpy/{tissue}/scanpy_pseudobulk_{counter}.benchmark.tsv",
    resources:
        mem_mb=rule_res("scanpy_pseudobulk").get("mem_mb", 32000),
    conda:
        "scanpy.yaml"
    container:
        sif("scanpy.yaml", "scanpy_pseudobulk")
    params:
        # 参数精确到 rule:Params.scanpy.<counter>.scanpy_pseudobulk.<key>
        # (原 pb_cfg 整段 dict 转发已展开逐 key,rule_param 见 common.smk)
        pb_min_cells=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "min_cells", 30),
        pb_min_counts=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "min_counts", 1000),
        pb_n_pseudo_reps=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "n_pseudo_reps", 3),
        pb_min_cells_per_rep=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "min_cells_per_rep", 10),
        pb_cell_types=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "cell_types", None),
        pb_compare=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "compare", []) or [],
        pb_top_cell_types=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "top_cell_types", []) or [],
        pb_highlight_genes=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "highlight_genes", []) or [],
        pb_deg_mode=lambda wildcards: rule_param(params.get(wildcards.counter), "scanpy_pseudobulk", "deg_mode", ["all"]),
        formats = lambda wildcards: _image_formats(wildcards, "scanpy_pseudobulk"),
        python = python,
        script = os.path.join(ROOT_DIR, "modules", "scanpy", "bin", "pseudobulk.py"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("scanpy_pseudobulk", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scanpy pseudobulk for tissue {wildcards.tissue} at {current_time}")
            out_root = os.path.dirname(str(output.report))
            os.makedirs(out_root, exist_ok=True)
            # dual-write receipts: timestamped copy in log/cmd, plain copy in outdir
            script = os.path.join(RECEIPT_DIR, f"scanpy_pseudobulk_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            plain_script = os.path.join(out_root, f"scanpy_pseudobulk_{wildcards.counter}_{wildcards.tissue}.sh")
            cmd = [params.python, params.script,
                   "--h5ad", f"{wildcards.tissue}={input.h5ad}",
                   "--counts-h5ad", f"{wildcards.tissue}={input.counts}",
                   "--out-dir", out_root,
                   "--min-cells", str(params.pb_min_cells),
                   "--min-counts", str(params.pb_min_counts),
                   "--n-pseudo-reps", str(params.pb_n_pseudo_reps),
                   "--min-cells-per-rep", str(params.pb_min_cells_per_rep)]
            # --cell-types is nargs='+': ONE flag carrying all values —
            # repeated flags overwrite each other in argparse.
            if params.pb_cell_types:
                cmd += ["--cell-types"] + [str(ct) for ct in params.pb_cell_types]
            # entries: tissue:ref_sample:treat_sample[=label] — every entry
            # is scoped to ONE tissue (pseudobulk.py rejects a compare whose
            # tissue has no --h5ad definition), so this per-tissue job only
            # takes the entries prefixed with its own wildcard.
            tissue_compares = [
                str(c) for c in (params.pb_compare or [])
                if str(c).split(":", 1)[0] == wildcards.tissue
            ]
            if not tissue_compares:
                raise RuntimeError(
                    f"Params.scanpy.{wildcards.counter}.scanpy_pseudobulk.compare "
                    f"has no entry for tissue '{wildcards.tissue}' — expected "
                    f"entries like '{wildcards.tissue}:ref_sample:"
                    f"treat_sample[=label]'")
            for cmp_ in tissue_compares:
                cmd += ["--compare", cmp_]
            for ct in (params.pb_top_cell_types or []):
                cmd += ["--top-cell-types", str(ct)]
            for g in (params.pb_highlight_genes or []):
                cmd += ["--highlight-genes", str(g)]
            # --deg-mode is nargs='+': ONE flag carrying all values —
            # repeated flags overwrite each other in argparse and only the
            # last mode would run.
            cmd += ["--deg-mode"] + [str(m) for m in (params.pb_deg_mode or ["all"])]
            # --format is action='append': one flag per format (modules.md #9)
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
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
