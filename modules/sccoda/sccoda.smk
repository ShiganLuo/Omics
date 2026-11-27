
include: "../common/common.smk"
import json

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")

tissue_samples = config.get("tissue_samples", {})
tissues = sorted(tissue_samples.keys())
counters = sorted(k for k, v in (config.get("Params", {}).get("sccoda", {}) or {}).items()
                  if isinstance(v, dict) and k != "tissue_samples")


def _sccoda_override(wildcards):
    """Path of the tissue's reannotate override JSON, [] when none."""
    entry = config.get("reannotate", {}) .get(wildcards.counter, {}) or {}
    if not isinstance(entry, dict):
        return []
    return entry.get(wildcards.tissue, "") or []


def _reannotate_marker(wildcards):
    """Marker file input enforcing auto -> reannotate -> sccoda ordering.

    The reannotate step overwrites the auto h5ad in place; sccoda depends on
    this marker so it never reads a pre-override auto.h5ad. Empty when no
    overrides are configured (the marker is then never produced)."""
    return (indir + f"/{wildcards.tissue}/{wildcards.tissue}_{wildcards.counter}_auto.h5ad.reannotated"
            if _sccoda_override(wildcards) else [])


# ---------------------------------------------------------------------------
# scCODA differential abundance over annotated cell types
# ---------------------------------------------------------------------------
rule sccoda_run:
    """scCODA differential abundance (pertpy): one job per tissue x counter,
    running the sample_params.design comparison design scoped to its tissue
    (multi-group entries as ordered 'groups' bars, pairwise as group_pairs).
    Consumes the scanpy_auto h5ad (obs: sample_id + cell_type); when cluster
    annotation overrides are configured for the tissue, the reannotate marker
    orders this rule after the in-place label fix."""
    input:
        h5ad = indir + "/{tissue}/{tissue}_{counter}_auto.h5ad",
        reannotated = _reannotate_marker
    output:
        report = outdir + "/{tissue}/{counter}/sccoda_report.html",
        da_dir = directory(outdir + "/{tissue}/{counter}/da")
    log:
        logdir + "/sccoda/{tissue}/sccoda_run_{counter}.log"
    threads: rule_res("sccoda_run").get("threads", 4)
    benchmark: logdir + "/sccoda/{tissue}/sccoda_run_{counter}.benchmark.tsv",
    resources:
        mem_mb=rule_res("sccoda_run").get("mem_mb", 32000),
    conda:
        "sccoda.yaml"
    container:
        sif("sccoda.yaml", "sccoda_run")
    params:
        python = config.get("Procedure", {}).get("python") or "python",
        script = os.path.join(ROOT_DIR, "modules", "sccoda", "bin", "sccoda.py"),
        ref_celltype=lambda wildcards: config.get("Params", {}).get("sccoda", {}).get(wildcards.counter, {}).get("sccoda_run", {}).get("ref_celltype", "automatic"),
        formula=lambda wildcards: config.get("Params", {}).get("sccoda", {}).get(wildcards.counter, {}).get("sccoda_run", {}).get("formula", "condition"),
        credible_threshold=lambda wildcards: config.get("Params", {}).get("sccoda", {}).get(wildcards.counter, {}).get("sccoda_run", {}).get("credible_threshold", 0.90),
        num_samples=lambda wildcards: config.get("Params", {}).get("sccoda", {}).get(wildcards.counter, {}).get("sccoda_run", {}).get("num_samples", 2000),
        num_warmup=lambda wildcards: config.get("Params", {}).get("sccoda", {}).get(wildcards.counter, {}).get("sccoda_run", {}).get("num_warmup", 500),
        sc_design=lambda wildcards: config.get("design", {}) or {},
        tissue_samples=lambda wildcards: config.get("tissue_samples", {}) or {},
        formats=lambda wildcards: (lambda fmts: [fmts] if isinstance(fmts, str) else [str(f) for f in fmts])(
            config.get("Params", {}).get("sccoda", {}).get(wildcards.counter, {}).get("sccoda_run", {}).get("formats") or ["png"]),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sccoda_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scCODA differential abundance for tissue {wildcards.tissue} at {current_time}")
            out_root = os.path.dirname(str(output.report))
            os.makedirs(out_root, exist_ok=True)
            script = receipt_script(f"sccoda_run_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(out_root, f"sccoda_run_{wildcards.counter}_{wildcards.tissue}.sh")
            cmd = [params.python, params.script,
                   "--h5ad", f"{wildcards.tissue}={input.h5ad}",
                   "--out-dir", out_root,
                   "--ref-celltype", str(params.ref_celltype or "automatic"),
                   "--formula", str(params.formula or "condition"),
                   "--credible-threshold", str(params.credible_threshold),
                   "--num-samples", str(params.num_samples),
                   "--num-warmup", str(params.num_warmup)]
            # 比较设计:sample_params.design → comparison_design(native 视图:
            # 两两 entry → group_pairs,多组 entry → groups 有序形状原样,
            # 首组 = 基线)→ compare_design.json(回执/溯源)→
            # --compare-json。转化集中在 common.smk::write_compare_design,
            # 语义错误(缺组织/组名悬空/样本互斥/多组 ref 非首组)在
            # 那里快速失败。
            design_json = os.path.join(out_root, "compare_design.json")
            write_compare_design(params.sc_design, wildcards.tissue,
                                 "native", design_json,
                                 tissue_samples=params.tissue_samples)
            cmd += ["--compare-json", design_json]
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
                f.write(f'echo "scCODA differential abundance for {wildcards.tissue} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scCODA differential abundance for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scCODA differential abundance for tissue {wildcards.tissue}: {e}")

rule sccoda_result:
    """Panel rule aggregating the module's main outputs."""

    input:
        expand(outdir + "/{tissue}/{counter}/sccoda_report.html",
               tissue=tissues, counter=counters)
