include: "../common/common.smk"


outdir = config.get("outdir", "output/QC/qc_summary")
logdir = config.get("logdir", "log")
samples = config.get("samples", []) or []
qc_dirs = config.get("qc_dirs", {}) or {}
qc_anchors = config.get("qc_anchors", {}) or {}
qc_anchor_samples = config.get("qc_anchor_samples", {}) or {}
report_title = config.get("Params", {}).get("qc_summary", {}).get("title") or "QC Summary"
_formats = (config.get("Params", {}).get("qc_summary", {}) or {}).get("formats") or ["png"]
if isinstance(_formats, str):
    _formats = [_formats]

QC_SCRIPT = os.path.join(ROOT_DIR, "modules", "qc_summary", "bin", "qc_summary.py")


def _qc_stage_inputs():
    """Concrete anchor files for the anchored QC stages.

    Per-sample templates (containing '{sample}') expand against each
    stage's directory (from qc_dirs) and its sample list
    (qc_anchor_samples override or all samples). Templates WITHOUT the
    placeholder are run-level anchors (one table for the whole run, e.g.
    sexcheck.tsv) and are joined as a single file. Entries for stages
    without a dir are skipped defensively.
    """
    files = []
    for stage, template in qc_anchors.items():
        stage_dir = qc_dirs.get(stage)
        if not stage_dir:
            continue
        template = str(template)
        if "{sample}" not in template:
            files.append(os.path.join(stage_dir, template))
            continue
        stage_samples = qc_anchor_samples.get(stage, samples)
        files += [os.path.join(stage_dir, template.format(sample=s))
                  for s in stage_samples]
    return files


rule qc_summary:
    """Aggregate per-sample QC artifacts from all pipeline stages into one
    sample x metric table (qc_summary.tsv) plus an HTML visual summary."""
    input:
        anchors = _qc_stage_inputs()
    output:
        tsv = outdir + "/qc_summary.tsv",
        html = outdir + "/qc_summary.html",
    log:
        logdir + "/qc_summary/qc_summary.log"
    threads: rule_res("qc_summary").get("threads", 1)
    benchmark: logdir + "/qc_summary/qc_summary.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "qc_summary", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "qc_summary.yaml"
    container:
        sif("qc_summary.yaml", "qc_summary")
    params:
        script = QC_SCRIPT,
        title = report_title,
        formats = _formats,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("qc_summary", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start qc_summary at {current_time}")
            summary_dir = os.path.dirname(str(output.tsv))
            os.makedirs(summary_dir, exist_ok=True)

            cmd = [
                "python3", params.script,
                "--outdir", summary_dir,
                "--title", params.title,
            ]
            for s in samples:
                cmd += ["--sample", s]
            for stage, stage_dir in qc_dirs.items():
                cmd += ["--qc-dir", f"{stage}={stage_dir}"]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            cmd += ["--log", log_path]

            script = receipt_script(f"qc_summary_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(summary_dir, "qc_summary.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write('echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"qc_summary completed: {output.tsv}")
        except Exception as e:
            rule_logger.error(f"Error occurred during qc_summary: {e}")
            raise RuntimeError(f"Error occurred during qc_summary: {e}, you can check the log file {log_path} for more details.")


rule qc_summary_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        tsv = outdir + "/qc_summary.tsv",
        html = outdir + "/qc_summary.html",
