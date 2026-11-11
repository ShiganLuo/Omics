include: "../common/common.smk"


outdir = config.get("outdir", "output/QC/qc_summary")
logdir = config.get("logdir", "log")
samples = config.get("samples", []) or []
# stage name -> QC artifact directory; subworkflow passes upstream
# xxx_config["outdir"] values here (subworkflow.md chaining rule)
qc_dirs = config.get("qc_dirs", {}) or {}
# stage -> per-sample anchor filename template ('{sample}' placeholder),
# ONLY for stages whose producing rules are guaranteed to run (registered
# outputs or unavoidable intermediates). Anchors become tracked rule inputs
# so Snakemake orders qc_summary after the producers AND re-runs it when an
# upstream artifact is refreshed. Stages listed in qc_dirs but absent here
# (e.g. fastqc outputs the workflow never registers) stay best-effort:
# the script scans the directory and skips it with a warning when missing.
qc_anchors = config.get("qc_anchors", {}) or {}
# optional per-stage sample list override (default: all samples) — e.g. MSI
# anchors exist only for tumor samples in the design pairs
qc_anchor_samples = config.get("qc_anchor_samples", {}) or {}
report_title = config.get("Params", {}).get("qc_summary", {}).get("title") or "QC Summary"

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
        # anchor the aggregation AFTER per-sample QC rules when the stage
        # directories are known outputs (DAG ordering only; the script
        # re-discovers files inside each directory)
        anchors = _qc_stage_inputs()
    output:
        tsv = outdir + "/qc_summary.tsv",
        html = outdir + "/qc_summary.html",
    log:
        logdir + "/qc_summary/qc_summary.log"
    threads: rule_res("qc_summary").get("threads", 1)
    resources:
        mem_mb = rule_res("qc_summary").get("mem_mb", 8000),
    conda:
        "qc_summary.yaml"
    container:
        sif("qc_summary.yaml", "qc_summary")
    params:
        script = QC_SCRIPT,
        title = report_title,
    run:
        log_path = str(log)
        open(log_path, "w").close()
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
            cmd += ["--log", log_path]

            script = os.path.join(RECEIPT_DIR, f"qc_summary_{current_time}.sh")
            plain_script = os.path.join(summary_dir, "qc_summary.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write('echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
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
    resources:
        mem_mb = rule_res("qc_summary_result").get("mem_mb", 8000),
