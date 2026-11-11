include: "../common/common.smk"


outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", []) or []
# tumor -> normal design pairs (same shape as Params.workflow.paired_design
# consumed by gatk_somatic/msi: {"experimental": t, "control": n})
paired_design = config.get("paired_design", []) or []
# Registered workflow outputs (node.py outfiles, passed through by the
# subworkflow). Used as DAG anchors: every registered file under a consumed
# directory becomes a tracked input so this report runs AFTER its producers
# and re-runs when a section output is refreshed (qc_summary's qc_anchors
# semantics). Optional stages that were never registered contribute nothing,
# so no MissingInputException.
outfiles = config.get("outfiles") or []

# Upstream result directories — subworkflow passes xxx_config["outdir"]
# values (subworkflow.md chaining rule); the script tolerates absent dirs.
qc_dir = config.get("qc_summary") or os.path.join(outdir, "QC/qc_summary")
somatic_dir = config.get("somatic_dir") or os.path.join(outdir, "variation/somatic_snv_indel")
spectrum_dir = config.get("spectrum_dir") or os.path.join(outdir, "stability/spectrum")
burden_dir = config.get("burden_dir") or os.path.join(outdir, "stability/burden")
msi_dir = config.get("msi_dir") or os.path.join(outdir, "stability/msi")
fragment_dir = config.get("fragment_dir") or os.path.join(outdir, "results/fragment_size/fragment")

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "Mutation_report", "bin", "generate_report.py")


def _report_anchors():
    """Registered outputs under the consumed dirs -> DAG anchor files.

    The report renders everything it finds under these directories (qc
    tables, spectrum/signature figures and CSVs, burden, msi, fragment);
    anchoring on their registered outputs both orders this rule after the
    producers and invalidates the report when any of them is regenerated.
    Files outside the consumed dirs (e.g. this report's own outputs) are
    ignored. Only registered outputs enter the list, so optional stages
    that never ran cannot raise MissingInputException.
    """
    dirs = [d for d in (qc_dir, somatic_dir, spectrum_dir, burden_dir,
                        msi_dir, fragment_dir) if d]
    return [f for f in outfiles
            if any(str(f).startswith(str(d).rstrip("/") + "/") for d in dirs)]

rule generate_report:
    """Generate the Mutation summary PPT report (QC + somatic calls +
    spectrum + signatures + MSI + burden + fragment size) and the results
    workbook (one sheet per module + file inventory). Optional sections
    degrade to 'not run' notes when their upstream dirs are absent."""
    input:
        # QC summary table + every registered output under the consumed
        # dirs anchor this report after their producing rules (DAG
        # ordering + invalidation; the script itself re-discovers files
        # on disk and degrades to 'not run' notes for absent stages)
        anchors = _report_anchors(),
        qc_tsv = qc_dir + "/qc_summary.tsv",
    output:
        report = outdir + "/Mutation_report.pptx",
        excel = outdir + "/Mutation_report.xlsx",
    log:
        logdir + "/Mutation_report/generate_report.log"
    threads: rule_res("generate_report").get("threads", 1)
    resources:
        mem_mb = rule_res("generate_report").get("mem_mb", 8000),
    conda:
        "Mutation_report.yaml"
    container:
        sif("Mutation_report.yaml", "generate_report")
    params:
        script = REPORT_SCRIPT,
        title = config.get("Params", {}).get("report", {}).get("title")
            or "Mutation 分析报告 / Mutation Report",
        pipeline = config.get("Params", {}).get("report", {}).get("pipeline") or "Mutation",
        genome = config.get("Params", {}).get("report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("date") or "",
        pairs = [
            (rec.get("experimental"), rec.get("control"))
            for rec in paired_design
            if rec.get("experimental") and rec.get("control")
        ],
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("Mutation_generate_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start Mutation generate_report at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            if report_dir:
                os.makedirs(report_dir, exist_ok=True)

            cmd = [
                "python3", params.script,
                "--qc-dir", qc_dir,
                "--somatic-dir", somatic_dir,
                "--spectrum-dir", spectrum_dir,
                "--burden-dir", burden_dir,
                "--msi-dir", msi_dir,
                "--fragment-dir", fragment_dir,
                "--output", str(output.report),
                "--excel-output", str(output.excel),
                "--title", params.title,
                "--pipeline", params.pipeline,
                "--genome", params.genome,
                "--date", params.date,
                "--log", log_path,
            ]
            for s in samples:
                cmd += ["--sample", s]
            for tumor, normal in params.pairs:
                cmd += ["--pair", f"{tumor}:{normal}"]

            script = os.path.join(RECEIPT_DIR, f"generate_report_{current_time}.sh")
            plain_script = os.path.join(report_dir or ".", "generate_report.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write('echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"Mutation report completed: {output.report}")
        except Exception as e:
            rule_logger.error(f"Error occurred during Mutation generate_report: {e}")
            raise RuntimeError(f"Error occurred during Mutation generate_report: {e}, you can check the log file {log_path} for more details.")

rule report_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        report = outdir + "/Mutation_report.pptx",
        excel = outdir + "/Mutation_report.xlsx",
    resources:
        mem_mb = rule_res("report_result").get("mem_mb", 8000),
