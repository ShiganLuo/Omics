include: "../common/common.smk"


outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", []) or []
dna_samples = config.get("dna_samples", []) or []
rna_samples = config.get("rna_samples", []) or []
paired_samples = config.get("paired_samples", []) or []
single_samples = config.get("single_samples", []) or []
# tumor -> normal design pairs (same shape as Params.workflow.paired_design
# consumed by gatk_somatic/msi: {"experimental": t, "control": n})
paired_design = config.get("paired_design", []) or []
contrasts = config.get("contrasts", []) or []
# Registered workflow outputs (node.py outfiles, passed through by the
# subworkflow). Used as DAG anchors (Mutation_report 同款): every registered
# file under a consumed root becomes a tracked input so this report runs
# AFTER its producers and re-runs when a section output is refreshed.
# Optional stages that were never registered contribute nothing, so no
# MissingInputException.
outfiles = config.get("outfiles") or []

# Upstream roots — subworkflow passes xxx_config["outdir"] values
# (subworkflow.md chaining rule); the script tolerates absent dirs/stages.
dna_dir = config.get("dna_dir") or os.path.join(outdir, "DNA")
rna_dir = config.get("rna_dir") or os.path.join(outdir, "RNA")
drclass_dir = config.get("drclass_dir") or os.path.join(outdir, "DRclass")
teffect_dir = config.get("teffect_dir") or os.path.join(outdir, "TransEffect")
burden_dir = config.get("burden_dir") or os.path.join(dna_dir, "stability", "burden")

# RNA counts/DE/function trees are genome-labelled (polygenomes layout):
# the report covers the default genome's tree. The label is passed
# explicitly by the subworkflow (node.py genome.default), so multi-genome
# runs still resolve the intended RNA tree; Params.report.generate_report
# .genome (subtitle) is only a fallback.
genome_label = config.get("genome_label") or config.get("Params", {}).get("report", {}).get("generate_report", {}).get("genome") or "all"

# 参数精确到 rule:generate_report 自己的参数挂 Params.report.generate_report.<key>
# (无 module 级共享回退,读点写 config 全量路径字面链)。
_formats = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_formats, str):
    _formats = [_formats]

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "DR_report", "bin", "generate_report.py")


def _report_anchors():
    """Registered outputs under the consumed roots -> DAG anchor files.

    The report renders everything it finds under these roots (DNA QC /
    variation / stability, RNA counts/DE/function/fusion/transcripts, and
    the DR-unique DRclass / TransEffect / burden outputs); anchoring on
    their registered outputs both orders this rule after the producers and
    invalidates the report when any of them is regenerated. Files outside
    the consumed roots (e.g. this report's own outputs) are ignored. Only
    registered outputs enter the list, so optional stages that never ran
    cannot raise MissingInputException.
    """
    dirs = [d for d in (dna_dir, rna_dir, drclass_dir, teffect_dir, burden_dir) if d]
    return [f for f in outfiles
            if any(str(f).startswith(str(d).rstrip("/") + "/") for d in dirs)]


rule generate_report:
    """Generate the DR joint summary PPT report (DNA chain + RNA chain +
    DR-unique analyses: burden / DRclass / TransEffect) + results workbook."""
    input:
        anchors = _report_anchors(),
        dna_qc_tsv = dna_dir + "/QC/qc_summary/qc_summary.tsv",
        rna_qc_tsv = rna_dir + "/QC/qc_summary/qc_summary.tsv",
    output:
        report = outdir + "/DR_report.pptx",
        file_inventory = outdir + "/DR_report_files.xlsx",
        # 报告自有图 figures/<name>.<fmt> 由脚本落盘，DAG 用文件清单
        # manifest.tsv 锚定（no directory() outputs）
        figures_manifest = outdir + "/figures/manifest.tsv",
    log:
        logdir + "/DR_report.log"
    threads: rule_res("generate_report").get("threads", 1)
    benchmark: logdir + "/DR_report.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "generate_report", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "DR_report.yaml"
    container:
        sif("DR_report.yaml", "generate_report")
    params:
        samples = samples,
        dna_samples = dna_samples,
        rna_samples = rna_samples,
        paired_samples = paired_samples,
        single_samples = single_samples,
        contrasts = contrasts,
        pairs = [
            f"{p.get('experimental')}:{p.get('control')}"
            for p in paired_design if p.get("control") and p.get("experimental")
        ],
        title = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("title") or "DNA×RNA 联合分析报告 / DR Joint Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("genome") or "",
        genome_label = genome_label,
        date = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("date") or "",
        lang = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("lang") or "zh",
        formats = _formats,
        img_dir = outdir + "/figures",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("DR_generate_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start DR report generation at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            os.makedirs(params.img_dir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"DR_report_{current_time}.sh")
            os.makedirs(os.path.dirname(script_path), exist_ok=True)
            plain_script = os.path.join(report_dir, f"DR_report.sh")
            cmd = [
                "python3", params.script,
                "--analysis-dir", outdir,
                "--dna-dir", dna_dir,
                "--rna-dir", rna_dir,
                "--drclass-dir", drclass_dir,
                "--teffect-dir", teffect_dir,
                "--burden-dir", burden_dir,
                "--output", str(output.report),
                "--file-inventory", str(output.file_inventory),
                "--title", params.title,
                "--subtitle", params.subtitle,
                "--pipeline", params.pipeline,
                "--genome", params.genome,
                "--genome-label", params.genome_label,
                "--date", params.date,
                "--lang", params.lang,
                "--img-dir", params.img_dir,
            ]
            if params.samples:
                for sid in params.samples:
                    cmd.extend(["--sample", str(sid)])
            if params.dna_samples:
                for sid in params.dna_samples:
                    cmd.extend(["--dna-sample", str(sid)])
            if params.rna_samples:
                for sid in params.rna_samples:
                    cmd.extend(["--rna-sample", str(sid)])
            if params.paired_samples:
                for sid in params.paired_samples:
                    cmd.extend(["--paired-sample", str(sid)])
            if params.single_samples:
                for sid in params.single_samples:
                    cmd.extend(["--single-sample", str(sid)])
            if params.contrasts:
                for contrast in params.contrasts:
                    cmd.extend(["--contrast", str(contrast)])
            for pair in params.pairs:
                cmd.extend(["--pair", str(pair)])
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join([shlex.quote(str(x)) for x in cmd]) + "\n")
                f.write("echo 'DR report generation completed successfully.'\n")
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"DR report generation failed: {e}")
            raise RuntimeError(f"DR report generation failed: {e}, please check the log file: {log_path} for details.")

rule report_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        report = outdir + "/DR_report.pptx",
        file_inventory = outdir + "/DR_report_files.xlsx",
        figures_manifest = outdir + "/figures/manifest.tsv",
