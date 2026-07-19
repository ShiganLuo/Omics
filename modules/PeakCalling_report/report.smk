include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
input_samples = config.get("input_samples", [])

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "PeakCalling_report", "bin", "generate_report.py")

rule generate_report:
    """
    Generate ChIP-seq Peak Calling Report (PPT).
    All content is parameterized via config — no hardcoded sample names or species.
    """
    input:
        peaks = expand(outdir + "/peaks/{sample}/{sample}_peaks.narrowPeak", sample=samples),
        frip = expand(outdir + "/QC/3_frip_score/{sample}/{sample}.FRiP.txt", sample=samples),
        annotations = expand(outdir + "/annotation/{sample}/{sample}_peaks.annotatePeaks.txt", sample=samples),
        bowtie2_logs = expand(outdir + "/log/{sample}/bowtie2_align.log", sample=samples + input_samples),
        markdup_metrics = expand(outdir + "/common/4_markdup_bam/{sample}/{sample}.Markdup-metrics.txt", sample=samples + input_samples),
    output:
        report = outdir + "/PeakCalling_report.pptx"
    log:
        logdir + "/report.log"
    threads: 1
    conda:
        "report.yaml"
    params:
        samples = " ".join(samples),
        input_samples = " ".join(input_samples),
        peaks_dir = outdir + "/peaks",
        annotation_dir = outdir + "/annotation",
        qc_dir = outdir + "/QC/3_frip_score",
        log_dir = outdir + "/log",
        markdup_dir = outdir + "/common/4_markdup_bam",
        img_dir = outdir + "/ppt_results",
        title = config.get("Params", {}).get("report", {}).get("title") or "",
        subtitle = config.get("Params", {}).get("report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("date") or "",
        top_n = config.get("Params", {}).get("report", {}).get("top_n") or 5,
        lang = config.get("Params", {}).get("report", {}).get("lang") or "zh",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("generate_report", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start generate_report at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            script = os.path.join(report_dir, f"generate_report_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"python3 {params.script} \\\n")
                f.write(f"    --samples {params.samples} \\\n")
                f.write(f"    --input-samples {params.input_samples} \\\n")
                f.write(f"    --peaks-dir {params.peaks_dir} \\\n")
                f.write(f"    --annotation-dir {params.annotation_dir} \\\n")
                f.write(f"    --qc-dir {params.qc_dir} \\\n")
                f.write(f"    --log-dir {params.log_dir} \\\n")
                f.write(f"    --markdup-dir {params.markdup_dir} \\\n")
                f.write(f"    --output {output.report} \\\n")
                f.write(f"    --img-dir {params.img_dir} \\\n")
                f.write(f"    --title \"{params.title}\" \\\n")
                f.write(f"    --subtitle \"{params.subtitle}\" \\\n")
                f.write(f"    --pipeline \"{params.pipeline}\" \\\n")
                f.write(f"    --genome \"{params.genome}\" \\\n")
                f.write(f"    --date \"{params.date}\" \\\n")
                f.write(f"    --top-n {params.top_n} \\\n")
                f.write(f"    --lang {params.lang}\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during generate_report: {e}\n")
            logger.error(f"Error occurred during generate_report: {e}")
            raise e


rule report_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        report = outdir + "/PeakCalling_report.pptx"
