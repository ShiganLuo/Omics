include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
input_samples = config.get("input_samples", [])

REPORT_SCRIPT = os.path.join(os.path.dirname(workflow.snakefile), "bin", "generate_report.py")

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
    shell:
        """
        python3 {params.script} \
            --samples {params.samples} \
            --input-samples {params.input_samples} \
            --peaks-dir {params.peaks_dir} \
            --annotation-dir {params.annotation_dir} \
            --qc-dir {params.qc_dir} \
            --log-dir {params.log_dir} \
            --markdup-dir {params.markdup_dir} \
            --output {output.report} \
            --img-dir {params.img_dir} \
            --title "{params.title}" \
            --subtitle "{params.subtitle}" \
            --pipeline "{params.pipeline}" \
            --genome "{params.genome}" \
            --date "{params.date}" \
            --top-n {params.top_n} \
            --lang {params.lang} \
            > {log} 2>&1
        """


rule report_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        report = outdir + "/PeakCalling_report.pptx"
