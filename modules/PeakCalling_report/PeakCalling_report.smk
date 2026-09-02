include: "../common/common.smk"
import shlex
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ip_samples = config.get("samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})

# Input directories (allow overriding for subworkflow import where
# peaks/annotation live under a different subdirectory than outdir)
peaks_dir = config.get("peaks_dir", outdir + "/peaks")
annotation_dir = config.get("annotation_dir", outdir + "/annotation")
qc_dir = config.get("qc_dir", outdir + "/QC/3_frip_score")
log_sample_dir = config.get("log_sample_dir", logdir)
markdup_dir = config.get("markdup_dir", outdir + "/common/4_markdup_bam")
trim_dir = outdir + "/common/2_trimmed_fastq"
metrics_dir = outdir + "/common/3_raw_bam"
te_dir = outdir + "/results/te_overlap"

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "PeakCalling_report", "bin", "generate_report.py")

all_samples = ip_samples + input_samples

rule generate_report:
    """
    Generate ChIP-seq Peak Calling Report (PPT + Excel).
    All content is parameterized via config — no hardcoded sample names or species.
    """
    input:
        # MACS3 peaks
        narrow_peaks = expand(peaks_dir + "/{sample}/{sample}_peaks.narrowPeak", sample=ip_samples),
        broad_peaks = expand(peaks_dir + "/{sample}/{sample}_broad_peaks.broadPeak", sample=ip_samples),
        xls_peaks = expand(peaks_dir + "/{sample}/{sample}_peaks.xls", sample=ip_samples),
        summits = expand(peaks_dir + "/{sample}/{sample}_summits.bed", sample=ip_samples),
        cutoff = expand(peaks_dir + "/{sample}/{sample}_cutoff_analysis.txt", sample=ip_samples),
        # FRiP QC
        frip = expand(qc_dir + "/{sample}/{sample}.FRiP.txt", sample=ip_samples),
        # HOMER annotation
        annotations = expand(annotation_dir + "/{sample}/{sample}_peaks.annotatePeaks.txt", sample=ip_samples),
        # Bowtie2 alignment logs
        bowtie2_logs = expand(log_sample_dir + "/{sample}/bowtie2_align.log", sample=all_samples),
        # Bowtie2 metrics
        bowtie2_metrics = expand(metrics_dir + "/{sample}/{sample}_bowtie2_metrics.txt", sample=all_samples),
        # MarkDuplicates
        markdup_metrics = expand(markdup_dir + "/{sample}/{sample}.Markdup-metrics.txt", sample=all_samples),
        # TrimGalore stats
        trim_stats_r1 = expand(trim_dir + "/{sample}/trimming_statistics_1.txt", sample=all_samples),
        trim_stats_r2 = expand(trim_dir + "/{sample}/trimming_statistics_2.txt", sample=all_samples),
        # TE overlap
        te_subfamily = expand(te_dir + "/{sample}/{sample}_te_subfamily_overlap.tsv", sample=ip_samples),
        peak_centric_te = expand(te_dir + "/{sample}/{sample}_peak_centric_te.tsv", sample=ip_samples),
        te_enrichment = expand(te_dir + "/{sample}_vs_{input}_enrichment.png",
                               zip, sample=ip_samples,
                               input=[sample_ip_input_map.get(s, "") for s in ip_samples]),
    output:
        report = outdir + "/PeakCalling_report.pptx",
        excel = outdir + "/PeakCalling_report.xlsx",
    log:
        logdir + "/PeakCalling_report.log"
    threads: 1
    conda:
        "PeakCalling_report.yaml"
    container:
        sif("PeakCalling_report.yaml")
    params:
        script = REPORT_SCRIPT,
        title = config.get("Params", {}).get("report", {}).get("title") or "ChIP-seq Peak Calling Report",
        subtitle = config.get("Params", {}).get("report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("date") or "",
        top_n = config.get("Params", {}).get("report", {}).get("top_n") or 5,
        sample_ip_input_map = sample_ip_input_map,
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("generate_report", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start generate_report at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)

            cmd = [
                "python3", params.script,
            ]
            for s in ip_samples:
                cmd += ["--samples", s]
            for s in input_samples:
                cmd += ["--input-samples", s]
            # IP-input mapping for enrichment figures
            for ip, inp in params.sample_ip_input_map.items():
                cmd += ["--ip-input-pair", f"{ip}:{inp}"]
            cmd += [
                "--peaks-dir", peaks_dir,
                "--annotation-dir", annotation_dir,
                "--qc-dir", qc_dir,
                "--log-dir", log_sample_dir,
                "--markdup-dir", markdup_dir,
                "--trim-dir", trim_dir,
                "--metrics-dir", metrics_dir,
                "--te-dir", te_dir,
                "--output", output.report,
                "--excel-output", output.excel,
                "--title", params.title,
                "--subtitle", params.subtitle,
                "--pipeline", params.pipeline,
                "--genome", params.genome,
                "--date", params.date,
                "--top-n", str(params.top_n),
            ]

            script = os.path.join(report_dir, f"generate_report_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during generate_report: {e}\n")
            logger.error(f"Error occurred during generate_report: {e}")
            raise e


rule report_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        report = outdir + "/PeakCalling_report.pptx",
        excel = outdir + "/PeakCalling_report.xlsx",
