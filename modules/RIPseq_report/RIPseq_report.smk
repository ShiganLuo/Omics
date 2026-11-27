include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ip_samples = config.get("samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})

peaks_dir = config.get("peaks_dir", outdir + "/peaks")
annotation_dir = config.get("annotation_dir", outdir + "/annotation")
qc_dir = config.get("qc_dir", outdir + "/QC/3_frip_score")
log_sample_dir = config.get("log_sample_dir", logdir)
markdup_dir = config.get("markdup_dir", outdir + "/common/4_markdup_bam")
trim_dir = outdir + "/common/2_trimmed_fastq"
metrics_dir = outdir + "/common/3_raw_bam"
te_dir = outdir + "/results/te_overlap"
heatmap_dir = outdir + "/results/heatmap"


_formats = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_formats, str):
    _formats = [_formats]

_RASTER_FMT_ORDER = ("png", "jpg", "jpeg", "tif", "tiff")


def _producer_fmt(formats):
    """Raster-first extension among one producer's configured formats."""
    fmts = formats or ["png"]
    if isinstance(fmts, dict):
        fmts = fmts.get("default") or ["png"]
    if isinstance(fmts, str):
        fmts = [fmts]
    for fmt in _RASTER_FMT_ORDER:
        if fmt in fmts:
            return fmt
    return fmts[0]


_pto_formats = config.get("Params", {}).get("peak_te_overlap", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_pto_formats, str):
    _pto_formats = [_pto_formats]
_ph_formats = config.get("Params", {}).get("plotHeatmap", {}).get("generate_report", {}).get("formats") or ["png"]
if isinstance(_ph_formats, str):
    _ph_formats = [_ph_formats]

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "RIPseq_report", "bin", "generate_report.py")

all_samples = ip_samples + input_samples

rule generate_report:
    """
    Generate RIP-seq Report (PPT + Excel).
    All content is parameterized via config — no hardcoded sample names or species.
    """
    input:
        narrow_peaks = expand(peaks_dir + "/{sample}/{sample}_peaks.narrowPeak", sample=ip_samples),
        broad_peaks = expand(peaks_dir + "/{sample}/{sample}_broad_peaks.broadPeak", sample=ip_samples),
        xls_peaks = expand(peaks_dir + "/{sample}/{sample}_peaks.xls", sample=ip_samples),
        summits = expand(peaks_dir + "/{sample}/{sample}_summits.bed", sample=ip_samples),
        cutoff = expand(peaks_dir + "/{sample}/{sample}_cutoff_analysis.txt", sample=ip_samples),
        frip = expand(qc_dir + "/{sample}/{sample}.FRiP.txt", sample=ip_samples),
        annotations = expand(annotation_dir + "/{sample}/{sample}_peaks.annotatePeaks.txt", sample=ip_samples),
        hisat2_logs = expand(log_sample_dir + "/{sample}/hisat2_align.log", sample=all_samples),
        hisat2_anchors = expand(metrics_dir + "/{sample}/{sample}.splice.txt", sample=all_samples),
        markdup_metrics = expand(markdup_dir + "/{sample}/{sample}.Markdup-metrics.txt", sample=all_samples),
        trim_stats_r1 = expand(trim_dir + "/{sample}/trimming_statistics_1.txt", sample=all_samples),
        trim_stats_r2 = expand(trim_dir + "/{sample}/trimming_statistics_2.txt", sample=all_samples),
        te_subfamily = expand(te_dir + "/{sample}/{sample}_te_subfamily_overlap.tsv", sample=ip_samples),
        peak_centric_te = expand(te_dir + "/{sample}/{sample}_peak_centric_te.tsv", sample=ip_samples),
        te_enrichment = expand(te_dir + "/{sample}_vs_{input}_enrichment." + _producer_fmt(_pto_formats),
                               zip, sample=ip_samples,
                               input=[sample_ip_input_map.get(s, "") for s in ip_samples]),
        heatmaps = expand(heatmap_dir + "/{sample}/{sample}_{te}_heatmap." + _producer_fmt(_ph_formats),
                          sample=ip_samples,
                          te=config.get("Params", {}).get("computeMatrix", {}).get("regions", {}).get("genes", []) if isinstance(config.get("Params", {}).get("computeMatrix", {}).get("regions"), dict) else []),
    output:
        report = outdir + "/RIPseq_report.pptx",
        excel = outdir + "/RIPseq_report.xlsx",
        figures_manifest = outdir + "/figures/manifest.tsv",
    log:
        logdir + "/RIPseq_report.log"
    threads: rule_res("generate_report").get("threads", 1)
    benchmark: logdir + "/RIPseq_report.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "generate_report", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "RIPseq_report.yaml"
    container:
        sif("RIPseq_report.yaml", "generate_report")
    params:
        script = REPORT_SCRIPT,
        title = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("title") or "RIP-seq Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("date") or "",
        top_n = config.get("Params", {}).get("report", {}).get("generate_report", {}).get("top_n") or 5,
        sample_ip_input_map = sample_ip_input_map,
        formats = _formats,
        img_dir = outdir + "/figures",
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("generate_report", log_file=log_path)
        try:
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
                "--heatmap-dir", heatmap_dir,
                "--img-dir", params.img_dir,
                "--output", output.report,
                "--excel-output", output.excel,
                "--title", params.title,
                "--subtitle", params.subtitle,
                "--pipeline", params.pipeline,
                "--genome", params.genome,
                "--date", params.date,
                "--top-n", str(params.top_n),
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]

            script = receipt_script(f"generate_report_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(report_dir, f"generate_report.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during generate_report: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule report_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        report = outdir + "/RIPseq_report.pptx",
        excel = outdir + "/RIPseq_report.xlsx",
        figures_manifest = outdir + "/figures/manifest.tsv",
