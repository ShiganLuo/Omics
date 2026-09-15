include: "../common/common.smk"
import shlex
import json as _json
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
sample_groups = config.get("sample_groups", {})
short_names = config.get("short_names", {})
aligner = config.get("aligner", "star_3pass_gene")

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "ncRNAseq_report", "bin", "generate_report.py")

# ── Conditional inputs based on aligner ──────────────────────────────────
has_per_gene = (aligner == "star_3pass_gene")
has_3pass = aligner in ("star_3pass", "star_3pass_gene")
has_tailer = aligner in ("star_3pass", "star_3pass_gene")

def _report_inputs(wildcards):
    """Return only the inputs that exist for the current aligner branch."""
    inp = {
        "trimming_stats": expand(outdir + "/common/2_trimmed_dedup_fastq/final_trimmed_fastq/{sample}/trimming_statistics_1.txt", sample=samples),
    }
    if has_per_gene:
        inp["per_gene_bams"] = expand(outdir + "/common/4_per_gene_bam/{sample}/{sample}.bam", sample=samples)
        inp["per_gene_tails"] = expand(outdir + "/common/4_per_gene_bam/{sample}/{sample}_tail.csv", sample=samples)
    if has_tailer and not has_per_gene:
        inp["tailer_csvs"] = expand(outdir + "/results/tailer/{sample}/{sample}_tail.csv", sample=samples)
    if has_per_gene or has_3pass:
        inp["smallrna_bed"] = outdir + "/genome/smallrna/smallrna_genes.bed"
    return inp

rule generate_report:
    input:
        unpack(_report_inputs)
    output:
        report = outdir + "/ncRNAseq_report.pptx",
        file_inventory = outdir + "/ncRNAseq_report_files.xlsx",
    log:
        logdir + "/ncRNAseq_report.log"
    threads: 1
    conda:
        "ncRNAseq_report.yaml"
    container:
        sif("ncRNAseq_report.yaml")
    params:
        samples = samples,
        paired_samples = paired_samples,
        single_samples = single_samples,
        sample_groups = sample_groups,
        short_names = short_names,
        aligner = aligner,
        title = config.get("Params", {}).get("report", {}).get("title") or "ncRNAseq Analysis Report",
        subtitle = config.get("Params", {}).get("report", {}).get("subtitle") or "",
        pipeline = config.get("Params", {}).get("report", {}).get("pipeline") or "",
        genome = config.get("Params", {}).get("report", {}).get("genome") or "",
        date = config.get("Params", {}).get("report", {}).get("date") or "",
        lang = config.get("Params", {}).get("report", {}).get("lang") or "zh",
        img_dir = outdir + "/ppt_results",
        script = REPORT_SCRIPT,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("ncRNAseq_generate_report", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ncRNAseq report generation at {current_time}")
            report_dir = os.path.dirname(str(output.report))
            os.makedirs(report_dir, exist_ok=True)
            os.makedirs(params.img_dir, exist_ok=True)
            script_path = os.path.join(report_dir, f"ncRNAseq_report_{current_time}.sh")
            cmd = [
                "python3", params.script,
                "--analysis-dir", outdir,
                "--output", str(output.report),
                "--file-inventory", str(output.file_inventory),
                "--title", params.title,
                "--subtitle", params.subtitle,
                "--pipeline", params.pipeline,
                "--genome", params.genome,
                "--date", params.date,
                "--lang", params.lang,
                "--img-dir", params.img_dir,
                "--aligner", params.aligner,
            ]
            if params.samples:
                cmd.extend(["--samples", *params.samples])
            if params.paired_samples:
                cmd.extend(["--paired-samples", *params.paired_samples])
            if params.single_samples:
                cmd.extend(["--single-samples", *params.single_samples])
            if params.sample_groups:
                cmd.extend(["--sample-groups", _json.dumps(params.sample_groups)])
            if params.short_names:
                cmd.extend(["--short-names", _json.dumps(params.short_names)])
            with open(script_path, "w") as fh:
                fh.write("#!/bin/bash\n")
                fh.write("set -euo pipefail\n")
                fh.write(" ".join([shlex.quote(str(x)) for x in cmd]) + "\n")
            shell(f"bash {script_path} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as fh:
                fh.write(f"ncRNAseq report generation failed: {e}\n")
            logger.error(f"ncRNAseq report generation failed: {e}")
            raise


rule report_result:
    input:
        report = outdir + "/ncRNAseq_report.pptx",
        file_inventory = outdir + "/ncRNAseq_report_files.xlsx",
