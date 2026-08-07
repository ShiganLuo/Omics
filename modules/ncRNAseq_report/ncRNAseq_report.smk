include: "../common/common.smk"
import shlex
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])

REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "ncRNAseq_report", "bin", "generate_report.py")

rule generate_report:
    input:
        per_gene_bams = expand(outdir + "/common/4_per_gene_bam/{sample}/{sample}.bam", sample=samples),
        per_gene_tails = expand(outdir + "/common/4_per_gene_bam/{sample}/{sample}_tail.csv", sample=samples),
        per_gene_manifests = expand(outdir + "/common/4_per_gene_bam/{sample}/genes.tsv", sample=samples),
        per_gene_overlaps = expand(outdir + "/common/4_per_gene_bam/{sample}/read_gene_overlaps.tsv", sample=samples),
        star_logs = expand(outdir + "/common/3_raw_bam/{sample}/star.Log.final.out", sample=samples),
        trimming_stats = expand(outdir + "/common/2_trimmed_dedup_fastq/final_trimmed_fastq/{sample}/trimming_statistics_1.txt", sample=samples),
        smallrna_bed = outdir + "/genome/smallrna/smallrna_genes.bed",
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
            ]
            if params.samples:
                cmd.extend(["--samples", *params.samples])
            if params.paired_samples:
                cmd.extend(["--paired-samples", *params.paired_samples])
            if params.single_samples:
                cmd.extend(["--single-samples", *params.single_samples])
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
