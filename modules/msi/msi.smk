include: "../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", logdir)
pairs = config.get("pairs", {})
samples = list(pairs)
bam_substring = config.get("bam_substring") or ""
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai")
msisensorpro = config.get("Procedure", {}).get("msisensorpro") or "msisensorpro"
BAR_PLOT = os.path.join(ROOT_DIR, "modules", "msi", "bin", "bar_plot.r")

def _bam(sample_id):
    if bam_substring != "":
        return os.path.join(indir, f"{sample_id}/{sample_id}." + bam_substring + ".bam")
    return os.path.join(indir, f"{sample_id}/{sample_id}.bam")

rule msi_scan:
    """Genome-wide microsatellite scan (msisensor-pro scan), run once per
    genome; every msi_call shares the scan table."""
    input:
        fasta = fasta,
        fai = fai,
    output:
        scan = outdir + "/msi_scan/msisensor_scan.tsv",
    log:
        logdir_combine + "/msi/msi_scan.log"
    threads: rule_res("msi_scan").get("threads", 1)
    resources:
        mem_mb = rule_res("msi_scan").get("mem_mb", 16000),
    params:
        msisensorpro = msisensorpro,
    conda:
        "msi.yaml"
    container:
        sif("msi.yaml", "msi_scan")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("msi_scan", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start msi scan at {current_time}")
            scan_dir = os.path.dirname(str(output.scan))
            os.makedirs(scan_dir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"msi_scan_{current_time}.sh")
            plain_script = os.path.join(scan_dir, "msi_scan.sh")
            cmd = [
                params.msisensorpro, "scan",
                "-d", str(input.fasta),
                "-o", str(output.scan),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "msi scan was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during msi scan: {e}\n")
            raise RuntimeError(f"Error occurred during msi scan: {e}, you can check the log file {log_path} for more details.\n")

rule msi_call:
    """MSI score per tumor sample against the matched normal
    (msisensor-pro msi)."""
    input:
        scan = outdir + "/msi_scan/msisensor_scan.tsv",
        tumor_bam = lambda wildcards: _bam(wildcards.sample_id),
        tumor_bai = lambda wildcards: _bam(wildcards.sample_id)[:-4] + ".bai",
        normal_bam = lambda wildcards: _bam(pairs[wildcards.sample_id]),
        normal_bai = lambda wildcards: _bam(pairs[wildcards.sample_id])[:-4] + ".bai",
    output:
        msi = outdir + "/{sample_id}/{sample_id}_msi_pro.tsv",
    log:
        logdir + "/{sample_id}/msi_call.log"
    threads: rule_res("msi_call").get("threads", 4)
    resources:
        mem_mb = rule_res("msi_call").get("mem_mb", 16000),
    params:
        msisensorpro = msisensorpro,
    conda:
        "msi.yaml"
    container:
        sif("msi.yaml", "msi_call")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("msi_call", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start msi call for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.msi))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"msi_call_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"msi_call_{wildcards.sample_id}.sh")
            cmd = [
                params.msisensorpro, "msi",
                "-d", str(input.scan),
                "-n", str(input.normal_bam),
                "-t", str(input.tumor_bam),
                "-o", str(output.msi),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "msi call for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during msi call for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during msi call for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

rule msi_plot:
    """Aggregate per-sample MSI scores into msi_input.txt and render the
    cross-sample bar plot (bin/bar_plot.r)."""
    input:
        msi = expand(outdir + "/{sample_id}/{sample_id}_msi_pro.tsv",
                     sample_id=samples),
    output:
        bar_input = outdir + "/msi_input.txt",
        barplot = outdir + "/msi_barplot.jpg",
    log:
        logdir_combine + "/msi/msi_plot.log"
    threads: rule_res("msi_plot").get("threads", 1)
    resources:
        mem_mb = rule_res("msi_plot").get("mem_mb", 4000),
    params:
        rscript = config.get("Procedure", {}).get("Rscript") or "Rscript",
        bar_plot = BAR_PLOT,
    conda:
        "msi.yaml"
    container:
        sif("msi.yaml", "msi_plot")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("msi_plot", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start msi plot at {current_time}")
            with open(str(output.bar_input), "w") as f:
                f.write("sample\tscore\n")
                for msi_tsv in input.msi:
                    sample = os.path.basename(msi_tsv).replace("_msi_pro.tsv", "")
                    with open(msi_tsv) as fh:
                        for i, line in enumerate(fh):
                            if i == 0 or not line.strip():
                                continue
                            f.write(f"{sample}\t{line.split()[2]}\n")
            script = os.path.join(RECEIPT_DIR, f"msi_plot_{current_time}.sh")
            plain_script = os.path.join(outdir, "msi_plot.sh")
            cmd = [
                params.rscript, params.bar_plot,
                "-i", str(output.bar_input),
                "-o", str(output.barplot),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "msi plot was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during msi plot: {e}\n")
            raise RuntimeError(f"Error occurred during msi plot: {e}, you can check the log file {log_path} for more details.\n")

rule msi_result:
    input:
        scan = outdir + "/msi_scan/msisensor_scan.tsv",
        msi = expand(outdir + "/{sample_id}/{sample_id}_msi_pro.tsv",
                     sample_id=samples),
        bar_input = outdir + "/msi_input.txt",
        barplot = outdir + "/msi_barplot.jpg",
