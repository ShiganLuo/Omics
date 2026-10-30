include: "../../common/common.smk"

from snakemake.logging import logger
import os
import shutil
import time

indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/annotate")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
p = config.get("Params", {}).get("fiberseq_annotate", {})

def get_input_for_fiberseq_annotate_run(wildcards):
    """Get input files for fiberseq_annotate_run rule."""
    peaks = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.fire_peaks.bed"
    gtf = genome_cfg.get("gtf")
    if not gtf:
        raise ValueError("Genome GTF file is not specified in the configuration.")
    repeat_bed = genome_cfg.get("repeat_bed")
    if not repeat_bed or not os.path.exists(repeat_bed):
        logger.warning(f"Repeat BED file is not specified or does not exist: {repeat_bed}. Repeat overlap will be skipped.")
        repeat_bed = None
    return {"peaks": peaks, "gtf": gtf, "repeat_bed": repeat_bed}

rule fiberseq_annotate_run:
    """Annotate FIRE peaks with gene context and repeat overlap."""
    input:
        unpack(get_input_for_fiberseq_annotate_run)
    output:
        annotated = outdir + "/{sample_id}/{sample_id}.peaks_annotated.tsv",
        summary = outdir + "/{sample_id}/{sample_id}.peaks_summary.tsv",
    log:
        logdir + "/{sample_id}/fiberseq_annotate_run.log"
    threads: rule_res("fiberseq_annotate").get("threads", 4)
    resources:
        mem_mb = rule_res("fiberseq_annotate").get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/fire_annotate.py",
        tss_flank = p.get("tss_flank", 1000),
        repeat_frac = p.get("repeat_frac", 0.25),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_annotate_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start annotate for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.annotated))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(cmd_dir(logdir, str(wildcards.sample_id)), f"fiberseq_annotate_run_{wildcards.sample_id}_{current_time}.sh")
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_annotate_run_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-p", str(input.peaks),
                "-g", str(input.gtf),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-t", str(params.tss_flank),
                "-f", str(params.repeat_frac),
                "-l", log_path,
            ]
            if input.repeat_bed:
                cmd += ["-r", str(input.repeat_bed)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_annotate_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_annotate_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_annotate_run for {wildcards.sample_id}: {e}\n")


rule fiberseq_annotate_result:
    input:
        annotated = expand(outdir + "/{sample_id}/{sample_id}.peaks_annotated.tsv",
                           sample_id=samples),
        summary = expand(outdir + "/{sample_id}/{sample_id}.peaks_summary.tsv",
                         sample_id=samples),
