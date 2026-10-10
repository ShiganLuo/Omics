include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_motif")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
p = config.get("Params", {}).get("fiberseq_motif", {})


rule fiberseq_motif_scan:
    """Scan a MEME motif file in FIRE peak sequences (any species)."""
    input:
        peaks = indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        fasta = lambda wildcards: _require(genome_cfg.get("fasta"), "genome.fasta", "fiberseq_motif"),
        motif_file = lambda wildcards: _require(genome_cfg.get("motif_file"), "genome.motif_file", "fiberseq_motif"),
    output:
        sites = outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
    log:
        logdir + "/{sample_id}/fiberseq_motif_scan.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/motif_scan.py",
        pwm_threshold = p.get("pwm_threshold", 0.8),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_motif_scan", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start motif scan for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.sites))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_motif_scan_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.motif_file),
                "-r", str(input.peaks),
                "-f", str(input.fasta),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-t", str(params.pwm_threshold),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_motif_scan for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_motif_scan for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_motif_scan for {wildcards.sample_id}: {e}\n")


rule fiberseq_motif_profile:
    """Aggregate m6A/accessibility/nucleosome profiles around motif sites."""
    input:
        sites = outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        profile = outdir + "/{sample_id}/{sample_id}.motif_profile.tsv",
        scores = outdir + "/{sample_id}/{sample_id}.motif_footprint_scores.tsv",
        pdf = outdir + "/{sample_id}/{sample_id}.motif_profile.pdf",
    log:
        logdir + "/{sample_id}/fiberseq_motif_profile.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/motif_profile.py",
        flank = p.get("flank", 250),
        min_fibers = p.get("min_fibers", 10),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_motif_profile", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start motif profile for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.profile))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_motif_profile_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.sites),
                "-a", str(input.m6a),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-f", str(params.flank),
                "-i", str(params.min_fibers),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_motif_profile for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_motif_profile for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_motif_profile for {wildcards.sample_id}: {e}\n")


rule fiberseq_motif_result:
    input:
        sites = expand(outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
                       sample_id=samples),
        profile = expand(outdir + "/{sample_id}/{sample_id}.motif_profile.tsv",
                         sample_id=samples),
