include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_footprint")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
p = config.get("Params", {}).get("fiberseq_footprint", {})

# upstream contract: motif sites come from modules/fiberseq_motif
MOTIF_DIR = config.get("motif_dir") or indir + "/fiberseq_motif"


rule fiberseq_footprint_run:
    """Single-fiber TF footprinting (3-category) around motif sites."""
    input:
        sites = MOTIF_DIR + "/{sample_id}/{sample_id}.motif_sites.bed",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        categories = outdir + "/{sample_id}/{sample_id}.footprint_categories.tsv",
        profiles = outdir + "/{sample_id}/{sample_id}.footprint_profiles.tsv",
        scores = outdir + "/{sample_id}/{sample_id}.footprint_scores.tsv",
        pdf = outdir + "/{sample_id}/{sample_id}.footprint_profile.pdf",
    log:
        logdir + "/{sample_id}/fiberseq_footprint_run.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/footprint_profile.py",
        flank = p.get("flank", 1000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_footprint_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start footprint for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.categories))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_footprint_run_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.sites),
                "-a", str(input.m6a),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-f", str(params.flank),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_footprint_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_footprint_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_footprint_run for {wildcards.sample_id}: {e}\n")


rule fiberseq_footprint_result:
    input:
        categories = expand(outdir + "/{sample_id}/{sample_id}.footprint_categories.tsv",
                            sample_id=samples),
        scores = expand(outdir + "/{sample_id}/{sample_id}.footprint_scores.tsv",
                        sample_id=samples),
