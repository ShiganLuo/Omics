include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_coactuation")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
p = config.get("Params", {}).get("fiberseq_coactuation", {})


rule fiberseq_coactuation_run:
    """Same-fiber co-actuation of region pairs (Fisher + actual vs expected)."""
    input:
        peaks = indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        fire = indir + "/{sample_id}/{sample_id}.fire.bed",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        pairs = outdir + "/{sample_id}/{sample_id}.coactuation.tsv",
        summary = outdir + "/{sample_id}/{sample_id}.coactuation_summary.tsv",
    log:
        logdir + "/{sample_id}/fiberseq_coactuation_run.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/coactuation.py",
        max_distance = p.get("max_distance", 50000),
        min_fibers = p.get("min_fibers", 10),
        max_pairs = p.get("max_pairs", 200000),
        alpha = p.get("alpha", 0.05),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_coactuation_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start co-actuation for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.pairs))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_coactuation_run_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-p", str(input.peaks),
                "-f", str(input.fire),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-d", str(params.max_distance),
                "-i", str(params.min_fibers),
                "-m", str(params.max_pairs),
                "-a", str(params.alpha),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_coactuation_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_coactuation_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_coactuation_run for {wildcards.sample_id}: {e}\n")


rule fiberseq_coactuation_result:
    input:
        pairs = expand(outdir + "/{sample_id}/{sample_id}.coactuation.tsv",
                       sample_id=samples),
