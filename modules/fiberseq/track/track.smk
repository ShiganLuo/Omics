include: "../../common/common.smk"

from snakemake.logging import logger
import os
import time

indir = config.get("indir", "results")
outdir = config.get("outdir", "results/fiberseq_track")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
p = config.get("Params", {}).get("fiberseq_track", {})


rule fiberseq_track_run:
    """bedGraph/bigWig accessibility and m6A density tracks + trackDb."""
    input:
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
        chrom_sizes = lambda wildcards: _require(genome_cfg.get("fai"), "genome.fai", "fiberseq_track"),
    output:
        acc_bw = outdir + "/{sample_id}/{sample_id}.accessibility.bw",
        m6a_bw = outdir + "/{sample_id}/{sample_id}.m6a_density.bw",
        trackdb = outdir + "/{sample_id}/trackDb.txt",
    log:
        logdir + "/{sample_id}/fiberseq_track_run.log"
    threads: 1
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/fiberseq_track.py",
        bin = p.get("bin", 10),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_track_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start track generation for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.acc_bw))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"fiberseq_track_run_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "-s", str(input.msp),
                "-a", str(input.m6a),
                "-u", str(input.nuc),
                "-c", str(input.chrom_sizes),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-b", str(params.bin),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_track_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_track_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_track_run for {wildcards.sample_id}: {e}\n")


rule fiberseq_track_result:
    input:
        tracks = expand(outdir + "/{sample_id}/{sample_id}.accessibility.bw",
                        sample_id=samples),
