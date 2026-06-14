from snakemake.logging import logger
import time
import os

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""

def get_input_for_telogator2(wildcards):
    logger.info(f"telogator2_run called with {wildcards}")
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bai")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bai")
    return in_dict
rule telogator2_run:
    """Run Telogator2 telomere length analysis on a single sample."""
    input:
        unpack(get_input_for_telogator2)
    output:
        tsv = outdir + "/{sample_id}/telomere_lengths.tsv",
        dir = directory(outdir + "/{sample_id}")
    log:
        logdir + "/{sample_id}/telogator2.log"
    params:
        telogator2 = config.get("Procedure", {}).get("telogator2") or "telogator2"
    threads: 16
    conda:
        "telomere.yaml"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start telogator2 for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"{wildcards.sample_id}/telogator2_{current_time}.sh")
        cmd = [
            "telogator2", "-i", input.bam,
            "-o", output.dir,
            "-r", "hifi",
            "-p", str(threads)
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell(f"bash {script} > {log} 2>&1")



rule telomere_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        tsv = outdir + "/{sample_id}/telomere_lengths.tsv"
