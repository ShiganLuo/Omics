from snakemake.logging import logger
import time
import os
import subprocess

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""
species = config.get("Params", {}).get("telogator2", {}).get("species", "human")

def get_telogator2_ref():
    """Get species-specific subtelomere reference path."""
    if species == "human":
        return ""  # use default
    # Find telogator2 installation and use non-human ref
    try:
        result = subprocess.run(["which", "telogator2"], capture_output=True, text=True)
        telogator2_bin = result.stdout.strip()
        # Resources are in site-packages/source/resources/
        resource_dir = os.path.join(os.path.dirname(telogator2_bin), "..",
            "lib", f"python{'.'.join(map(str, __import__('sys').version_info[:2]))}",
            "site-packages", "source", "resources")
        ref_file = os.path.join(resource_dir, "non-human", f"telogator-ref-{species}.fa.gz")
        if os.path.exists(ref_file):
            return ref_file
        logger.warning(f"No telogator2 reference found for species '{species}', using default (human)")
    except Exception as e:
        logger.warning(f"Error finding telogator2 ref: {e}")
    return ""

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
        telogator2 = config.get("Procedure", {}).get("telogator2") or "telogator2",
        ref = get_telogator2_ref()
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
        if params.ref:
            cmd.extend(["-t", params.ref])
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell(f"bash {script} > {log} 2>&1")



rule telomere_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        tsv = outdir + "/{sample_id}/telomere_lengths.tsv"
