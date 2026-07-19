include: "../../common/common.smk"

shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
import time

indir = config.get("indir", "data/psm_fdr")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
protein_inference_params = config.get("Params", {}).get("protein_inference", {})
method = protein_inference_params.get("method", "epifany")
greedy_group_resolution = protein_inference_params.get("greedy_group_resolution", "none")

# Get executables
epifany = config.get("Procedure", {}).get("epifany") or "Epifany"

rule protein_inference:
    input:
        idxml = indir + "/{sample_id}/{sample_id}_filtered.idXML"
    output:
        protein_idxml = outdir + "/{sample_id}/{sample_id}_protein.idXML"
    log:
        logdir + "/{sample_id}/protein_inference.log"
    conda:
        "openms.yaml"
    params:
        epifany = epifany,
        method = method,
        greedy_group_resolution = greedy_group_resolution
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start protein inference for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/protein_inference_{current_time}.sh")
        cmd = [
            params.epifany,
            "-in", input.idxml,
            "-out", output.protein_idxml,
            "-algorithm:method", params.method,
            "-algorithm:greedy_group_resolution", params.greedy_group_resolution
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule protein_inference_result:
    input:
        expand(outdir + "/{sid}/{sid}_protein.idXML", sid=samples)
