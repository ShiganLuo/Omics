shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
import time

indir = config.get("indir", "data/quantification")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
quantification_method = config.get("quantification_method", "lfq")
skip_post_msstats = config.get("Params", {}).get("skip_post_msstats", False)

# Get parameters
msstats_params = config.get("Params", {}).get("msstats", {})
reference_condition = msstats_params.get("reference_condition", "")
reference_samples = msstats_params.get("reference_samples", [])

# Get executables
msstats = config.get("Procedure", {}).get("msstats") or "MSstatsConverter"

rule msstats:
    input:
        mztab = indir + "/{quant_method}_quantification.mzTab".format(quant_method=quantification_method)
    output:
        csv = outdir + "/msstats_results.csv"
    log:
        logdir + "/msstats.log"
    conda:
        "openms.yaml"
    params:
        msstats = msstats,
        reference_condition = reference_condition,
        reference_samples = reference_samples
    run:
        if skip_post_msstats:
            logger.info("Skipping MSstats analysis")
            # Create empty output file
            with open(output.csv, "w") as f:
                f.write("MSstats analysis skipped\n")
            return
        
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start MSstats analysis at {current_time}")
        script = os.path.join(outdir, f"msstats_{current_time}.sh")
        cmd = [
            params.msstats,
            "-in", input.mztab,
            "-out", output.csv
        ]
        if params.reference_condition:
            cmd.extend(["-reference_condition", params.reference_condition])
        if params.reference_samples:
            cmd.extend(["-reference_samples", ",".join(params.reference_samples)])
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule msstats_result:
    input:
        outdir + "/msstats_results.csv"
