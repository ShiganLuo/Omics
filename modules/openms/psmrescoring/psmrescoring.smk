include: "../../common/common.smk"

shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
import time

indir = config.get("indir", "data/search_engine")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
psm_rescoring_params = config.get("Params", {}).get("psm_rescoring", {})
percolator_params = psm_rescoring_params.get("percolator", {})

# Get executables
percolator = config.get("Procedure", {}).get("percolator") or "PercolatorAdapter"

rule psm_rescoring:
    input:
        idxml = indir + "/{sample_id}/{sample_id}_comet.idXML"
    output:
        scored_idxml = outdir + "/{sample_id}/{sample_id}_scored.idXML"
    log:
        logdir + "/{sample_id}/psm_rescoring.log"
    conda:
        "openms.yaml"
    params:
        percolator = percolator,
        train_FDR = percolator_params.get("train_FDR", 0.05),
        test_FDR = percolator_params.get("test_FDR", 0.05),
        feature = percolator_params.get("feature", "top_psm")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start PSM rescoring for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/psm_rescoring_{current_time}.sh")
        cmd = [
            params.percolator,
            "-in", input.idxml,
            "-out", output.scored_idxml,
            "-train_FDR", str(params.train_FDR),
            "-test_FDR", str(params.test_FDR),
            "-feature", params.feature
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule psm_rescoring_result:
    input:
        expand(outdir + "/{sid}/{sid}_scored.idXML", sid=samples)
