include: "../../common/common.smk"

shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
import time

indir = config.get("indir", "data/psm_rescoring")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
psm_fdr_params = config.get("Params", {}).get("psm_fdr_control", {})
fdr = psm_fdr_params.get("fdr", 0.01)
method = psm_fdr_params.get("method", "BH")

# Get executables
openms = config.get("Procedure", {}).get("openms") or "FalseDiscoveryRate"

rule psm_fdr:
    input:
        idxml = indir + "/{sample_id}/{sample_id}_scored.idXML"
    output:
        filtered_idxml = outdir + "/{sample_id}/{sample_id}_filtered.idXML"
    log:
        logdir + "/{sample_id}/psm_fdr.log"
    conda:
        "openms.yaml"
    params:
        openms = openms,
        fdr = fdr,
        method = method
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start PSM FDR control for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/psm_fdr_{current_time}.sh")
        cmd = [
            params.openms,
            "-in", input.idxml,
            "-out", output.filtered_idxml,
            "-FDR:PSM", str(params.fdr),
            "-method", params.method
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule psm_fdr_result:
    input:
        expand(outdir + "/{sid}/{sid}_filtered.idXML", sid=samples)
