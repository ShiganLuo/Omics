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
percolator_params = config.get("Params", {}).get("psm_rescoring", {}).get("percolator", {})

search_engine = config.get("Params", {}).get("search_engine")
# Get executables
percolator = config.get("Procedure", {}).get("percolator") or "PercolatorAdapter"
def get_input_for_psm_rescoring(wildcards):
    """Get input files for PSM rescoring."""
    if search_engine == "comet":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_comet.idXML"}
    elif search_engine == "msgf":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_msgf.idXML"}
    elif search_engine == "sage":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_sage.idXML"}
    else:
        raise ValueError("No search engine specified in the configuration.")

rule psm_rescoring:
    input:
        unpack(get_input_for_psm_rescoring)
    output:
        scored_idxml = outdir + "/{sample_id}/{sample_id}_scored.idXML"
    log:
        logdir + "/{sample_id}/psm_rescoring.log"
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml")
    params:
        percolator = percolator,
        train_FDR = percolator_params.get("train_FDR", 0.05),
        test_FDR = percolator_params.get("test_FDR", 0.05),
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("psm_rescoring", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start PSM rescoring for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/psm_rescoring_{current_time}.sh")
            # OpenMS 3.x PercolatorAdapter uses -trainFDR/-testFDR (no
            # underscores); the legacy 2.x -feature option no longer exists.
            cmd = [
                params.percolator,
                "-in", input.idxml,
                "-out", output.scored_idxml,
                "-trainFDR", str(params.train_FDR),
                "-testFDR", str(params.test_FDR),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'PSM rescoring for {wildcards.sample_id} was completed successfully'\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule psm_rescoring was call failed,error: {e}")
            raise RuntimeError(f"rule psm_rescoring was call failed,error: {e}")

rule psm_rescoring_result:
    input:
        expand(outdir + "/{sid}/{sid}_scored.idXML", sid=samples)
