include: "../../common/common.smk"

shell.prefix("set -x; set -e;")
from snakemake.logging import logger
import os
import time

indir = config.get("indir", "data/fasta")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
fasta = config.get("genome", {}).get("fasta")

# Get parameters
decoy_params = config.get("Params", {}).get("decoy_database", {})
decoy_string = decoy_params.get("decoy_string", "DECOY_")
decoy_string_position = decoy_params.get("decoy_string_position", "prefix")
method = decoy_params.get("method", "shuffle")
shuffle_max_attempts = decoy_params.get("shuffle_max_attempts", 30)
shuffle_sequence_identity_threshold = decoy_params.get("shuffle_sequence_identity_threshold", 0.5)

# Get OpenMS executable
openms = config.get("Procedure", {}).get("openms") or "DecoyDatabase"

rule decoy_database:
    input:
        fasta = fasta
    output:
        decoy_fasta = outdir + "/" + os.path.basename(fasta) + "_decoy.fasta"
    log:
        logdir + "/decoy_database.log"
    conda:
        "openms.yaml"
    params:
        openms = openms,
        decoy_string = decoy_string,
        decoy_string_position = decoy_string_position,
        method = method,
        shuffle_max_attempts = shuffle_max_attempts,
        shuffle_sequence_identity_threshold = shuffle_sequence_identity_threshold
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start decoy database generation at {current_time}")
        script = os.path.join(outdir, f"decoy_database_{current_time}.sh")
        cmd = [
            params.openms,
            "-in", input.fasta,
            "-out", output.decoy_fasta,
            "-decoy_string", params.decoy_string,
            "-decoy_string_position", params.decoy_string_position,
            "-method", params.method,
            "-shuffle_max_attempts", str(params.shuffle_max_attempts),
            "-shuffle_sequence_identity_threshold", str(params.shuffle_sequence_identity_threshold)
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule decoy_database_result:
    input:
        decoy_fasta = outdir + "/" + os.path.basename(fasta) + "_decoy.fasta"
