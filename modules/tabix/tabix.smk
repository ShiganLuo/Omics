include: "../common/common.smk"

from snakemake.logging import logger
import time
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule tabix_bgzip:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf"
    output:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
    log:
        logdir + "/{sample_id}/bgzip.log"
    params:
        bgzip = config.get("Procedure", {}).get("bgzip") or "bgzip"
    conda:
        "tabix.yaml"
    container:
        sif("tabix.yaml")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start bgzip for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(RECEIPT_DIR, f"bgzip_{wildcards.sample_id}_{current_time}.sh")
        plain_script = os.path.join(outdir, f"bgzip.sh")
        cmd1 = [
            params.bgzip, input.vcf,
            "-o", output.vcf_gz
        ]
        cmd2 = [
            params.tabix, "-p", "vcf", output.vcf_gz
        ]
        with open(script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(f'echo "script completed successfully"\n')
        shutil.copyfile(script, plain_script)
        shell("bash {plain_script} > {log} 2>&1")

rule tabix_result:
    input:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
