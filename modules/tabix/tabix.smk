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
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start bgzip for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"bgzip_{current_time}.sh")
        cmd1 = [
            params.bgzip, input.vcf,
            "-o", output.vcf_gz
        ]
        cmd2 = [
            params.tabix, "-p", "vcf", output.vcf_gz
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
        shell("bash {script} > {log} 2>&1")

rule tabix_result:
    input:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
