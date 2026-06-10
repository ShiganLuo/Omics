from snakemake.logging import logger
import time
import os
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")

rule pbsv_discover:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        bai = indir + "/{sample_id}/{sample_id}.bam.bai"
    output:
        svsig = outdir + "/discover/{sample_id}/{sample_id}.svsig.gz"
    log:
        logdir + "/{sample_id}/pbsv_discover.log"
    threads: 4
    conda:
        "pbsv.yaml"
    params:
        pbsv = config.get("Procedure", {}).get("pbsv") or "pbsv"
    run:
        current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        logger.info(f"Start pbsv discover for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"pbsv_discover_{current_time}.sh")
        cmd = [
            params.pbsv, "discover",
            "--num-threads", str(threads),
            input.bam,
            output.svsig
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")


rule pbsv_call:
    input:
        svsig = outdir + "/discover/{sample_id}/{sample_id}.svsig.gz",
        fasta = fasta
    output:
        vcf = outdir + "/call/{sample_id}/{sample_id}.vcf"
    log:
        logdir + "/{sample_id}/pbsv_call.log"
    threads: 8
    conda:
        "pbsv.yaml"
    params:
        pbsv = config.get("Procedure", {}).get("pbsv") or "pbsv"
    run:
        current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        logger.info(f"Start pbsv call for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"pbsv_call_{current_time}.sh")
        cmd = [
            params.pbsv, "call",
            "--num-threads", str(threads),
            input.fasta,
            input.svsig,
            output.vcf
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")


rule pbsv_result:
    input:
        vcf = outdir + "/call/{sample_id}/{sample_id}.vcf"
