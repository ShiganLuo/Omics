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
        svsig = outdir + "/{sample_id}/{sample_id}.svsig.gz"
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
        svsig = outdir + "/{sample_id}/{sample_id}.svsig.gz",
        fasta = fasta
    output:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.sv.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.sv.vcf.gz.csi"
    log:
        logdir + "/{sample_id}/pbsv_call.log"
    threads: 8
    conda:
        "pbsv.yaml"
    params:
        pbsv = config.get("Procedure", {}).get("pbsv") or "pbsv",
        bgzip = config.get("Procedure", {}).get("bgzip") or "bgzip",
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools",
        vcf = outdir + "/{sample_id}/{sample_id}.sv.vcf"
    run:
        current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        logger.info(f"Start pbsv call for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"pbsv_call_{current_time}.sh")
        cmd1 = [
            params.pbsv, "call",
            "--num-threads", str(threads),
            input.fasta,
            input.svsig,
            params.vcf
        ]
        cmd2 = [
            params.bgzip, params.vcf,
            "-o", output.vcf_gz
        ]
        cmd3 = [
            params.bcftools, "index", output.vcf_gz
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(" ".join(cmd3) + "\n")
        shell("bash {script} > {log} 2>&1")


rule pbsv_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.sv.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.sv.vcf.gz.csi"

