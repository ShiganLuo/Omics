from snakemake.logging import logger
import time
import os
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
bam_dir = config.get("bam_dir", indir)
vcf_dir = config.get("vcf_dir", indir)

rule hiphase_phase:
    input:
        bam = bam_dir + "/{sample_id}/{sample_id}.bam",
        bai = bam_dir + "/{sample_id}/{sample_id}.bam.bai",
        vcf = vcf_dir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = vcf_dir + "/{sample_id}/{sample_id}.vcf.gz.tbi",
        fasta = fasta
    output:
        vcf = outdir + "/{sample_id}/{sample_id}.phased.vcf.gz",
        bam = outdir + "/{sample_id}/{sample_id}.phased.bam"
    log:
        logdir + "/{sample_id}/hiphase.log"
    threads: 8
    conda:
        "hiphase.yaml"
    params:
        hiphase = config.get("Procedure", {}).get("hiphase") or "hiphase"
    run:
        current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        logger.info(f"Start hiphase for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"hiphase_{current_time}.sh")
        cmd = [
            params.hiphase,
            "--num-threads", str(threads),
            "--reference", input.fasta,
            "--input-bam", input.bam,
            "--output-bam", output.bam,
            "--input-vcf", input.vcf,
            "--output-vcf", output.vcf
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

rule hiphase_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.phased.vcf.gz",
        bam = outdir + "/{sample_id}/{sample_id}.phased.bam"
