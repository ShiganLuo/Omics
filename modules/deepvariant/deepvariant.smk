from snakemake.logging import logger
import time
import os
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai")
bam_substring = config.get("bam_substring") or ""

rule deepvariant_run:
    input:
        bam = indir + "/{sample_id}/{sample_id}." + bam_substring + ".bam",
        bai = indir + "/{sample_id}/{sample_id}." + bam_substring + ".bam.bai",
        fasta = fasta,
        fai = fai
    output:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi",
        gvcf = outdir + "/{sample_id}/{sample_id}.g.vcf.gz"
    log:
        logdir + "/{sample_id}/deepvariant.log"
    threads: 8
    container:
        "docker://google/deepvariant:latest"
    params:
        deepvariant = config.get("Procedure", {}).get("deepvariant") or "run_deepvariant",
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools",
        model_type = config.get("Params", {}).get("deepvariant", {}).get("model_type") or "PACBIO",
        outdir_sample = outdir + "/{sample_id}"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start DeepVariant for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"deepvariant_{current_time}.sh")
        cmd1 = [
            params.deepvariant,
            "--num_shards", str(threads),
            "--model_type", params.model_type,
            "--ref", input.fasta,
            "--reads", input.bam,
            "--output_vcf", output.vcf,
            "--output_gvcf", output.gvcf
        ]
        cmd2 = [
            params.bcftools, "index", output.vcf
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
        shell("bash {script} > {log} 2>&1")

rule deepvariant_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
