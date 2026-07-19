include: "../common/common.smk"

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
input_vcf_substring = config.get("input_vcf_substring") or ""
input_bam_substring = config.get("input_bam_substring") or ""
output_substring = config.get("output_substring") or ""
def get_input_for_hiphase_phase(wildcards):
    logger.info(f"hiphase_phase called with {wildcards}")
    in_dict = {}
    if input_bam_substring != "":
        in_dict["bam"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_bam_substring}.bam")
        in_dict["bai"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_bam_substring}.bai")
    else:
        in_dict["bam"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(bam_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam.bai")
    if input_vcf_substring != "":
        in_dict["vcf"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_vcf_substring}.vcf.gz")
        in_dict["csi"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.{input_vcf_substring}.vcf.gz.csi")
    else:
        in_dict["vcf"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.vcf.gz")
        in_dict["csi"] = os.path.join(vcf_dir, f"{wildcards.sample_id}/{wildcards.sample_id}.vcf.gz.csi")
    in_dict["fasta"] = fasta
    logger.info(f"hiphase_phase input for sample {wildcards.sample_id}: {in_dict}")
    return in_dict

rule hiphase_phase:
    input:
        unpack(get_input_for_hiphase_phase)
    output:
        vcf = outdir + "/{sample_id}/{sample_id}.phased.vcf.gz" if output_substring == "" else outdir + "/{sample_id}/{sample_id}." + output_substring + ".phased.vcf.gz",
        bam = outdir + "/{sample_id}/{sample_id}.phased.bam" if output_substring == "" else outdir + "/{sample_id}/{sample_id}." + output_substring + ".phased.bam"
    log:
        logdir + "/{sample_id}/hiphase.log"
    threads: 8
    params:
        hiphase = config.get("Procedure", {}).get("hiphase") or "hiphase",
    conda:
        "hiphase.yaml"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start hiphase for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"{wildcards.sample_id}/hiphase_{current_time}.sh")
        cmd1 = [
            params.hiphase,
            "--threads", str(threads),
            "--reference", input.fasta,
            "--bam", input.bam,
            "--output-bam", output.bam,
            "--vcf", input.vcf,
            "--output-vcf", output.vcf
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
        shell("bash {script} > {log} 2>&1")

rule hiphase_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.phased.vcf.gz" if output_substring == "" else outdir + "/{sample_id}/{sample_id}." + output_substring + ".phased.vcf.gz",
        bam = outdir + "/{sample_id}/{sample_id}.phased.bam" if output_substring == "" else outdir + "/{sample_id}/{sample_id}." + output_substring + ".phased.bam"
