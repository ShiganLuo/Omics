from snakemake.logging import logger
import time
import os
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai")
karyotype = config.get("Params", {}).get("trgt", {}).get("karyotype") or "XX"
repeat_bed = config.get("genome", {}).get("repeat_bed")
bam_substring = config.get("bam_substring") or ""

def get_input_for_trgt_genotype(wildcards):
    logger.info(f"trgt_genotype called with {wildcards}")
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam.bai")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam.bai")
    in_dict["fasta"] = fasta
    in_dict["fai"] = fai
    in_dict["bed"] = repeat_bed
    return in_dict
rule trgt_genotype:
    input:
        unpack(get_input_for_trgt_genotype)
    output:
        vcf = outdir + "/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        bam = outdir + "/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam"
    log:
        logdir + "/{sample_id}/trgt_genotype.log"
    threads: 4
    conda:
        "trgt.yaml"
    params:
        trgt = config.get("Procedure", {}).get("trgt") or "trgt",
        karyotype = karyotype,
        prefix = outdir + "/genotype/{sample_id}/{sample_id}.trgt"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start trgt genotype for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"{wildcards.sample_id}/trgt_genotype_{current_time}.sh")
        cmd = [params.trgt, "genotype",
            "--num-threads", str(threads),
            "--genome", input.fasta,
            "--repeats", input.bed,
            "--karyotype", params.karyotype,
            "--reads", input.bam,
            "--output-prefix", params.prefix
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell(f"bash {script} > {log} 2>&1")

rule trgt_plot:
    input:
        vcf = outdir + "/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        bam = outdir + "/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam",
        fasta = fasta,
        fai = fai,
        bed = repeat_bed
    output:
        png = outdir + "/plot/{sample_id}/{sample_id}.trgt.repeat.png"
    log:
        logdir + "/{sample_id}/trgt_plot.log"
    conda:
        "trgt.yaml"
    params:
        trgt = config.get("Procedure", {}).get("trgt") or "trgt",
        repeat_id = config.get("Params", {}).get("trgt", {}).get("repeat_id") or "HTT"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start trgt plot for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"{wildcards.sample_id}/trgt_plot_{current_time}.sh")
        cmd = [
            {params.trgt} , "plot",
            "--genome", input.fasta,
            "--repeats", input.bed,
            "--vcf", input.vcf,
            "--spanning-reads", input.bam,
            "--repeat-id", params.repeat_id,
            "--image", output.png
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell(f"bash {script} > {log} 2>&1")

rule trgt_result:
    input:
        vcf = outdir + "/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        png = outdir + "/plot/{sample_id}/{sample_id}.trgt.repeat.png"
