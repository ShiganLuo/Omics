from snakemake.logging import logger
include: "../common/common.smk"
import time
import os

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", samples)

def get_fumitools_input(wildcards):
    """Dynamically determine paired-end or single-end input."""
    sid = wildcards.sample_id
    if sid in paired_samples:
        r1 = f"{indir}/{sid}/{sid}_1.fq.gz"
        r2 = f"{indir}/{sid}/{sid}_2.fq.gz"
        return {"r1": r1, "r2": r2}
    return {"r1": f"{indir}/{sid}/{sid}.fq.gz"}

# ========================
# copy_umi — extract UMI from read into header
# ========================

rule fumitools_copy_umi:
    """Copy UMI from read sequence into FASTQ header."""
    input:
        unpack(get_fumitools_input)
    output:
        r1 = outdir + "/{sample_id}/{sample_id}_1.umi.fq.gz",
    log:
        logdir + "/{sample_id}/fumitools_copy_umi.log"
    conda:
        "fumitools.yaml"
    params:
        fumi_tools = config.get("Procedure", {}).get("fumitools") or "fumi_tools",
        umi_length = config.get("Params", {}).get("fumitools", {}).get("umi_length", 12),
        tag_umi = config.get("Params", {}).get("fumitools", {}).get("tag_umi", False),
    threads: config.get("Params", {}).get("fumitools", {}).get("threads", 4)
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"{wildcards.sample_id}/fumitools_copy_umi_{current_time}.sh")

            cmd = [
                params.fumi_tools, "copy_umi",
                "-i", input.r1,
                "-o", output.r1,
                "--umi-length", str(params.umi_length),
                "--threads", str(threads),
            ]
            if params.tag_umi:
                cmd += ["--tag-umi"]

            # Paired-end: also process R2
            if "r2" in input.keys():
                r2_out = outdir + f"/{wildcards.sample_id}/{wildcards.sample_id}_2.umi.fq.gz"
                cmd += ["-I", input.r2, "-O", r2_out]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"fumitools_copy_umi failed for {wildcards.sample_id}: {e}\n")
            raise

# ========================
# dedup — UMI-based deduplication on sorted BAM
# ========================

rule fumitools_dedup:
    """Deduplicate BAM reads based on UMI and mapping position."""
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted.bam",
    output:
        bam = outdir + "/{sample_id}/{sample_id}.dedup.bam",
    log:
        logdir + "/{sample_id}/fumitools_dedup.log"
    conda:
        "fumitools.yaml"
    params:
        fumi_tools = config.get("Procedure", {}).get("fumitools") or "fumi_tools",
        memory = config.get("Params", {}).get("fumitools", {}).get("memory", "3G"),
        start_only = config.get("Params", {}).get("fumitools", {}).get("start_only", False),
    threads: config.get("Params", {}).get("fumitools", {}).get("threads", 4)
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script_path = os.path.join(outdir, f"{wildcards.sample_id}/fumitools_dedup_{current_time}.sh")

            cmd = [
                params.fumi_tools, "dedup",
                "-i", input.bam,
                "-o", output.bam,
                "--threads", str(threads),
                "--memory", params.memory,
            ]
            if wildcards.sample_id in paired_samples:
                cmd += ["--paired"]
            if params.start_only:
                cmd += ["--start-only"]

            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"fumitools_dedup failed for {wildcards.sample_id}: {e}\n")
            raise

# ========================
# Result rules
# ========================

rule fumitools_copy_umi_result:
    input:
        r1 = outdir + "/{sample_id}/{sample_id}_1.umi.fq.gz",

rule fumitools_dedup_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.dedup.bam",
