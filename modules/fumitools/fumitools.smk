from snakemake.logging import logger
include: "../common/common.smk"
import time
import os
import sys

# Import UMI detection helper from bin/
_fumitools_bin_dir = os.path.join(ROOT_DIR, "modules", "fumitools", "bin")
if _fumitools_bin_dir not in sys.path:
    sys.path.insert(0, _fumitools_bin_dir)
from detect_umi import has_umi_in_header
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
paired_samples = config.get("paired_samples", [])


rule fumitools_copy_umi_paired:
    """Copy UMI from read sequence into FASTQ header."""
    input:
        r1 = indir + "/{sample_id}_1.fq.gz",
        r2 = indir + "/{sample_id}_2.fq.gz",
    output:
        r1 = outdir + "/{sample_id}/{sample_id}_1.umi.fq.gz",
        r2 = outdir + "/{sample_id}/{sample_id}_2.umi.fq.gz",
    log:
        logdir + "/{sample_id}/fumitools_copy_umi.log"
    params:
        fumi_tools = config.get("Procedure", {}).get("fumitools") or "fumi_tools",
        umi_length = config.get("Params", {}).get("fumitools", {}).get("umi_length") or None,
        tag_umi = config.get("Params", {}).get("fumitools", {}).get("tag_umi") or False,
    threads: 4
    conda:
        "fumitools.yaml"
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger = setup_logger(logger_name="fumitools_copy_umi_paired", log_file=log_path)
            # Check if UMI is already present in the read header (check R1)
            if has_umi_in_header(input.r1):
                logger.info(f"UMI already detected in header for {wildcards.sample_id}, symlinking inputs to outputs.")
                os.makedirs(os.path.dirname(output.r1), exist_ok=True)
                for src, dst in [(input.r1, output.r1), (input.r2, output.r2)]:
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.remove(dst)
                    os.symlink(os.path.abspath(src), dst)
                logger.info(f"Symlinked {input.r1} -> {output.r1}")
                logger.info(f"Symlinked {input.r2} -> {output.r2}")
            else:
                script_path = os.path.join(outdir, f"{wildcards.sample_id}/fumitools_copy_umi_{current_time}.sh")
                cmd = [
                    params.fumi_tools, "copy_umi",
                    "-i", input.r1,
                    "-I", input.r2,
                    "-o", output.r1,
                    "-O", output.r2,
                    "--threads", str(threads),
                ]
                if params.umi_length:
                    cmd += ["--umi-length", str(params.umi_length)]
                else:
                    logger.error(f"UMI length must be specified for fumitools_copy_umi_paired for sample {wildcards.sample_id}.")
                    raise ValueError(f"UMI length must be specified for fumitools_copy_umi_paired for sample {wildcards.sample_id}.")
                if params.tag_umi:
                    cmd += ["--tag-umi"]

                with open(script_path, "w") as f:
                    f.write("#!/bin/bash\n")
                    f.write(" ".join(cmd) + "\n")
                shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"fumitools_copy_umi failed for {wildcards.sample_id}: {e}\n")
            raise

rule fumitools_copy_umi_single:
    input:
        r1 = indir + "/{sample_id}.single.fq.gz",
    output:
        r1 = outdir + "/{sample_id}/{sample_id}.umi.single.fq.gz",
    log:
        logdir + "/{sample_id}/fumitools_copy_umi.log"
    params:
        fumi_tools = config.get("Procedure", {}).get("fumitools") or "fumi_tools",
        umi_length = config.get("Params", {}).get("fumitools", {}).get("umi_length") or None,
        tag_umi = config.get("Params", {}).get("fumitools", {}).get("tag_umi") or False,
    threads: 4
    conda:
        "fumitools.yaml"
    run:
        try:
            log_path = str(log)
            open(log_path, "w").close()
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger = setup_logger(logger_name="fumitools_copy_umi_single", log_file=log_path)
            # Check if UMI is already present in the read header
            if has_umi_in_header(input.r1):
                logger.info(f"UMI already detected in header for {wildcards.sample_id}, symlinking input to output.")
                os.makedirs(os.path.dirname(output.r1), exist_ok=True)
                # Remove existing output if it exists (symlink or file)
                if os.path.exists(output.r1) or os.path.islink(output.r1):
                    os.remove(output.r1)
                os.symlink(os.path.abspath(input.r1), output.r1)
                logger.info(f"Symlinked {input.r1} -> {output.r1}")
            else:
                script_path = os.path.join(outdir, f"{wildcards.sample_id}/fumitools_copy_umi_{current_time}.sh")
                cmd = [
                    params.fumi_tools, "copy_umi",
                    "-i", input.r1,
                    "-o", output.r1,
                    "--threads", str(threads),
                ]
                if params.umi_length:
                    cmd += ["--umi-length", str(params.umi_length)]
                else:
                    logger.error(f"UMI length must be specified for fumitools_copy_umi_single for sample {wildcards.sample_id}.")
                    raise ValueError(f"UMI length must be specified for fumitools_copy_umi_single for sample {wildcards.sample_id}.")
                if params.tag_umi:
                    cmd += ["--tag-umi"]
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
    params:
        fumi_tools = config.get("Procedure", {}).get("fumitools") or "fumi_tools",
        memory = config.get("Params", {}).get("fumitools", {}).get("memory", "3G"),
        start_only = config.get("Params", {}).get("fumitools", {}).get("start_only", False),
    threads: config.get("Params", {}).get("fumitools", {}).get("threads", 4)
    conda:
        "fumitools.yaml"
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
