from snakemake.logging import logger
import time
import os
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")

rule pbmm2_align:
    input:
        bam = indir + "/{sample_id}.bam",
        fasta = fasta
    output:
        bam = outdir + "/{sample_id}/{sample_id}.sorted.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted.bai"
    log:
        logdir + "/{sample_id}/pbmm2_align.log"
    threads: 16
    conda:
        "pbmm2.yaml"
    params:
        pbmm2 = config.get("Procedure", {}).get("pbmm2") or "pbmm2",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start pbmm2 alignment for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"pbmm2_align_{current_time}.sh")
        raw_align_bam = outdir + f"/{wildcards.sample_id}/{wildcards.sample_id}.raw.bam"
        cmd1 = [
            params.pbmm2, "align",
            "--sort",
            "--num-threads", str(threads),
            input.fasta,
            input.bam,
            raw_align_bam
        ]
        cmd2 = [
            params.samtools, "sort",
            "-@", str(threads),
            raw_align_bam,
            "-o", output.bam
        ]
        cmd3 = [
            params.samtools, "index",
            "-@", str(threads),
            output.bam,
            "-o", output.bai
        ]
        cmd4 = [
            "rm", raw_align_bam
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(" ".join(cmd3) + "\n")
            f.write(" ".join(cmd4) + "\n")
        shell(f"bash {script} > {log} 2>&1")

rule pbmm2_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted.bai"
