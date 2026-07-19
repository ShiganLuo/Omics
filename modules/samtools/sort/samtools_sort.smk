include: "../../common/common.smk"

from snakemake.logging import logger
import time
import os
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
fasta = config.get("genome", {}).get("fasta")

rule bam_sort:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        bai = outdir + "/{sample_id}/{sample_id}.bam.bai"
    log:
        logdir + "/{sample_id}/samtools_sort.log"
    threads: 8
    conda:
        "../samtools.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start samtools sort for sample {wildcards.sample_id} at {current_time}")
        cmd1 = [
            params.samtools,
            "sort",
            "-@",
            str(threads),
            "-o",
            output.bam,
            input.bam
        ]
        cmd2 = [
            params.samtools,
            "index",
            "-@",
            str(threads),
            output.bam,
            output.bai
        ]
        script = os.path.join(outdir,f"samtools_sort_{current_time}.sh")
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
        shell("bash {script} > {log} 2>&1")

rule samtools_sort_index_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        bai = outdir + "/{sample_id}/{sample_id}.bam.bai"
