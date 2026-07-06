from snakemake.logging import logger
import time
import os

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])


rule seqprep_merge:
    input:
        read1 = indir + "/{sample_id}/{sample_id}_1.fastq.gz",
        read2 = indir + "/{sample_id}/{sample_id}_2.fastq.gz"
    output:
        read1 = outdir + "/{sample_id}/{sample_id}_unmerge_1.fastq.gz",
        read2 = outdir + "/{sample_id}/{sample_id}_unmerge_2.fastq.gz",
        read3 = outdir + "/{sample_id}/{sample_id}_reject_1.fastq.gz",
        read4 = outdir + "/{sample_id}/{sample_id}_reject_2.fastq.gz",
        merge = outdir + "/{sample_id}/{sample_id}_merge.fastq.gz",
        report = outdir + "/{sample_id}/{sample_id}_merge.report"
    log:
        logdir + "/{sample_id}/seqprep_merge.log"
    params:
        seqprep = config.get("Procedure", {}).get("SeqPrep") or "SeqPrep",
        seqkit = config.get("Procedure", {}).get("seqkit") or "seqkit"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start SeqPrep merge for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"{wildcards.sample_id}/seqprep_merge_{current_time}.sh")
        cmd_merge = [
            params.seqprep,
            "-f", input.read1,
            "-r", input.read2,
            "-1", output.read1,
            "-2", output.read2,
            "-3", output.read3,
            "-4", output.read4,
            "-s", output.merge,
        ]
        cmd_stats = [
            params.seqkit, "stats", output.merge, ">", output.report
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd_merge) + " 2> {log}\n")
            f.write(" ".join(cmd_stats) + " 2>> {log}\n")
        shell("bash {script} > {log} 2>&1")


rule seqprep_result:
    input:
        merge = outdir + "/{sample_id}/{sample_id}_merge.fastq.gz",
        report = outdir + "/{sample_id}/{sample_id}_merge.report"
