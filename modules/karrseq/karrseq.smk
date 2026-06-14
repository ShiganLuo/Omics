from snakemake.logging import logger
import time
import os

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")


rule karrseq_chimeric_to_pairs:
    input:
        aligned = indir + "/{sample_id}/{sample_id}_Aligned.sortedByCoord.out.bam",
        chimeric = indir + "/{sample_id}/{sample_id}_Chimeric.out.sam"
    output:
        pairs = outdir + "/{sample_id}/{sample_id}.txt.gz"
    log:
        logdir + "/{sample_id}/karrseq_chimeric_to_pairs.log"
    conda:
        "karrseq.yaml"
    params:
        mapq = config.get("Params", {}).get("karrseq", {}).get("mapq") or 1,
        span = config.get("Params", {}).get("karrseq", {}).get("span") or 0,
        script = os.path.join(ROOT_DIR, "modules/karrseq/bin/get_STAR_reads.py"),
        tmp1 = outdir + "/{sample_id}/{sample_id}.tmp1",
        tmp2 = outdir + "/{sample_id}/{sample_id}.tmp2"
    threads: 4
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start KARRseq chimeric to pairs for sample {wildcards.sample_id} at {current_time}")
        script_file = os.path.join(outdir, f"{wildcards.sample_id}/karrseq_chimeric_to_pairs_{current_time}.sh")
        cmd_view = [
            "samtools", "view", input.aligned,
            "|", "python", params.script, "Aligned", str(params.mapq), str(params.span),
            ">", params.tmp1
        ]
        cmd_chimeric = [
            "cat", input.chimeric,
            "|", "python", params.script, "Chimeric", str(params.mapq), str(params.span),
            ">", params.tmp2
        ]
        cmd_merge = [
            "cat", params.tmp1, params.tmp2,
            "|", "sort", "-k2,2", "-k4,4", "-k3,3n", "-k5,5n", "-k10,10n", "-k11,11n",
            f"--parallel={threads}",
            "|", "gzip", "-c", ">", output.pairs
        ]
        cmd_cleanup = ["rm", params.tmp1, params.tmp2]
        with open(script_file, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd_view) + "\n")
            f.write(" ".join(cmd_chimeric) + "\n")
            f.write(" ".join(cmd_merge) + "\n")
            f.write(" ".join(cmd_cleanup) + "\n")
        shell("bash {script_file} > {log} 2>&1")


rule karrseq_remove_duplicates:
    input:
        pairs = outdir + "/{sample_id}/{sample_id}.txt.gz"
    output:
        dedup = outdir + "/{sample_id}/{sample_id}.dedup.txt.gz",
        bed = outdir + "/{sample_id}/{sample_id}.dedup.bed",
        pairs_gz = outdir + "/{sample_id}/{sample_id}.dedup.pairs.gz"
    log:
        logdir + "/{sample_id}/karrseq_remove_duplicates.log"
    conda:
        "karrseq.yaml"
    params:
        mapq = config.get("Params", {}).get("karrseq", {}).get("mapq") or 1,
        dedup_script = os.path.join(ROOT_DIR, "modules/karrseq/bin/remove_duplicates.py"),
        bed_script = os.path.join(ROOT_DIR, "modules/karrseq/bin/pairs_to_bed.py"),
        todedup = "dedup"
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start KARRseq dedup for sample {wildcards.sample_id} at {current_time}")
        script_file = os.path.join(outdir, f"{wildcards.sample_id}/karrseq_dedup_{current_time}.sh")
        cmd_dedup = [
            "zcat", input.pairs,
            "|", "python", params.dedup_script, params.todedup,
            "|", "gzip", "-c", ">", output.dedup
        ]
        cmd_bed = [
            "zcat", output.dedup,
            "|", "python", params.bed_script, str(params.mapq),
            ">", output.bed
        ]
        cmd_pairs = [
            "zcat", output.dedup,
            "|", "cut", "-f1-7",
            "|", "bgzip", "-c", ">", output.pairs_gz
        ]
        cmd_index = ["pairix", output.pairs_gz]
        with open(script_file, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd_dedup) + "\n")
            f.write(" ".join(cmd_bed) + "\n")
            f.write(" ".join(cmd_pairs) + "\n")
            f.write("sleep 120\n")
            f.write(" ".join(cmd_index) + "\n")
        shell("bash {script_file} > {log} 2>&1")


rule karrseq_ligation:
    input:
        dedup = outdir + "/{sample_id}/{sample_id}.dedup.txt.gz"
    output:
        pairs_gz = outdir + "/{sample_id}/{sample_id}.dedup.ligation.pairs.gz"
    log:
        logdir + "/{sample_id}/karrseq_ligation.log"
    conda:
        "karrseq.yaml"
    threads: 4
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start KARRseq ligation for sample {wildcards.sample_id} at {current_time}")
        script_file = os.path.join(outdir, f"{wildcards.sample_id}/karrseq_ligation_{current_time}.sh")
        cmd_awk1 = [
            "zcat", input.dedup,
            "|", "awk", "'{print $1\"\\t\"$2\"\\t\"$3+$10\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}'"
        ]
        cmd_awk2 = [
            "|", "awk", "'{if($2==$4){if($3>$5){print $1\"\\t\"$4\"\\t\"$5\"\\t\"$2\"\\t\"$3\"\\t\"$7\"\\t\"$6}else{print $0}}else{print $0}}'"
        ]
        cmd_sort = [
            "|", "sort", "-k2,2", "-k4,4", "-k3,3n", "-k5,5n", f"--parallel={threads}"
        ]
        cmd_bgzip = ["|", "bgzip", "-c", ">", output.pairs_gz]
        cmd_index = ["pairix", output.pairs_gz]
        with open(script_file, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd_awk1 + cmd_awk2 + cmd_sort + cmd_bgzip) + "\n")
            f.write("sleep 60\n")
            f.write(" ".join(cmd_index) + "\n")
        shell("bash {script_file} > {log} 2>&1")


rule karrseq_result:
    input:
        pairs_gz = outdir + "/{sample_id}/{sample_id}.dedup.ligation.pairs.gz"
