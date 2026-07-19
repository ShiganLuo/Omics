include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
smallrna_bed = config.get("genome", {}).get("smallrna_bed")

# ── Extract reads overlapping small RNA genes (used after pass 1 & pass 3a) ─
rule star_3p_extract_smallrna:
    input:
        bam = outdir + "/pass1/{sample_id}/{sample_id}.Aligned.sortedByCoord.out.bam",
        bai = outdir + "/pass1/{sample_id}/{sample_id}.Aligned.sortedByCoord.out.bam.bai",
        bed = smallrna_bed,
    output:
        bam = outdir + "/pass1_extract/{sample_id}/{sample_id}.smallrna.bam",
        fq  = outdir + "/pass1_extract/{sample_id}/{sample_id}.smallrna.fq.gz",
    log:
        logdir + "/star3p/{sample_id}/extract_smallrna.log"
    threads: 4
    conda:
        "star_3pass.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("star_3p_extract_smallrna", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3p_extract_smallrna for sample {wildcards.sample_id} at {current_time}")
            outdir_bam = os.path.dirname(str(output.bam))
            os.makedirs(outdir_bam, exist_ok=True)
            script = os.path.join(outdir_bam, f"extract_smallrna_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"mkdir -p $(dirname {output.bam})\n")
                f.write(f"{BEDTOOLS} intersect -abam {input.bam} -b {input.bed} -u \\\n")
                f.write(f"    | {SAMTOOLS} sort -n -@ {threads} -T $(dirname {output.bam})/tmp \\\n")
                f.write(f"    -o {output.bam}\n")
                f.write(f"{SAMTOOLS} fastq -@ {threads} {output.bam} \\\n")
                f.write(f"    | gzip > {output.fq}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_3p_extract_smallrna for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during star_3p_extract_smallrna for sample {wildcards.sample_id}: {e}")
            raise e

# ── Extract canonical reads from pass 3a (still overlapping small RNA genes) ─
rule star_3p_pass3a_extract:
    input:
        bam = outdir + "/pass3a/{sample_id}/{sample_id}.Aligned.sortedByCoord.out.bam",
        bed = smallrna_bed,
    output:
        bam = outdir + "/pass3a_extract/{sample_id}/{sample_id}.canonical.bam",
        bai = outdir + "/pass3a_extract/{sample_id}/{sample_id}.canonical.bam.bai",
    log:
        logdir + "/star3p/{sample_id}/pass3a_extract.log"
    threads: 4
    conda:
        "star_3pass.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("star_3p_pass3a_extract", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3p_pass3a_extract for sample {wildcards.sample_id} at {current_time}")
            outdir_bam = os.path.dirname(str(output.bam))
            os.makedirs(outdir_bam, exist_ok=True)
            script = os.path.join(outdir_bam, f"pass3a_extract_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"mkdir -p $(dirname {output.bam})\n")
                f.write(f"{BEDTOOLS} intersect -abam {input.bam} -b {input.bed} -u \\\n")
                f.write(f"    > {output.bam}\n")
                f.write(f"{SAMTOOLS} index -@ {threads} {output.bam}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_3p_pass3a_extract for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during star_3p_pass3a_extract for sample {wildcards.sample_id}: {e}")
            raise e

# ── Merge canonical + non-canonical reads ────────────────────────────────
rule star_3p_merge:
    input:
        canonical = outdir + "/pass3a_extract/{sample_id}/{sample_id}.canonical.bam",
        noncanonical = outdir + "/pass3b/{sample_id}/{sample_id}.Aligned.sortedByCoord.out.bam",
    output:
        bam = outdir + "/{sample_id}/{sample_id}.Aligned.sortedByCoord.out.bam",
        bai = outdir + "/{sample_id}/{sample_id}.Aligned.sortedByCoord.out.bam.bai",
    log:
        logdir + "/star3p/{sample_id}/merge.log"
    threads: 4
    conda:
        "star_3pass.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("star_3p_merge", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3p_merge for sample {wildcards.sample_id} at {current_time}")
            outdir_bam = os.path.dirname(str(output.bam))
            os.makedirs(outdir_bam, exist_ok=True)
            script = os.path.join(outdir_bam, f"merge_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"mkdir -p $(dirname {output.bam})\n")
                f.write(f"{SAMTOOLS} merge -@ {threads} -f {output.bam} \\\n")
                f.write(f"    {input.canonical} {input.noncanonical}\n")
                f.write(f"{SAMTOOLS} index -@ {threads} {output.bam}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_3p_merge for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during star_3p_merge for sample {wildcards.sample_id}: {e}")
            raise e
