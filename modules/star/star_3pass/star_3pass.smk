"""Canonical three-pass STAR alignment for small non-coding RNA reads."""

include: "../../common/common.smk"

import os
import shlex

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
STAR = config.get("Procedure", {}).get("STAR") or "STAR"
SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
smallrna_bed = config.get("genome", {}).get("smallrna_bed")
three_pass_params = config.get("Params", {}).get("star_3pass", {})
HELPER = os.path.join(config["ROOT_DIR"], "modules/star/star_3pass/bin/three_pass_align.py")


def fq_inputs(wildcards):
    """Resolve paired- or single-end preprocessed FASTQ inputs."""
    indir = config.get("indir", "input")
    if wildcards.sample_id in config.get("paired_samples", []):
        return [
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz",
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz",
        ]
    return [f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"]


rule star_3p_align:
    input:
        fastq=fq_inputs,
        genome_index=config.get("genome", {}).get("genome_index"),
        smallrna_index=config.get("genome", {}).get("smallrna_index"),
        smallrna_bed=smallrna_bed,
    output:
        pass1_bam=outdir + "/{sample_id}/{sample_id}_pass1.bam",
        pass1_bai=outdir + "/{sample_id}/{sample_id}_pass1.bam.bai",
        pass1_fq1=outdir + "/{sample_id}/{sample_id}_pass1_extract_1.fq.gz",
        pass1_fq2=outdir + "/{sample_id}/{sample_id}_pass1_extract_2.fq.gz",
        pass2_bam=outdir + "/{sample_id}/{sample_id}_pass2.bam",
        pass2_bai=outdir + "/{sample_id}/{sample_id}_pass2.bam.bai",
        pass2_unmapped1=outdir + "/{sample_id}/{sample_id}_pass2_unmapped_fq_1.fq.gz",
        pass2_unmapped2=outdir + "/{sample_id}/{sample_id}_pass2_unmapped_fq_2.fq.gz",
        pass3a_bam=outdir + "/{sample_id}/{sample_id}_pass3a.bam",
        pass3a_bai=outdir + "/{sample_id}/{sample_id}_pass3a.bam.bai",
        pass3b_bam=outdir + "/{sample_id}/{sample_id}_pass3b.bam",
        pass3b_bai=outdir + "/{sample_id}/{sample_id}_pass3b.bam.bai",
        bam=outdir + "/{sample_id}/{sample_id}.bam",
        bai=outdir + "/{sample_id}/{sample_id}.bam.bai",
    log:
        logdir + "/{sample_id}/star3p.log"
    threads: 12
    conda:
        "../star.yaml"
    container:
        sif("../star.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("star_3p_align", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3p_align for sample {wildcards.sample_id} at {current_time}")

            sample_outdir = os.path.dirname(str(output.bam))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            script = os.path.join(sample_outdir, f"three_pass_{current_time}.sh")

            p1 = three_pass_params.get("pass1", {})
            p2 = three_pass_params.get("pass2", {})
            p3 = three_pass_params.get("pass3", {})
            cmd = [
                "python", HELPER,
                "--fastq", *[str(p) for p in input.fastq],
                "--genome-index", str(input.genome_index),
                "--smallrna-index", str(input.smallrna_index),
                "--smallrna-bed", str(input.smallrna_bed),
                "--pass1-bam", str(output.pass1_bam),
                "--pass1-fq1", str(output.pass1_fq1),
                "--pass1-fq2", str(output.pass1_fq2),
                "--pass2-bam", str(output.pass2_bam),
                "--pass2-unmapped1", str(output.pass2_unmapped1),
                "--pass2-unmapped2", str(output.pass2_unmapped2),
                "--pass3a-bam", str(output.pass3a_bam),
                "--pass3b-bam", str(output.pass3b_bam),
                "--bam", str(output.bam),
                "--log", log_path,
                "--threads", str(threads),
                "--star", STAR,
                "--samtools", SAMTOOLS,
                "--bedtools", BEDTOOLS,
            ]
            # Pass 1 parameters
            if "outFilterMultimapNmax" in p1:
                cmd += ["--pass1-outFilterMultimapNmax", str(p1["outFilterMultimapNmax"])]
            if "outFilterMultimapScoreRange" in p1:
                cmd += ["--pass1-outFilterMultimapScoreRange", str(p1["outFilterMultimapScoreRange"])]
            if "outFilterMismatchNoverLmax" in p1:
                cmd += ["--pass1-outFilterMismatchNoverLmax", str(p1["outFilterMismatchNoverLmax"])]
            if "alignIntronMin" in p1:
                cmd += ["--pass1-alignIntronMin", str(p1["alignIntronMin"])]
            if "alignMatesGapMax" in p1:
                cmd += ["--pass1-alignMatesGapMax", str(p1["alignMatesGapMax"])]
            # Pass 2 parameters
            if "outFilterMultimapNmax" in p2:
                cmd += ["--pass2-outFilterMultimapNmax", str(p2["outFilterMultimapNmax"])]
            if "outFilterMultimapScoreRange" in p2:
                cmd += ["--pass2-outFilterMultimapScoreRange", str(p2["outFilterMultimapScoreRange"])]
            if "outFilterMismatchNoverLmax" in p2:
                cmd += ["--pass2-outFilterMismatchNoverLmax", str(p2["outFilterMismatchNoverLmax"])]
            if "outFilterMismatchNoverReadLmax" in p2:
                cmd += ["--pass2-outFilterMismatchNoverReadLmax", str(p2["outFilterMismatchNoverReadLmax"])]
            if "alignIntronMin" in p2:
                cmd += ["--pass2-alignIntronMin", str(p2["alignIntronMin"])]
            if "alignMatesGapMax" in p2:
                cmd += ["--pass2-alignMatesGapMax", str(p2["alignMatesGapMax"])]
            if "alignEndsType" in p2:
                cmd += ["--pass2-alignEndsType", str(p2["alignEndsType"])]
            if "clip5pNbases" in p2:
                cmd += ["--pass2-clip5pNbases", str(p2["clip5pNbases"])]
            if "clip3pNbases" in p2:
                cmd += ["--pass2-clip3pNbases", str(p2["clip3pNbases"])]
            # Pass 3 parameters
            if "outFilterMultimapNmax" in p3:
                cmd += ["--pass3-outFilterMultimapNmax", str(p3["outFilterMultimapNmax"])]
            if "outFilterMultimapScoreRange" in p3:
                cmd += ["--pass3-outFilterMultimapScoreRange", str(p3["outFilterMultimapScoreRange"])]
            if "outFilterMismatchNoverLmax" in p3:
                cmd += ["--pass3-outFilterMismatchNoverLmax", str(p3["outFilterMismatchNoverLmax"])]
            if "alignIntronMin" in p3:
                cmd += ["--pass3-alignIntronMin", str(p3["alignIntronMin"])]
            if "alignMatesGapMax" in p3:
                cmd += ["--pass3-alignMatesGapMax", str(p3["alignMatesGapMax"])]
            if "alignEndsType" in p3:
                cmd += ["--pass3-alignEndsType", str(p3["alignEndsType"])]

            with open(script, "w") as f:
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"star_3p_align for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"star_3p_align failed for sample {wildcards.sample_id}: {e}\n")
            raise e
