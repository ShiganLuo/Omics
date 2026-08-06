"""Strict gene-specific three-pass local alignment and Tailer after canonical star_3pass."""

include: "../../common/common.smk"

import os
import shlex

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
final_bam_dir = config.get("final_bam_dir", "")
STAR = config.get("Procedure", {}).get("STAR") or "STAR"
SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
TAILER = config.get("Procedure", {}).get("tailer") or "Tailer"
smallrna_bed = config.get("genome", {}).get("smallrna_bed")
genome_fasta = config.get("genome", {}).get("genome_fasta")
params = config.get("Params", {}).get("star_3pass_gene", {})
passes = params.get("passes", {})
p1 = passes.get("pass1", {})
p2 = passes.get("pass2", {})
p3 = passes.get("pass3", {})
ambiguous = params.get("ambiguous", "exclude")
flank = int(params.get("flank", 50))
tailer_params = config.get("Params", {}).get("tailer", {})
HELPER = os.path.join(config["ROOT_DIR"], "modules/star/star_3pass/bin/gene_specific_align.py")
PREPARE_HELPER = os.path.join(config["ROOT_DIR"], "modules/star/star_3pass/bin/prepare_gene_inputs.py")

if not final_bam_dir or not smallrna_bed or not genome_fasta:
    raise ValueError("star_3pass_gene requires final_bam_dir, smallrna_bed, and genome_fasta")


def final_bam(wildcards):
    return [f"{final_bam_dir}/{wildcards.sample_id}/{wildcards.sample_id}.bam",
            f"{final_bam_dir}/{wildcards.sample_id}/{wildcards.sample_id}.bam.bai"]


rule star_3pg_gene_specific:
    input:
        bam=final_bam,
        bed=smallrna_bed,
        genome_fasta=genome_fasta,
    output:
        bam=outdir + "/{sample_id}/{sample_id}.bam",
        bai=outdir + "/{sample_id}/{sample_id}.bam.bai",
        tail=outdir + "/{sample_id}/{sample_id}_tail.csv",
    log:
        logdir + "/star3pg/{sample_id}/gene_specific.log"
    threads: 2
    conda:
        "../star.yaml"
    container:
        sif("../star.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("star_3pg_gene_specific", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3pg_gene_specific for sample {wildcards.sample_id} at {current_time}")

            sample_outdir = os.path.dirname(str(output.bam))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            input_dir = os.path.join(outdir, wildcards.sample_id)
            script = os.path.join(sample_outdir, f"gene_specific_{current_time}.sh")

            cmd = [
                "python", HELPER,
                "--input-bam", str(input.bam[0]),
                "--input-bed", str(input.bed),
                "--input-genome-fasta", str(input.genome_fasta),
                "--output-bam", str(output.bam),
                "--output-tail", str(output.tail),
                "--input-dir", input_dir,
                "--helper", PREPARE_HELPER,
                "--log", log_path,
                "--threads", str(threads),
                "--star", STAR,
                "--samtools", SAMTOOLS,
                "--bedtools", BEDTOOLS,
                "--tailer", TAILER,
                "--ambiguous", ambiguous,
                "--flank", str(flank),
            ]
            # Pass 1 parameters (includes genomeSAindexNbases for genomeGenerate)
            if "genomeSAindexNbases" in p1:
                cmd += ["--pass1-genome-sa-index-nbases", str(p1["genomeSAindexNbases"])]
            if "outFilterMultimapNmax" in p1:
                cmd += ["--pass1-out-filter-multimap-nmax", str(p1["outFilterMultimapNmax"])]
            if "outFilterMultimapScoreRange" in p1:
                cmd += ["--pass1-out-filter-multimap-score-range", str(p1["outFilterMultimapScoreRange"])]
            if "outFilterMismatchNoverLmax" in p1:
                cmd += ["--pass1-out-filter-mismatch-nover-lmax", str(p1["outFilterMismatchNoverLmax"])]
            if "alignIntronMin" in p1:
                cmd += ["--pass1-align-intron-min", str(p1["alignIntronMin"])]
            if "alignMatesGapMax" in p1:
                cmd += ["--pass1-align-mates-gap-max", str(p1["alignMatesGapMax"])]
            if "alignEndsType" in p1:
                cmd += ["--pass1-align-ends-type", str(p1["alignEndsType"])]
            # Pass 2 parameters
            if "outFilterMultimapNmax" in p2:
                cmd += ["--pass2-out-filter-multimap-nmax", str(p2["outFilterMultimapNmax"])]
            if "outFilterMultimapScoreRange" in p2:
                cmd += ["--pass2-out-filter-multimap-score-range", str(p2["outFilterMultimapScoreRange"])]
            if "outFilterMismatchNoverLmax" in p2:
                cmd += ["--pass2-out-filter-mismatch-nover-lmax", str(p2["outFilterMismatchNoverLmax"])]
            if "alignIntronMin" in p2:
                cmd += ["--pass2-align-intron-min", str(p2["alignIntronMin"])]
            if "alignMatesGapMax" in p2:
                cmd += ["--pass2-align-mates-gap-max", str(p2["alignMatesGapMax"])]
            if "alignEndsType" in p2:
                cmd += ["--pass2-align-ends-type", str(p2["alignEndsType"])]
            # Pass 3 parameters
            if "outFilterMultimapNmax" in p3:
                cmd += ["--pass3-out-filter-multimap-nmax", str(p3["outFilterMultimapNmax"])]
            if "outFilterMultimapScoreRange" in p3:
                cmd += ["--pass3-out-filter-multimap-score-range", str(p3["outFilterMultimapScoreRange"])]
            if "outFilterMismatchNoverLmax" in p3:
                cmd += ["--pass3-out-filter-mismatch-nover-lmax", str(p3["outFilterMismatchNoverLmax"])]
            if "alignIntronMin" in p3:
                cmd += ["--pass3-align-intron-min", str(p3["alignIntronMin"])]
            if "alignMatesGapMax" in p3:
                cmd += ["--pass3-align-mates-gap-max", str(p3["alignMatesGapMax"])]
            if "alignEndsType" in p3:
                cmd += ["--pass3-align-ends-type", str(p3["alignEndsType"])]
            # Tailer parameters
            if "threshold" in tailer_params:
                cmd += ["--tailer-threshold", str(int(tailer_params["threshold"]))]
            if tailer_params.get("rev_comp", False):
                cmd.append("--tailer-rev-comp")

            with open(script, "w") as f:
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"star_3pg_gene_specific for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"star_3pg_gene_specific failed for sample {wildcards.sample_id}: {e}\n")
            raise e
