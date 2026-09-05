"""Strict gene-specific three-pass local alignment and Tailer after canonical star_3pass (polygenomes)."""

include: "../../../common/common.smk"

import os
import shlex

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
final_bam_dir = config.get("final_bam_dir", "")
STAR = config.get("Procedure", {}).get("STAR") or "STAR"
SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
TAILER = config.get("Procedure", {}).get("tailer") or "Tailer"
params = config.get("Params", {}).get("star_3pass_gene", {})
passes = params.get("passes", {})
p1 = passes.get("pass1", {})
p2 = passes.get("pass2", {})
p3 = passes.get("pass3", {})
ambiguous = params.get("ambiguous", "exclude")
flank = int(params.get("flank", 50))
index_params = params.get("index", {})
genome_sa_index_nbases = int(index_params.get("genomeSAindexNbases", 3))
tailer_params = config.get("Params", {}).get("tailer", {})
HELPER = os.path.join(config["ROOT_DIR"], "modules/star/star_3pass/bin/gene_specific_align.py")
PREPARE_HELPER = os.path.join(config["ROOT_DIR"], "modules/star/star_3pass/bin/prepare_gene_inputs.py")


def final_bam(wildcards):
    return [f"{final_bam_dir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.bam",
            f"{final_bam_dir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.bam.bai"]


rule star_3pg_gene_specific:
    input:
        bam=final_bam,
        bed=lambda wildcards: config["genome"][wildcards.genome]["smallrna_bed"],
        genome_fasta=lambda wildcards: config["genome"][wildcards.genome]["genome_fasta"],
    output:
        bam=outdir + "/{genome}/{sample_id}/{sample_id}.bam",
        bai=outdir + "/{genome}/{sample_id}/{sample_id}.bam.bai",
        tail=outdir + "/{genome}/{sample_id}/{sample_id}_tail.csv",
    log:
        logdir + "/{sample_id}/{genome}/star3p_gene_specific.log"
    threads: 2
    conda:
        "../../star.yaml"
    container:
        sif("../../star.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("star_3pg_gene_specific", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3pg_gene_specific for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")

            sample_outdir = os.path.dirname(str(output.bam))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            input_dir = os.path.join(outdir, wildcards.genome, wildcards.sample_id)
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
            # Index parameters (genomeGenerate)
            cmd += ["--index-genome-sa-index-nbases", str(genome_sa_index_nbases)]
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
            # Tailer parameters
            if "threshold" in tailer_params:
                cmd += ["--tailer-threshold", str(int(tailer_params["threshold"]))]
            if tailer_params.get("rev_comp", False):
                cmd.append("--tailer-rev-comp")

            with open(script, "w") as f:
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"star_3pg_gene_specific for sample {wildcards.sample_id} genome {wildcards.genome} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"star_3pg_gene_specific failed for sample {wildcards.sample_id} genome {wildcards.genome}: {e}\n")
            raise e
