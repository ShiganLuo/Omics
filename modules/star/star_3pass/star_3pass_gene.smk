"""Strict gene-specific three-pass local alignment and Tailer after canonical star_3pass."""

include: "../../common/common.smk"


outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
final_bam_dir = config.get("final_bam_dir", "")
STAR = config.get("Procedure", {}).get("STAR") or "STAR"
SAMTOOLS = config.get("Procedure", {}).get("samtools") or "samtools"
BEDTOOLS = config.get("Procedure", {}).get("bedtools") or "bedtools"
TAILER = config.get("Procedure", {}).get("tailer") or "Tailer"
smallrna_bed = config.get("genome", {}).get("smallrna_bed")
genome_fasta = config.get("genome", {}).get("genome_fasta")

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
        logdir + "/{sample_id}/star3p_gene_specific.log"
    threads: rule_res("star_3pg_gene_specific").get("threads", 2)
    benchmark: logdir + "/{sample_id}/star3p_gene_specific.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "star_3pg_gene_specific", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../star.yaml"
    container:
        sif("../star.yaml", "star_3pg_gene_specific")
    params:
        # 工具/参数挂 Params 的 rule 子键(config 全量路径直读,run 只碰 params)
        HELPER = os.path.join(config["ROOT_DIR"], "modules/star/star_3pass/bin/gene_specific_align.py"),
        PREPARE_HELPER = os.path.join(config["ROOT_DIR"], "modules/star/star_3pass/bin/prepare_gene_inputs.py"),
        STAR = STAR,
        SAMTOOLS = SAMTOOLS,
        BEDTOOLS = BEDTOOLS,
        TAILER = TAILER,
        ambiguous = config.get("Params", {}).get("star_3pass_gene", {}).get("star_3pg_gene_specific", {}).get("ambiguous") or "exclude",
        flank = int(config.get("Params", {}).get("star_3pass_gene", {}).get("star_3pg_gene_specific", {}).get("flank") or 50),
        genome_sa_index_nbases = int(config.get("Params", {}).get("star_3pass_gene", {}).get("index", {}).get("genomeSAindexNbases") or 3),
        pass1 = config.get("Params", {}).get("star_3pass_gene", {}).get("passes", {}).get("pass1", {}),
        pass2 = config.get("Params", {}).get("star_3pass_gene", {}).get("passes", {}).get("pass2", {}),
        pass3 = config.get("Params", {}).get("star_3pass_gene", {}).get("passes", {}).get("pass3", {}),
        # tailer 参数:Params.tailer.star_3pg_gene_specific.<key>(与 tailer_global 各自独立键)
        tailer_threshold = config.get("Params", {}).get("tailer", {}).get("star_3pg_gene_specific", {}).get("threshold", None),
        tailer_rev_comp = config.get("Params", {}).get("tailer", {}).get("star_3pg_gene_specific", {}).get("rev_comp", False),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("star_3pg_gene_specific", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_3pg_gene_specific for sample {wildcards.sample_id} at {current_time}")

            sample_outdir = os.path.dirname(str(output.bam))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            input_dir = os.path.join(outdir, wildcards.sample_id)
            script = os.path.join(RECEIPT_DIR, f"gene_specific_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"gene_specific.sh")
            p1, p2, p3 = params.pass1, params.pass2, params.pass3

            cmd = [
                "python", params.HELPER,
                "--input-bam", str(input.bam[0]),
                "--input-bed", str(input.bed),
                "--input-genome-fasta", str(input.genome_fasta),
                "--output-bam", str(output.bam),
                "--output-tail", str(output.tail),
                "--input-dir", input_dir,
                "--helper", params.PREPARE_HELPER,
                "--log", log_path,
                "--threads", str(threads),
                "--star", params.STAR,
                "--samtools", params.SAMTOOLS,
                "--bedtools", params.BEDTOOLS,
                "--tailer", params.TAILER,
                "--ambiguous", params.ambiguous,
                "--flank", str(params.flank),
            ]
            # Index parameters (genomeGenerate)
            cmd += ["--index-genome-sa-index-nbases", str(params.genome_sa_index_nbases)]
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
            # Tailer parameters (rule 级键,与 tailer_global 各自独立)
            if params.tailer_threshold is not None:
                cmd += ["--tailer-threshold", str(int(params.tailer_threshold))]
            if params.tailer_rev_comp:
                cmd.append("--tailer-rev-comp")

            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'star_3pg_gene_specific for sample {wildcards.sample_id} completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"star_3pg_gene_specific failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"star_3pg_gene_specific failed for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")
