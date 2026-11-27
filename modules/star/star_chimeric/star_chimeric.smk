"""KARRseq chimeric STAR alignment (Kethoxal-Assisted RNA-RNA interaction).

Sub-module of modules/star: single-end STAR run tuned for chimeric (ligation)
junction recovery. Unlike star_align (modules/star/star.smk) the output contract
keeps STAR's raw names — the downstream karrseq module consumes
``<sid>_Chimeric.out.sam`` directly — and forces ``--chimOutType SeparateSAMold``
so that SAM is actually produced (STAR's default Junctions writes
``Chimeric.out.junctions`` instead).
"""

include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
star_index_dir = config.get("genome", {}).get("star_index_dir")
sample_fastqs = config.get("sample_fastqs") or {}

# KARRseq_star_align 参数段:Params.star.KARRseq_star_align.<key>(config 全量
# 路径直读;无 module 级共享回退,KARRseq_star_align 与 star_align 各自独立键)。


def star_fastq(wildcards):
    """Per-sample preprocessed FASTQ (seqprep-merged PE or raw SE, 形态三:
    subworkflow 传 sample_id → 完整路径字典,布局差异不进模块)。"""
    return sample_fastqs[wildcards.sample_id]


def star_genome_index(wildcards):
    """Prebuilt STAR index dir; _require runs at DAG build, not module parse."""
    return _require(star_index_dir, "genome.references.<genome>.star_index_dir", "star_chimeric")


rule KARRseq_star_align:
    input:
        fastq = star_fastq,
        genome_index = star_genome_index,
    output:
        # unregistered stage files consumed by karrseq + sexcheck QC; kept by
        # default, temp when Params.star.keep_intermediates is false
        aligned = maybe_temp(outdir + "/{sample_id}/{sample_id}_Aligned.sortedByCoord.out.bam", "star"),
        chimeric = maybe_temp(outdir + "/{sample_id}/{sample_id}_Chimeric.out.sam", "star"),
        # samtools index -o 命名遵循 sexcheck 的 splitext(bam)+.bai 约定
        bai = maybe_temp(outdir + "/{sample_id}/{sample_id}_Aligned.sortedByCoord.out.bai", "star"),
    log:
        logdir + "/{sample_id}/STAR.log"
    threads: rule_res("KARRseq_star_align").get("threads", 15)
    resources:
        mem_mb = rule_res("KARRseq_star_align").get("mem_mb", 8000),
    conda:
        "../star.yaml"
    container:
        sif("../star.yaml", "KARRseq_star_align")
    params:
        STAR = config.get("Procedure", {}).get("STAR") or "STAR",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
        prefix = outdir + "/{sample_id}/{sample_id}",
        outFilterMultimapNmax = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("outFilterMultimapNmax", 100),
        outSAMattributes = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("outSAMattributes", "All"),
        alignIntronMin = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("alignIntronMin", 1),
        scoreGapNoncan = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("scoreGapNoncan", -4),
        scoreGapATAC = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("scoreGapATAC", -4),
        chimSegmentMin = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("chimSegmentMin", 15),
        chimJunctionOverhangMin = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("chimJunctionOverhangMin", 15),
        limitOutSJcollapsed = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("limitOutSJcollapsed", 10000000),
        limitIObufferSize = config.get("Params", {}).get("star", {}).get("KARRseq_star_align", {}).get("limitIObufferSize", 1500000000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("KARRseq_star_align", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start STAR alignment for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.aligned))
            os.makedirs(sample_outdir, exist_ok=True)  # sample_id 可能含路径分隔符
            script = receipt_script(f"KARRseq_star_align_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "KARRseq_star_align.sh")
            cmd = [
                params.STAR,
                "--runMode", "alignReads",
                "--runThreadN", str(threads),
                "--genomeDir", str(input.genome_index),
                "--readFilesIn", str(input.fastq),
                "--readFilesCommand", "zcat",
                "--outFileNamePrefix", f"{params.prefix}_",
                "--outFilterMultimapNmax", str(params.outFilterMultimapNmax),
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outSAMattributes", params.outSAMattributes,
                "--alignIntronMin", str(params.alignIntronMin),
                "--scoreGapNoncan", str(params.scoreGapNoncan),
                "--scoreGapATAC", str(params.scoreGapATAC),
                "--chimSegmentMin", str(params.chimSegmentMin),
                "--chimJunctionOverhangMin", str(params.chimJunctionOverhangMin),
                # 声明产物 <sid>_Chimeric.out.sam 的前提:Junctions 只出 .junctions
                "--chimOutType", "SeparateSAMold",
                "--outReadsUnmapped", "Fastq",
                "--limitOutSJcollapsed", str(params.limitOutSJcollapsed),
                "--limitIObufferSize", str(params.limitIObufferSize),
            ]
            cmd_index = [params.samtools, "index", "-@", str(threads), "-o", str(output.bai), str(output.aligned)]
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd_index) + "\n")
                f.write(f'echo "STAR alignment for {wildcards.sample_id} was completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during STAR alignment for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(
                f"Error during STAR alignment for sample {wildcards.sample_id}: {e}, "
                f"you can check the log file {log_path} for more details.")


rule star_chimeric_result:
    input:
        aligned = expand(outdir + "/{sample_id}/{sample_id}_Aligned.sortedByCoord.out.bam",
                         sample_id=sample_fastqs.keys()),
        chimeric = expand(outdir + "/{sample_id}/{sample_id}_Chimeric.out.sam",
                          sample_id=sample_fastqs.keys()),
