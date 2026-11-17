include: "../common/common.smk"

import logging
logger = logging.getLogger(__name__)
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_index = config.get("logdir_index", logdir)
log_index_substring = config.get("log_index_substring", "star_index")
indir= config.get("indir", "input")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
fasta = config.get('genome',{}).get('fasta')
is_fasta_intermediate = config.get('genome',{}).get('is_fasta_intermediate', False)
omics_type = config.get('omics_type', 'RNA-seq')
if not is_fasta_intermediate and (not fasta or not os.path.exists(fasta)):
    raise ValueError(f"FASTA file not found: {fasta}. Please provide a valid genome FASTA file in the config.")
gtf = config.get('genome',{}).get('gtf')
fastq_sample_suffix = config.get('fastq_sample_suffix') or None

# 参数精确到 rule:Params.star.<rule_name>.<key>(rule_param 见 common.smk;
# 无 module 级共享回退,star_align 独立键;index 子段仍为 Params.star.index)。
_star_params = config.get("Params", {}).get("star", {})

rule star_index:
    input:
        fasta = fasta,
    output:
        index_file = directory(outdir + "/index")
    log:
        logdir_index + "/star/" + log_index_substring + ".log"
    threads: rule_res("star_index").get("threads", 12)
    benchmark: logdir_index + "/star/" + log_index_substring + ".benchmark.tsv",
    resources:
        mem_mb=rule_res("star_index").get("mem_mb", 64000),
    conda:
        "star.yaml"
    container:
        sif("star.yaml", "star_index")
    params:
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        index_dir = outdir + "/index",
        sjdbOverhang = config.get('Params',{}).get('star', {}).get('index', {}).get('sjdbOverhang') or 100,
        genomeSAindexNbases = config.get('Params',{}).get('star', {}).get('index', {}).get('genomeSAindexNbases'),
        gtf = gtf,
        outTmpDir = outdir + "/tmp_star"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("star_index", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_index at {current_time}")

            sample_outdir = os.path.dirname(str(output.index_file))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"star_index_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"star_index.sh")

            cmd = [
                params.STAR, "--runMode", "genomeGenerate",
                "--runThreadN", str(threads),
                "--genomeDir", params.index_dir,
                "--genomeFastaFiles", input.fasta]
            # --genomeSAindexNbases: explicit config wins; otherwise auto-compute
            # STAR's recommendation min(14, log2(genomeSize)/2 - 1) — required
            # for small genomes (default 14 breaks mapping).
            if params.genomeSAindexNbases:
                sa_index_bases = int(params.genomeSAindexNbases)
                rule_logger.info(f"Using genomeSAindexNbases from config: {sa_index_bases}")
            else:
                # chunked byte scan: count non-header bytes (fast, ~500MB/s+)
                genome_size = 0
                with open(input.fasta, "rb") as _fh:
                    for _chunk in iter(lambda: _fh.read(1 << 24), b""):
                        for _line in _chunk.split(b"\n"):
                            if not _line.startswith(b">"):
                                genome_size += len(_line.strip())
                import math
                sa_index_bases = max(2, min(14, int(math.log2(genome_size) / 2) - 1)) if genome_size > 0 else 14
                rule_logger.info(f"Genome size {genome_size} bp -> auto genomeSAindexNbases {sa_index_bases}")
            cmd.extend(["--genomeSAindexNbases", str(sa_index_bases)])
            if params.outTmpDir:
                shutil.rmtree(params.outTmpDir, ignore_errors=True)
                cmd.extend(["--outTmpDir", params.outTmpDir])
            if params.gtf:
                # STAR can't use --sjdbOverhang without an annotation file
                cmd.extend(["--sjdbGTFfile", str(params.gtf)])
                cmd.extend(["--sjdbOverhang", str(params.sjdbOverhang)])
                rule_logger.info(f"Using sjdbGTFfile: {params.gtf}")
            else:
                rule_logger.info("No GTF provided, skipping --sjdbGTFfile")

            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                # Run STAR from inside the genome dir: Log.out and other
                # run-time files land there and STAR's end-of-run move stays
                # on one filesystem (cross-mount rename fails under
                # singularity, emitting the "Could not move Log.out" warning).
                f.write(f'cd "{sample_outdir}"\n')
                f.write(" ".join(cmd) + "\n")
                f.write("echo 'STAR index generation completed successfully.'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during star_index: {e}")
            raise RuntimeError(f"Error occurred during star_index: {e}, you can check the log file {log_path} for more details.")

def get_star_index(wildcards):
    logger.debug(f"[get_star_index] called with wildcards: {wildcards}")
    star_index_dir = config.get('genome',{}).get('index_dir') or None
    if star_index_dir:
        logger.debug(f"[get_star_index] using provided star_index_dir: {star_index_dir}")
        return star_index_dir
    logger.debug(f"[get_star_index] using default star_index_dir")
    return outdir + f"/index"

def get_alignment_input(wildcards):
    """
    function: Dynamically determines the input file type: paired-end or single-end sequencing.
    Based on the paired_samples and single_samples lists.This function is called in the star_align rule.

    param: 
        wildcards: Snakemake wildcards object containing the sample_id.
        paired_samples = ['sample1', 'sample2', ...]
        single_samples = ['sample3', 'sample4', ...]
    These lists must be defined in the Snakefile or config file.

    return: A list of input file paths for the STAR alignment step. 
    """
    logger.debug(f"[get_alignment_input] called with wildcards: {wildcards}")
    # 构造可能的输入路径
    if fastq_sample_suffix:
        paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}_1.fq.gz"
        paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}_2.fq.gz"
        single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}.single.fq.gz"
    else:
        paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
        paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
        single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

    # 检查文件实际存在情况
    if wildcards.sample_id in paired_samples:
        logger.debug(f"双端测序：{[paired_r1, paired_r2]}")
        logger.debug(f"双端测序：{[paired_r1, paired_r2]}")
        if omics_type == "scRNAseq":
            logger.debug(f"Detected scRNA-seq data for sample {wildcards.sample_id}.R1 is expected to contain cell barcodes and UMIs, R2 contains the transcript sequence.")
            return [paired_r2, paired_r1]
        else:   
            return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        logger.debug(f"单端测序：{[single]}")
        logger.debug(f"单端测序：{[single]}")
        return [single]
    else:
        logger.error(f"样本 {wildcards.sample_id} 未在 paired_samples 或 single_samples 中定义")
        raise ValueError(f"Sample {wildcards.sample_id} not defined in paired_samples or single_samples")

rule star_align:
    input:
        fastq = get_alignment_input,
        genome_index = get_star_index
    output:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        bai = outdir + "/{sample_id}/{sample_id}.bam.bai",
        unmapped_r1 = outdir + "/{sample_id}/{sample_id}.Unmapped.out.mate1",
        unmapped_r2 = outdir + "/{sample_id}/{sample_id}.Unmapped.out.mate2",
    log:
        logdir + "/{sample_id}/star_align.log"
    threads: rule_res("star_align").get("threads", 12)
    benchmark: logdir + "/{sample_id}/star_align.benchmark.tsv",
    resources:
        mem_mb=lambda wildcards, input, attempt: rule_res(
            "star_align", input_size=input_size_named(input, 'fastq'),
            attempt=attempt).get("mem_mb", 32000),
    params:
        outPrefix = outdir + "/{sample_id}/{sample_id}.",
        input_params = lambda wildcards, input: \
            f"{input.fastq[0]} {input.fastq[1]}" if len(input.fastq) == 2 else f"{input.fastq[0]}",
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        SAMTOOLS = config.get('Procedure',{}).get('samtools') or 'samtools',
        twopassMode = rule_param(_star_params, "star_align", "twopassMode", None),
        alignEndsType = rule_param(_star_params, "star_align", "alignEndsType", "Local"),
        outFilterMismatchNoverReadLmax = rule_param(_star_params, "star_align", "outFilterMismatchNoverReadLmax", 1.0),
        outFilterMismatchNmax = rule_param(_star_params, "star_align", "outFilterMismatchNmax", 10),
        outFilterMultimapNmax = rule_param(_star_params, "star_align", "outFilterMultimapNmax", 10),
        winAnchorMultimapNmax = rule_param(_star_params, "star_align", "winAnchorMultimapNmax", 50),
        genomeLoad = rule_param(_star_params, "star_align", "genomeLoad", 'NoSharedMemory'),
        limitBAMsortRAM = rule_param(_star_params, "star_align", "limitBAMsortRAM", 0),
        outReadsUnmapped = rule_param(_star_params, "star_align", "outReadsUnmapped", None),
        outFilterMismatchNoverLmax = rule_param(_star_params, "star_align", "outFilterMismatchNoverLmax", 0.3),
        outFilterMatchNminOverLread = rule_param(_star_params, "star_align", "outFilterMatchNminOverLread", 0.66),
        alignSJoverhangMin = rule_param(_star_params, "star_align", "alignSJoverhangMin", 5),
        alignSJDBoverhangMin = rule_param(_star_params, "star_align", "alignSJDBoverhangMin", 3),
        chimSegmentMin = rule_param(_star_params, "star_align", "chimSegmentMin", 0),
        chimOutType = rule_param(_star_params, "star_align", "chimOutType", "Junctions"),
        chimJunctionOverhangMin = rule_param(_star_params, "star_align", "chimJunctionOverhangMin", 20),
        outSAMstrandField = rule_param(_star_params, "star_align", "outSAMstrandField", None),
        chimScoreMin = rule_param(_star_params, "star_align", "chimScoreMin", 0),
        chimScoreDropMax = rule_param(_star_params, "star_align", "chimScoreDropMax", 20),
        chimScoreJunctionNonGTAG = rule_param(_star_params, "star_align", "chimScoreJunctionNonGTAG", -1),
        chimScoreSeparation = rule_param(_star_params, "star_align", "chimScoreSeparation", 10),
        alignSJstitchMismatchNmax = rule_param(_star_params, "star_align", "alignSJstitchMismatchNmax", "0 -1 0 0"),
        chimSegmentReadGapMax = rule_param(_star_params, "star_align", "chimSegmentReadGapMax", 0),
        outSAMattributes = rule_param(_star_params, "star_align", "outSAMattributes", "NM"),
        outMultimapperOrder = rule_param(_star_params, "star_align", "outMultimapperOrder", "Old_2.4"),
        runRNGseed = rule_param(_star_params, "star_align", "runRNGseed", 777),
        outSAMmultNmax = rule_param(_star_params, "star_align", "outSAMmultNmax", -1),
        soloType = rule_param(_star_params, "star_align", "soloType", None),
        soloCBwhitelist = rule_param(_star_params, "star_align", "soloCBwhitelist", None),
        soloBarcodeReadLength = rule_param(_star_params, "star_align", "soloBarcodeReadLength", 1), # may be 0
        limitSjdbInsertNsj = rule_param(_star_params, "star_align", "limitSjdbInsertNsj", 1000000),
        outTmpDir = outdir + "/{sample_id}/tmp_star"
    conda:
        "star.yaml"
    container:
        sif("star.yaml", "star_align")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger(logger_name="star_align", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = os.path.join(RECEIPT_DIR, f"star_align_{wildcards.sample_id}_{current_time}.sh")
            plain_script = f"{outdir}/{wildcards.sample_id}/star_align.sh"
            rule_logger.info(f"Start STAR alignment for sample {wildcards.sample_id} at {current_time}")
            cmd1 = [
                params.STAR, "--runThreadN", str(threads),
                "--genomeDir", input.genome_index,
                "--readFilesCommand", "zcat",
                "--genomeLoad", params.genomeLoad,
                "--limitBAMsortRAM", str(params.limitBAMsortRAM),
                "--alignEndsType", params.alignEndsType,
                "--winAnchorMultimapNmax", str(params.winAnchorMultimapNmax),
                "--outFilterMismatchNmax", str(params.outFilterMismatchNmax),
                "--outFilterMultimapNmax", str(params.outFilterMultimapNmax),
                "--outFilterMismatchNoverLmax", str(params.outFilterMismatchNoverLmax),
                "--outFilterMatchNminOverLread", str(params.outFilterMatchNminOverLread),
                "--alignSJoverhangMin", str(params.alignSJoverhangMin),
                "--alignSJDBoverhangMin", str(params.alignSJDBoverhangMin),
                "--chimSegmentMin", str(params.chimSegmentMin),
                "--chimOutType", params.chimOutType,
                "--chimJunctionOverhangMin", str(params.chimJunctionOverhangMin),
                "--chimScoreMin", str(params.chimScoreMin),
                "--chimScoreDropMax", str(params.chimScoreDropMax),
                "--chimScoreJunctionNonGTAG", str(params.chimScoreJunctionNonGTAG),
                "--chimScoreSeparation", str(params.chimScoreSeparation),
                "--alignSJstitchMismatchNmax", params.alignSJstitchMismatchNmax,
                "--chimSegmentReadGapMax", str(params.chimSegmentReadGapMax),
                "--outSAMtype", "BAM SortedByCoordinate",
                "--outSAMattributes", str(params.outSAMattributes),
                "--outMultimapperOrder", str(params.outMultimapperOrder),
                "--runRNGseed", str(params.runRNGseed),
                "--outSAMmultNmax", str(params.outSAMmultNmax),
                "--soloType", str(params.soloType),
                "--soloCBwhitelist", str(params.soloCBwhitelist),
                "--soloBarcodeReadLength", str(params.soloBarcodeReadLength),
                "--limitSjdbInsertNsj", str(params.limitSjdbInsertNsj),
                "--outFileNamePrefix", params.outPrefix
            ]
            if params.twopassMode:
                cmd1.extend(["--twopassMode", params.twopassMode])
            if params.outTmpDir:
                shutil.rmtree(params.outTmpDir, ignore_errors=True)
                cmd1.extend(["--outTmpDir", params.outTmpDir])
            if params.outReadsUnmapped:
                cmd1.extend(["--outReadsUnmapped", params.outReadsUnmapped])
            if params.outSAMstrandField:
                cmd1.extend(["--outSAMstrandField", params.outSAMstrandField])
            read_files = params.input_params.split()
            cmd1.extend(["--readFilesIn"] + read_files)
            cmd2 = ["mv", f"{params.outPrefix}Aligned.sortedByCoord.out.bam", output.bam]
            cmd3 = [params.SAMTOOLS, "index", "-@", str(threads), output.bam]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                # Touch unmapped output files if STAR didn't produce them
                f.write(f"test -f {output.unmapped_r1} || touch {output.unmapped_r1}\n")
                f.write(f"test -f {output.unmapped_r2} || touch {output.unmapped_r2}\n")
                f.write(f"echo 'STAR alignment for {wildcards.sample_id} completed successfully.'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during STAR alignment for sample {wildcards.sample_id}: {str(e)}")
            raise RuntimeError(f"Error occurred during STAR alignment for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

rule star_result:
    input:
        star_align = outdir + "/{sample_id}/{sample_id}.bam"

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "star_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
