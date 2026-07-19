include: "../common/common.smk"

import logging
import time
logger = logging.getLogger(__name__)
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "input")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
fasta = config.get('genome',{}).get('fasta')
gtf = config.get('genome',{}).get('gtf')
fastq_sample_suffix = config.get('fastq_sample_suffix') or None
rule star_index:
    input:
        fasta = fasta,
    output:
        index_file = directory(outdir + "/index")
    log:
        logdir + "/index/star_index.log"
    threads: 12
    conda:
        "star.yaml"
    params:
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        index_dir = outdir + "/index",
        sjdbOverhang = config.get('Params',{}).get('STAR', {}).get('sjdbOverhang') or 100,
        gtf = gtf
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("star_index", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_index at {current_time}")

            sample_outdir = os.path.dirname(str(output.index_file))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"star_index_{current_time}.sh")

            cmd = [
                params.STAR, "--runMode", "genomeGenerate",
                "--runThreadN", str(threads),
                "--genomeDir", params.index_dir,
                "--genomeFastaFiles", input.fasta,
                "--sjdbOverhang", str(params.sjdbOverhang),
            ]
            if params.gtf:
                cmd.extend(["--sjdbGTFfile", params.gtf])
                rule_logger.info(f"Using sjdbGTFfile: {params.gtf}")
            else:
                rule_logger.info("No GTF provided, skipping --sjdbGTFfile")

            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} > {log_path} 2>&1")

            rule_logger.info(f"star_index completed successfully")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_index: {e}\n")
            logger.error(f"Error occurred during star_index: {e}")
            raise e
def get_star_index(wildcards):
    logger.info(f"[get_star_index] called with wildcards: {wildcards}")
    star_index_dir = config.get('genome',{}).get('index_dir') or None
    if star_index_dir:
        logger.info(f"[get_star_index] using provided star_index_dir: {star_index_dir}")
        first_file = os.path.join(star_index_dir, "Genome")
        if os.path.exists(first_file):
            return star_index_dir
    logger.info(f"[get_star_index] using default star_index_dir")
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
    logger.info(f"[get_alignment_input] called with wildcards: {wildcards}")
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
        logger.info(f"双端测序：{[paired_r1, paired_r2]}")
        logger.info(f"双端测序：{[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        logger.info(f"单端测序：{[single]}")
        logger.info(f"单端测序：{[single]}")
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
        bai = outdir + "/{sample_id}/{sample_id}.bam.bai"
    log:
        logdir + "/{sample_id}/star_align.log"
    threads: 12
    params:
        outPrefix = outdir + "/{sample_id}/{sample_id}.",
        input_params = lambda wildcards, input: \
            f"{input.fastq[0]} {input.fastq[1]}" if len(input.fastq) == 2 else f"{input.fastq[0]}",
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        SAMTOOLS = config.get('Procedure',{}).get('samtools') or 'samtools',
        alignEndsType = config.get('Params',{}).get('STAR', {}).get('alignEndsType') or "Local",
        outFilterMismatchNoverReadLmax = config.get('Params',{}).get('STAR', {}).get('outFilterMismatchNoverReadLmax') or 1.0,
        outFilterMismatchNmax = config.get('Params',{}).get('STAR', {}).get('outFilterMismatchNmax') or 10,
        outFilterMultimapNmax = config.get('Params',{}).get('STAR',{}).get('outFilterMultimapNmax') or 10,
        winAnchorMultimapNmax = config.get('Params',{}).get('STAR', {}).get('winAnchorMultimapNmax') or 50,
        genomeLoad = config.get('Params',{}).get('STAR', {}).get('genomeLoad') or 'NoSharedMemory',
        limitBAMsortRAM = config.get('Params',{}).get('STAR', {}).get('limitBAMsortRAM') or 0,
        outReadsUnmapped = config.get('Params',{}).get('STAR', {}).get('outReadsUnmapped') or None,
        outFilterMismatchNoverLmax = config.get('Params',{}).get('STAR', {}).get('outFilterMismatchNoverLmax') or 0.3,
        outFilterMatchNminOverLread = config.get('Params',{}).get('STAR', {}).get('outFilterMatchNminOverLread') or 0.66,
        alignSJoverhangMin = config.get('Params',{}).get('STAR', {}).get('alignSJoverhangMin') or 5,
        alignSJDBoverhangMin = config.get('Params',{}).get('STAR', {}).get('alignSJDBoverhangMin') or 3,
        chimSegmentMin = config.get('Params',{}).get('STAR', {}).get('chimSegmentMin') or 0,
        chimOutType = config.get('Params',{}).get('STAR', {}).get('chimOutType') or "Junctions",
        chimJunctionOverhangMin = config.get('Params',{}).get('STAR', {}).get('chimJunctionOverhangMin') or 20,
        outSAMstrandField = config.get('Params',{}).get('STAR', {}).get('outSAMstrandField') or None,
        chimScoreMin = config.get('Params',{}).get('STAR', {}).get('chimScoreMin') or 0,
        chimScoreDropMax = config.get('Params',{}).get('STAR', {}).get('chimScoreDropMax') or 20,
        chimScoreJunctionNonGTAG = config.get('Params',{}).get('STAR', {}).get('chimScoreJunctionNonGTAG') or -1,
        chimScoreSeparation = config.get('Params',{}).get('STAR', {}).get('chimScoreSeparation') or 10,
        alignSJstitchMismatchNmax = config.get('Params',{}).get('STAR', {}).get('alignSJstitchMismatchNmax') or "0 -1 0 0",
        chimSegmentReadGapMax = config.get('Params',{}).get('STAR', {}).get('chimSegmentReadGapMax') or 0
    conda:
        "star.yaml"
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/{wildcards.sample_id}/star_align.{current_time}.sh"
        cmd1 = [
            params.STAR, "--runThreadN", str(threads),
            "--genomeDir", input.genome_index,
            "--twopassMode", "Basic",
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
            "--outSAMattributes", "NM",
            "--outFileNamePrefix", params.outPrefix
        ]
        if params.outReadsUnmapped:
            cmd1.extend(["--outReadsUnmapped", params.outReadsUnmapped])
        if params.outSAMstrandField:
            cmd1.extend(["--outSAMstrandField", params.outSAMstrandField])
        read_files = params.input_params.split()
        cmd1.extend(["--readFilesIn"] + read_files)
        cmd2 = ["mv", f"{params.outPrefix}Aligned.sortedByCoord.out.bam", output.bam]
        cmd3 = [params.SAMTOOLS, "index", "-@", str(threads), output.bam]
        with open(script, 'w') as f:
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(" ".join(cmd3) + "\n")
        shell(f"bash {script} > {log} 2>&1")

rule star_result:
    input:
        star_align = outdir + "/{sample_id}/{sample_id}.bam"
