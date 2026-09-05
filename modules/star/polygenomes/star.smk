include: "../../common/common.smk"

import logging
import time
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_index = config.get("logdir_index", logdir)
log_index_substring = config.get("log_index_substring", "star_index")
indir= config.get("indir", "input")
genome_paired_samples = config.get("genome_paired_samples", {})
genome_single_samples = config.get("genome_single_samples", {})
omics_type = config.get('omics_type', 'RNA-seq')
fastq_sample_suffix = config.get('fastq_sample_suffix') or None
isGenomeSubdir = config.get('isGenomeSubdir', False)

def get_input_for_star_index(wildcards):
    """Dynamically determines the input fasta file for STAR index based on the genome."""
    logger.info(f"[get_input_for_star_index] called with wildcards: {wildcards}")
    fasta = config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('fasta')
    if not fasta:
        logger.error(f"Fasta file for genome {wildcards.genome} not found in config")
        raise ValueError(f"Fasta file for genome {wildcards.genome} not found in config")
    return fasta
rule star_index:
    input:
        fasta = get_input_for_star_index,
    output:
        index_file = directory(outdir + "/index/{genome}")
    log:
        logdir_index + "/star/{genome}/" + log_index_substring + ".log"
    threads: 12
    conda:
        "../star.yaml"
    container:
        sif("../star.yaml")
    params:
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        index_dir = lambda wildcards: outdir + f"/index/{wildcards.genome}",
        sjdbOverhang = config.get('Params',{}).get('star', {}).get('index', {}).get('sjdbOverhang') or 100,
        gtf = lambda wildcards: config.get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('gtf'),
        outTmpDir = lambda wildcards: outdir + f"/tmp_star_{wildcards.genome}"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("star_index", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start star_index for genome {wildcards.genome} at {current_time}")

            sample_outdir = os.path.dirname(str(output.index_file))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"star_index_{current_time}.sh")

            cmd = [
                params.STAR, "--runMode", "genomeGenerate",
                "--runThreadN", str(threads),
                "--genomeDir", params.index_dir,
                "--genomeFastaFiles", input.fasta]
            if params.outTmpDir:
                shutil.rmtree(params.outTmpDir, ignore_errors=True)
                cmd.extend(["--outTmpDir", params.outTmpDir])
            if params.gtf and os.path.exists(params.gtf):
                # STAR can't use --sjdbOverhang without an annotation file
                cmd.extend(["--sjdbGTFfile", str(params.gtf)])
                cmd.extend(["--sjdbOverhang", str(params.sjdbOverhang)])
                rule_logger.info(f"Using sjdbGTFfile: {params.gtf}")
            else:
                rule_logger.info("No GTF provided, skipping --sjdbGTFfile")

            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"star_index for genome {wildcards.genome} completed successfully")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during star_index for genome {wildcards.genome}: {e}\n")
            logger.error(f"Error occurred during star_index for genome {wildcards.genome}: {e}")
            raise e

def get_star_index(wildcards):
    logger.info(f"[get_star_index] called with wildcards: {wildcards}")
    config_index_dir = config.get('Params', {}).get('genome', {}).get('references', {}).get(wildcards.genome, {}).get('index_dir') or None
    if config_index_dir:
        logger.info(f"[get_star_index] using provided index_dir for genome {wildcards.genome}: {config_index_dir}")
        return config_index_dir
    logger.info(f"[get_star_index] using default index_dir for genome {wildcards.genome}")
    return outdir + f"/index/{wildcards.genome}"

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
        logger.info(f"[get_alignment_input] fastq_sample_suffix is set to {fastq_sample_suffix}")
        if not isGenomeSubdir:
            logger.info(f"[get_alignment_input] isGenomeSubdir is False. Assuming fastq files are in the sample directory.")
            paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}_1.fq.gz"
            paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}_2.fq.gz"
            single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}.single.fq.gz"
        else:
            logger.info(f"[get_alignment_input] isGenomeSubdir is True. Assuming fastq files are in the genome subdirectory.")
            paired_r1 = f"{indir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}_1.fq.gz"
            paired_r2 = f"{indir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}_2.fq.gz"
            single = f"{indir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}_{fastq_sample_suffix}.single.fq.gz"
    else:
        logger.info("[get_alignment_input] fastq_sample_suffix is not set. Using default naming convention.")
        if not isGenomeSubdir:
            logger.info(f"[get_alignment_input] isGenomeSubdir is False. Assuming fastq files are in the sample directory.")
            paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
            paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
            single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"
        else:
            logger.info(f"[get_alignment_input] isGenomeSubdir is True. Assuming fastq files are in the genome subdirectory.")
            paired_r1 = f"{indir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
            paired_r2 = f"{indir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
            single = f"{indir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

    # 检查文件实际存在情况
    if wildcards.sample_id in genome_paired_samples.get(wildcards.genome, []):
        logger.info(f"双端测序：{[paired_r1, paired_r2]}")
        if omics_type == "scRNAseq":
            logger.info(f"Detected scRNA-seq data for sample {wildcards.sample_id}.R1 is expected to contain cell barcodes and UMIs, R2 contains the transcript sequence.")
            return [paired_r2, paired_r1]
        else:   
            return [paired_r1, paired_r2]
    elif wildcards.sample_id in genome_single_samples.get(wildcards.genome, []):
        logger.info(f"单端测序：{[single]}")
        return [single]
    else:
        logger.error(f"样本 {wildcards.sample_id} 未在 {wildcards.genome} 的 paired_samples: {genome_paired_samples.get(wildcards.genome, [])} 或 single_samples: {genome_single_samples.get(wildcards.genome, [])} 中定义")
        raise ValueError(f"Sample {wildcards.sample_id} not defined in paired_samples or single_samples")

rule star_align:
    input:
        fastq = get_alignment_input,
        genome_index = get_star_index
    output:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.bam.bai",
        unmapped_r1 = outdir + "/{genome}/{sample_id}/{sample_id}.Unmapped.out.mate1",
        unmapped_r2 = outdir + "/{genome}/{sample_id}/{sample_id}.Unmapped.out.mate2",
    log:
        logdir + "/{sample_id}/{genome}/star_align.log"
    threads: 12
    params:
        outPrefix = outdir + "/{genome}/{sample_id}/{sample_id}.",
        input_params = lambda wildcards, input: \
            f"{input.fastq[0]} {input.fastq[1]}" if len(input.fastq) == 2 else f"{input.fastq[0]}",
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        SAMTOOLS = config.get('Procedure',{}).get('samtools') or 'samtools',
        alignEndsType = config.get('Params',{}).get('star', {}).get('alignEndsType') or "Local",
        outFilterMismatchNoverReadLmax = config.get('Params',{}).get('star', {}).get('outFilterMismatchNoverReadLmax') or 1.0,
        outFilterMismatchNmax = config.get('Params',{}).get('star', {}).get('outFilterMismatchNmax') or 10,
        outFilterMultimapNmax = config.get('Params',{}).get('star',{}).get('outFilterMultimapNmax') or 10,
        winAnchorMultimapNmax = config.get('Params',{}).get('star', {}).get('winAnchorMultimapNmax') or 50,
        genomeLoad = config.get('Params',{}).get('star', {}).get('genomeLoad') or 'NoSharedMemory',
        limitBAMsortRAM = config.get('Params',{}).get('star', {}).get('limitBAMsortRAM') or 0,
        outReadsUnmapped = config.get('Params',{}).get('star', {}).get('outReadsUnmapped') or None,
        outFilterMismatchNoverLmax = config.get('Params',{}).get('star', {}).get('outFilterMismatchNoverLmax') or 0.3,
        outFilterMatchNminOverLread = config.get('Params',{}).get('star', {}).get('outFilterMatchNminOverLread') or 0.66,
        alignSJoverhangMin = config.get('Params',{}).get('star', {}).get('alignSJoverhangMin') or 5,
        alignSJDBoverhangMin = config.get('Params',{}).get('star', {}).get('alignSJDBoverhangMin') or 3,
        chimSegmentMin = config.get('Params',{}).get('star', {}).get('chimSegmentMin') or 0,
        chimOutType = config.get('Params',{}).get('star', {}).get('chimOutType') or "Junctions",
        chimJunctionOverhangMin = config.get('Params',{}).get('star', {}).get('chimJunctionOverhangMin') or 20,
        outSAMstrandField = config.get('Params',{}).get('star', {}).get('outSAMstrandField') or None,
        chimScoreMin = config.get('Params',{}).get('star', {}).get('chimScoreMin') or 0,
        chimScoreDropMax = config.get('Params',{}).get('star', {}).get('chimScoreDropMax') or 20,
        chimScoreJunctionNonGTAG = config.get('Params',{}).get('star', {}).get('chimScoreJunctionNonGTAG') or -1,
        chimScoreSeparation = config.get('Params',{}).get('star', {}).get('chimScoreSeparation') or 10,
        alignSJstitchMismatchNmax = config.get('Params',{}).get('star', {}).get('alignSJstitchMismatchNmax') or "0 -1 0 0",
        chimSegmentReadGapMax = config.get('Params',{}).get('star', {}).get('chimSegmentReadGapMax') or 0,
        outSAMattributes = config.get('Params',{}).get('star', {}).get('outSAMattributes') or "NM",
        outMultimapperOrder = config.get('Params',{}).get('star', {}).get('outMultimapperOrder') or "Old_2.4",
        runRNGseed = config.get('Params',{}).get('star', {}).get('runRNGseed') or 777,
        outSAMmultNmax = config.get('Params',{}).get('star', {}).get('outSAMmultNmax') or -1,
        soloType = config.get('Params',{}).get('star', {}).get('soloType') or None,
        soloCBwhitelist = config.get('Params',{}).get('star', {}).get('soloCBwhitelist') or None,
        soloBarcodeReadLength = config.get('Params',{}).get('star', {}).get('soloBarcodeReadLength', 1), # may be 0
        limitSjdbInsertNsj = config.get('Params',{}).get('star', {}).get('limitSjdbInsertNsj') or 1000000,
        outTmpDir = outdir + "/{genome}/{sample_id}/tmp_star"
    conda:
        "../star.yaml"
    container:
        sif("../star.yaml")
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/{wildcards.genome}/{wildcards.sample_id}/star_align.{current_time}.sh"
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
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(" ".join(cmd3) + "\n")
            # Touch unmapped output files if STAR didn't produce them
            f.write(f"test -f {output.unmapped_r1} || touch {output.unmapped_r1}\n")
            f.write(f"test -f {output.unmapped_r2} || touch {output.unmapped_r2}\n")
        shell(f"bash {script} > {log} 2>&1")

rule star_result:
    input:
        star_align = outdir + "/{genome}/{sample_id}/{sample_id}.bam"
