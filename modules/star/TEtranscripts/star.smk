import logging
logger = logging.getLogger(__name__)
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "input")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])


def get_star_index(wildcards):
    logger.info(f"[get_star_index] called with wildcards: {wildcards}")
    star_index_dir = config.get('genome',{}).get(wildcards.genome,{}).get('star_index_dir') or None
    if star_index_dir:
        first_file = os.path.join(star_index_dir, "Genome")
        if os.path.exists(first_file):
            return star_index_dir

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
    paired_r1 = f"{indir}/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{indir}/{wildcards.sample_id}_2.fq.gz"
    single = f"{indir}/{wildcards.sample_id}.fq.gz"

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
rule TEtranscript_prepare_star:
    """
    STAR alignment to do TE analysis.
    --outFilterMultimapNmax 100 : allow each read to map to at most 100 genomic locations
    --winAnchorMultimapNmax 100 : allow each alignment seed (anchor) to map to at most 100 locations
    """
    input:
        get_alignment_input,
        genome_index = get_star_index
    output:
        outfile = temp(outdir + "/{sample_id}.bam")
    log:
        log = logdir + "/{sample_id}/star_align.log",
        STAR_log = logdir + "/{sample_id}Log.out",
        STAR_progress = logdir + "/{sample_id}Log.progress.out",
        STAR_final = logdir + "/{sample_id}Log.final.out"
    threads: 12
    params:
        outPrefix = outdir + "/{sample_id}.",
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        # 动态判断输入参数,加上genome_index，如果三个参数，即为双端测序，两个参数即为单端测序
        input_params = lambda wildcards, input: \
            f"{input[0]} {input[1]}" if len(input) == 3 else f"{input[0]}"
    shell:
        """
        # Allows for a large number of multi-maps 
        {params.STAR} --runThreadN {threads} \
            --genomeDir {input.genome_index} \
            --readFilesCommand zcat \
            --readFilesIn {params.input_params} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.outPrefix} \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 100  > {log.log} 2>&1
        cp {params.outPrefix}Log.out {log.STAR_log}
        cp {params.outPrefix}Log.progress.out {log.STAR_progress}
        cp {params.outPrefix}Log.final.out {log.STAR_final}
        """
