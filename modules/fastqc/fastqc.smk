import os
indir = config.get("indir", "output/raw_fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
log_suffix = config.get("log_suffix", "txt")
# need test
def get_fastqc_input(wildcards):
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
    logger.info(f"[get_fastqc_input] called with wildcards: {wildcards}")
    # 构造可能的输入路径
    paired_r1 = f"{indir}/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{indir}/{wildcards.sample_id}_2.fq.gz"
    single = f"{indir}/{wildcards.sample_id}.single.fq.gz"

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


rule fastqc:
    input:
        get_fastqc_input
    output:
        outdir = directory(outdir + "/{sample_id}"),
        flag = outdir + "/{sample_id}/fastqc." + log_suffix
    params:
        fastqc = config.get("Procedure", {}).get("fastqc") or "fastqc"
    threads: 2
    log:
        log = logdir + "/{sample_id}/fastqc." + log_suffix
    conda:
        "fastqc.yaml"
    shell:
        """
        {params.fastqc} \
            --threads {threads} \
            -o {output.outdir} \
            -t {threads} \
            {input} \
            > {log.log} 2>&1
        touch {output.flag}
        """

