import logging
logger = logging.getLogger(__name__)
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
fasta = config.get('genome',{}).get('fasta')
paired_samples = config.get('paired_samples', [])
single_samples = config.get('single_samples', [])
# wait for test, remove after test
rule hisat2_index:
    input:
        fasta = fasta
    output:
        index = expand(
            outdir + "/index/genome.{idx}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    threads: 8
    conda:
        "hisat2.yaml"
    params:
        prefix = outdir + "/index/genome",
        HISAT2_BUILD = config.get('Procedure',{}).get('hisat2-build') or 'hisat2-build'
    log:
        logdir + "/index/hisat2_build.log"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        {params.HISAT2_BUILD} -p {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """

def get_hisat2_index(wildcards):
    logger.info(f"[get_hisat2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('genome',{}).get('index_prefix') or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.1.ht2"
        if os.path.exists(first_file):
            return [f"{config_index_prefix}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]
    return [outdir + f"/index/genome.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]


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

rule hisat2_align:
    input:
        fastq = get_alignment_input,
        index = get_hisat2_index
    output:
        outfile = outdir + "/{sample_id}/{sample_id}.bam",
        splice = outdir + "/{sample_id}/splice.txt"
    log:
        logdir + "/{sample_id}/hisat2_align.log"
    threads: 12
    conda:
        "hisat2.yaml"
    params:
        HISAT2 = config.get('Procedure',{}).get('hisat2') or 'hisat2',
        SAMTOOLS = config.get('Procedure',{}).get('samtools') or 'samtools',
        index_prefix = lambda wildcards, input: input.index[0].rsplit('.', 2)[0],
        unmapped_prefix = lambda wildcards: f"{outdir}/{wildcards.sample_id}/unmapped",
        input_params = lambda wildcards, input: \
            f"-1 {input.fastq[0]} -2 {input.fastq[1]}" if len(input.fastq) == 2 else f"-U {input.fastq[0]}"
    shell:
        """
        mkdir -p $(dirname {output.outfile})
        {params.HISAT2} -x {params.index_prefix} \
            -k 1 \
            -p {threads} \
            -q \
            --no-unal \
            --dta \
            --novel-splicesite-outfile {output.splice} \
            --un-conc-gz {params.unmapped_prefix} \
            {params.input_params} \
            2> {log} | \
        {params.SAMTOOLS} sort -@ {threads} -@ {threads} -o {output.outfile}
        """


rule hisat2_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.bam"
