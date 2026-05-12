import os
indir = config.get("indir","input")
outdir = config.get("outdir","output")
logdir = config.get("logdir","log")
fasta = config.get("fasta")
paired_samples = config.get("paired_samples")
single_samples = config.get("single_samples")
#Mapping with bwa
rule bwaMem2_index:
    input:
        fa = fasta
    output:
        index = expand(
            outdir + "/genome/{{genome}}/index/{{genome}}.{ext}",
            ext = BWAMEM2_IDX_SUFFIX
        )
    log:
        logdir + "/index/{genome}/bwa_mem2-indexing.txt"
    threads: 15
    conda:
        "../bwa-mem2.yaml"
    params:
        bwa_mem2 = config.get("Procedure",{}).get("bwaMem2") or "bwa-mem2",
        index_prefix = lambda wildcards: outdir + f"/genome/{wildcards.genome}/index/{wildcards.genome}"
    shell:
        """
        mkdir -p $(dirname {params.index_prefix})
        {params.bwa_mem2} index -p {params.index_prefix} {input.fa} > {log} 2>&1
        """
def get_bwaMem2_index(wildcards):
    print(f"[get_bwaMem2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('bwaMem2_index_prefix') or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.123"
        if os.path.exists(first_file):
            return expand(
                config_index_prefix + ".{ext}",
                ext = BWAMEM2_IDX_SUFFIX
            )
        else:
            print(f"Config index prefix {config_index_prefix} provided, but does not exist. Falling back to default index path.")
    else:
        print(f"No config index prefix provided. Using default index path.")
    return expand(
        outdir + f"/genome/{wildcards.genome}/index/{wildcards.genome}.{{ext}}",
        ext = BWAMEM2_IDX_SUFFIX
    )

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
    paired_r1 = f"{indir}/cutadapt/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{indir}/cutadapt/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
    single = f"{indir}/cutadapt/{wildcards.sample_id}/{wildcards.sample_id}single.fq.gz"
    
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
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample_id}\n"
            f"Checked paths:\n- {paired_r1}\n- {paired_r2}\n- {single}"
        )


rule bwaMem2_alignment:
    input:
        fastq = get_alignment_input,
        index = get_bwaMem2_index
    output:
        bam = outdir + "/bam/{genome}/{sample_id}.bam",
    log:
        logdir + "/log/bwa-mem2/{genome}/{sample_id}/bwa-alignment.txt"
    threads: 15
    conda:
        "../bwa-mem2.yaml"
    params:
        bwa_mem2 = config.get("Procedure",{}).get("bwaMem2") or "bwa-mem2",
        samtools = config.get("Procedure",{}).get("samtools") or "samtools",
        index_prefix = lambda wildcards, input: input.index[0].rsplit(".", 1)[0],
        input_params = lambda wildcards, input: \
            f"{input.fastq[0]} {input.fastq[1]}" if len(input.fastq) == 2 else f"{input.fastq[0]}"
    shell:
        """
        {params.bwa_mem2} mem \
        -T 0 \
        -t {threads} \
        {params.index_prefix} \
        {params.input_params} \
        2>> {log} \
        | {params.samtools} view -b - 2>> {log} \
        > {output.bam}

        """
rule bwaMemm2_result:
    input:
        bam = outdir + "/bam/{genome}/{sample_id}.bam"