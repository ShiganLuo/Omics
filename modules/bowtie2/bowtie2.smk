from snakemake.logging import logger
import time
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "input")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
BOWTIE2_IDX_SUFFIX = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
rule bowtie2_index:
    input:
        fasta = config.get('genome',{}).get('fasta')
    output:
        index = expand(
            outdir + "/index/genome.{ext}",
            ext = BOWTIE2_IDX_SUFFIX
        )
    log:
        logdir + "/bowtie2_index.log"
    threads: 12
    params:
        bowtie2_build = config.get('Procedure',{}).get('bowtie2-build') or 'bowtie2-build',
        index_prefix = outdir + "/index/genome"
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/index/bowtie2_index.{current_time}.sh"
        cmd = [params.bowtie2_build, 
                "--threads", str(threads), 
                input.fasta, 
                params.index_prefix
                ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")


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
    paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
    single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}single.fq.gz"
    
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
def get_bowtie2_index(wildcards):
    logger.info(f"[get_bowtie2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('genome',{}).get("index_prefix") or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.1.bt2"
        if os.path.exists(first_file):
            return [f"{config_index_prefix}.{ext}" for ext in BOWTIE2_IDX_SUFFIX]
        else:
            logger.info(f"Config index prefix {config_index_prefix} provided, but does not exist. Falling back to default index path.")
    else:
        logger.info(f"No config index prefix provided. Using default index path.")
    return [f"{outdir}/index/genome.{ext}" for ext in BOWTIE2_IDX_SUFFIX]

rule bowtie2_align_paired:
    input:
        fastq1 = indir + "/{sample_id}/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}/{sample_id}_2.fq.gz",
        index = get_bowtie2_index
    output:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        metrics = outdir + "/{sample_id}/{sample_id}_bowtie2_metrics.txt",
        unmapped1 = outdir + "/{sample_id}/{sample_id}_unmapped_1.fq.gz",
        unmapped2 = outdir + "/{sample_id}/{sample_id}_unmapped_2.fq.gz"
    log:
        logdir + "/{sample_id}/bowtie2_align.log"
    threads: 8
    params:
        bowtie2 = config.get('Procedure',{}).get('bowtie2') or 'bowtie2',
        index_prefix = lambda wildcards, input: input.index[0].split(".")[0],
        unmapped_prefix = outdir + "/{sample_id}/{sample_id}_unmapped",
        sam_append_comment = config.get('Params',{}).get('bowtie2', {}).get('sam-append-comment') or False
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/{wildcards.sample_id}/bowtie2_align.{current_time}.sh"
        cmd1 = [params.bowtie2, 
                "-x", params.index_prefix, 
                "-1", input.fastq1, 
                "-2", input.fastq2,
                "--threads", str(threads),
                "--sam-append-comment" if params.sam_append_comment else "",
                "--met-file", output.metrics,
                "--un-conc-gz", params.unmapped_prefix,
                "|", "samtools", "view", "-@", str(threads), "-bS", "-", ">", output.bam
                ]
        cmd2 = ["mv", f"{params.unmapped_prefix}.1", output.unmapped1]
        cmd3 = ["mv", f"{params.unmapped_prefix}.2", output.unmapped2]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(" ".join(cmd3) + "\n")
        shell("bash {script} > {log} 2>&1")

rule bowtie2_align_single:
    input:
        fastq = indir + "/{sample_id}/{sample_id}.single.fq.gz",
        index = get_bowtie2_index
    output:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        metrics = outdir + "/{sample_id}/{sample_id}_bowtie2_metrics.txt",
        unmapped = outdir + "/{sample_id}/{sample_id}_unmapped.single.fq.gz"
    log:
        logdir + "/{sample_id}/bowtie2_align.log"
    threads: 8
    params:
        bowtie2 = config.get('Procedure',{}).get('bowtie2') or 'bowtie2',
        unmapped_prefix = outdir + "/{sample_id}/{sample_id}_unmapped",
        index_prefix = lambda wildcards, input: input.index[0].split(".")[0],
        sam_append_comment = config.get('Params',{}).get('bowtie2', {}).get('sam-append-comment') or False
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/{wildcards.sample_id}/bowtie2_align.{current_time}.sh"
        cmd1 = [params.bowtie2, 
                "-x", params.index_prefix, 
                "-U", input.fastq, 
                "--threads", str(threads),
                "--sam-append-comment" if params.sam_append_comment else "",
                "--met-file", output.metrics,
                "--un-gz", params.unmapped_prefix,
                "|", "samtools", "view", "-@", str(threads), "-bS", "-", ">", output.bam
                ]
        cmd2 = ["mv", params.unmapped_prefix, output.unmapped]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
        shell("bash {script} > {log} 2>&1")