SNAKEFILE_FULL_PATH_Align = workflow.snakefile
SNAKEFILE_DIR_Align = os.path.dirname(SNAKEFILE_FULL_PATH_Align)
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_Align ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
alignYaml = get_yaml_path("Align")
configfile: alignYaml
logging.info(f"Include Align config: {alignYaml}")
rule trimming_Paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fastq.gz",
        fastq2 = indir + "/{sample_id}_2.fastq.gz"
    output:
        fastq1 = temp(outdir + "/cutadapt/{sample_id}_1.fq.gz"),
        fastq2 = temp(outdir + "/cutadapt/{sample_id}_2.fq.gz"),
        report1 = outdir + "/log/Align/{sample_id}/trimming_statistics_1.txt",
        report2 = outdir + "/log/Align/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir = outdir + "/cutadapt",
        quality = 30,
        trim_galore = config['Procedure']['trim_galore']
    threads: 6
    log:
        log = outdir + "/log/Align/{sample_id}/trimming.txt"
    shell:
        """
        # trim_galore can automatically judge the fq quality scoring system,it's no need to add such as --phred33 --phred64
        {params.trim_galore} --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fastq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fastq.gz_trimming_report.txt {output.report2}
        """
rule trimming_Single:
    input:
        fastq = indir + "/{sample_id}.fastq.gz"
    output:
        fastq = temp(outdir + "/cutadapt/{sample_id}Single.fq.gz"),
        report = outdir + "/log/Align/{sample_id}/trimming_statistics.txt"
    params:
        outdir = outdir + "/cutadapt",
        quality = 30,
        trim_galore = config['Procedure']['trim_galore']
    threads: 6
    log:
        log = outdir + "/log/Align/{sample_id}/trimming.txt"
    shell:
        """
        {params.trim_galore} --phred33  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_trimmed.fq.gz {output.fastq}
        mv {params.outdir}/{wildcards.sample_id}.fastq.gz_trimming_report.txt {output.report}
        """

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
    logging.info(f"[get_alignment_input] called with wildcards: {wildcards}")
    # 构造可能的输入路径
    paired_r1 = f"{outdir}/cutadapt/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{outdir}/cutadapt/{wildcards.sample_id}_2.fq.gz"
    single = f"{outdir}/cutadapt/{wildcards.sample_id}Single.fq.gz"
    
    # 检查文件实际存在情况
    if wildcards.sample_id in paired_samples:
        logging.info(f"双端测序：{[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        logging.info(f"单端测序：{[single]}")
        return [single]
    else:
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample_id}\n"
            f"Checked paths:\n- {paired_r1}\n- {paired_r2}\n- {single}"
        )

rule star_align:
    input:
        fastq = get_alignment_input,
        genome_index = lambda wildcards: config['genome'][wildcards.genome]['genome_index']
    output:
        outfile = outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        outdir + "/log/Align/{sample_id}/{genome}/star_align.log"
    threads: 25
    params:
        outPrefix = outdir + "/2pass/{sample_id}/{genome}/{sample_id}",
        # 动态判断输入参数,加上genome_index，如果三个参数，即为双端测序，两个参数即为单端测序
        input_params = lambda wildcards, input: \
            f"{input[0]} {input[1]}" if len(input) == 3 else f"{input[0]}",
        STAR = config['Procedure']['STAR']
    shell:
        """
        mkdir -p $(dirname {params.outPrefix})
        {params.STAR} --runThreadN {threads} \
            --genomeDir {input.genome_index} \
            --twopassMode Basic \
            --readFilesCommand zcat \
            --readFilesIn {params.input_params} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NM \
            --outFileNamePrefix {params.outPrefix} > {log} 2>&1
        """
def get_bams_for_featureCounts_single(wildcards):
    bams = []
    for sample_id, genome in single_sample_genome_pairs:
        if genome == wildcards.genome:
            bams.append(f"{outdir}/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_single_noMultiple didn't get any input bams, genome: {wildcards.genome},\nsingle_sample_genome_pairs: {single_sample_genome_pairs}")
    return bams

rule featureCounts_single_noMultiple:
    input:
        bams = get_bams_for_featureCounts_single
    output:
        outfile = outdir + "/counts/featureCounts/{genome}/{genome}_single_count.tsv"
    log:
        outdir + "/log/Align/{genome}_featureCounts_single_noMultiple.log"
    threads:
        20
    params:
        featureCounts = config['Procedure']['featureCounts'],
        gtf = lambda wildcards: config["genome"][wildcards.genome]["gtf"]
    shell:
        """
        {params.featureCounts} -T {threads} -t exon -g gene_id -a {params.gtf} -o {output.outfile} {input.bams} > {log} 2>&1
        """

def get_bams_for_featureCounts_paired(wildcards):
    bams = []
    for sample_id, genome in paired_sample_genome_pairs:
        if genome == wildcards.genome:
            bams.append(f"{outdir}/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_paired_noMultiple didn't get any input bams, genome: {wildcards.genome},\npaired_sample_genome_pairs:{paired_sample_genome_pairs}")
    return bams

rule featureCounts_paired_noMultiple:
    input:
        bams = get_bams_for_featureCounts_paired
    output:
        outfile = outdir + "/counts/featureCounts/{genome}/{genome}_paired_count.tsv"
    log:
        outdir + "/log/Align/{genome}_featureCounts_paired_noMultiple.log"
    threads:
        20
    params:
        featureCounts = config['Procedure']['featureCounts'],
        gtf = lambda wildcards: config["genome"][wildcards.genome]["gtf"]
    shell:
        """
        # for multiple -M -O
        {params.featureCounts} -T {threads} -B -p --countReadPairs -t exon -g gene_id -a {params.gtf} -o {output.outfile} {input.bams} > {log} 2>&1
        """
