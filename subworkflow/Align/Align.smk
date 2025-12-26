SNAKEFILE_FULL_PATH_Align = workflow.snakefile
SNAKEFILE_DIR_Align = os.path.dirname(SNAKEFILE_FULL_PATH_Align)
alignYaml = get_yaml_path("Align",SNAKEFILE_DIR_Align)
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
        trim_galore = config.get('Procedure',{}).get('trim_galore') or 'trim_galore'
    threads: 6
    conda:
        config['conda']['run']
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
        trim_galore = config.get('Procedure',{}).get('trim_galore') or 'trim_galore'
    threads: 6
    conda:
        config['conda']['run']
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


rule alignment_result:
    input:
        outdir + "/Align/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"


if config["Procedure"]["aligner"] == "star":
    include: "align_star.smk"
    logging.info("aligner: star, load align_star.smk")

elif config["Procedure"]["aligner"] == "hisat2":
    include: "align_hisat2.smk"
    logging.info("aligner: hisat2, load align_hisat2.smk")
else:
    # 默认使用star比对
    include: "align_star.smk"
    logging.info("default: load align_star.smk")