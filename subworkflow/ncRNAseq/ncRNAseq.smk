SNAKEFILE_FULL_PATH_ncRNAseq = workflow.snakefile
SNAKEFILE_DIR_ncRNAseq = os.path.dirname(SNAKEFILE_FULL_PATH_ncRNAseq)
ncRNAseqYaml = get_yaml_path("ncRNAseq",SNAKEFILE_DIR_ncRNAseq)
configfile: ncRNAseqYaml 
logging.info(f"Include ncRNAseq config: {ncRNAseqYaml}")
# from https://doi.org/10.1093/nar/gkae1288
rule fastx_trimmer:
    """
    Trim the first N bases from the 5' end of reads using fastx_trimmer.
    N is specified by the parameter 'f_keep_base' in the configuration file.
    only for single-end reads.
    """
    input:
        fastq = indir + "/{sample_id}.fastq.gz",
    output:
        fastq_trimmed = outdir + "/ncRNAseq/fastx_trimmer/{sample_id}_fastx1_trimmed.fq.gz",
    log:
        outdir + "/log/ncRNAseq/Trim/fastx_trimmer/{sample_id}_fastx_trimmer.log"
    conda:
        config['conda']['smallRNAseq']
    params:
        fastx_trimmer = config.get('Procedure',{}).get('fastx_trimmer',{}) or "fastx_trimmer",
        f_keep_base = config.get('Parameters',{}).get('fastx_trimmer',{}).get('f_keep_base',{}) or 5
    shell:
        """
        zcat {input.fastq} \
        | {params.fastx_trimmer} -f {params.f_keep_base} -z \
          -o {output.fastq_trimmed} \
        > {log} 2>&1
        """

rule cutadapt_ncRNAseq_single:
    """
    Cutadapt (v1.8.3) for single-end ncRNA-seq reads after fastx_trimmer.
    """
    input:
        fastq_trimmed = outdir + "/ncRNAseq/fastx_trimmer/{sample_id}_fastx1_trimmed.fq.gz",
    output:
        cutadapt_fastq = outdir + "/ncRNAseq/cutadapt/{sample_id}_cutadapt2_trimmed.fq.gz",
    log:
        outdir + "/log/ncRNAseq/Trim/cutadapt/{sample_id}_cutadapt_ncRNAseq_single.log"
    conda:
        config['conda']['smallRNAseq']
    params:
        cutadapt = config.get('Procedure', {}).get('cutadapt', "cutadapt") or "cutadapt"
    shell:
        """
        {params.cutadapt} \
            -g GTTCAGAGTTCTACAGTCCGACGATCNNNN \
            -a NNNNTGGAATTCTCGGGTGCCAAGG \
            -e 0.075 \
            -n 2 \
            -O 14 \
            -m 12 \
            --max-n 0.1 \
            --trim-n \
            --match-read-wildcards \
            -o {output.cutadapt_fastq} \
            {input.fastq_trimmed} \
            > {log} 2>&1
        """



rule ncRNAseq_result:
    input:
        bam = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam"


if config["Procedure"]["aligner"] == "star":
    include: "ncRNAseq_star.smk"
    logging.info("aligner: star, load ncRNAseq_star.smk")

elif config["Procedure"]["aligner"] == "hisat2":
    include: "ncRNAseq_hisat2.smk"
    logging.info("aligner: hisat2, load ncRNAseq_hisat2.smk")
else:
    # 默认使用star比对
    include: "ncRNAseq_star.smk"
    logging.info("default: load ncRNAseq_star.smk")
