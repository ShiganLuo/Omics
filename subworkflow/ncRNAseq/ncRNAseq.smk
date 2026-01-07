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

rule star_index_ncRNAseq:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta'],
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    output:
        index_dir = directory(outdir + "/genome/{genome}/index/star")
    log:
        outdir + "/log/genome/{genome}/star_index.log"
    threads: 12
    conda:
        config['conda']['run']
    params:
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        # 索引目录路径
        
    shell:
        """
        mkdir -p {output.index_dir}
        {params.STAR} --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100 > {log} 2>&1
        """

def get_star_index_ncRNAseq(wildcards):
    logging.info(f"[get_star_index_ncRNAseq] called with wildcards: {wildcards}")
    star_index_dir = config.get('genome',{}).get(wildcards.genome,{}).get('star_index_dir') or None
    if star_index_dir:
        logging.info(f"[get_star_index_ncRNAseq] using provided star_index_dir: {star_index_dir}")
        first_file = os.path.join(star_index_dir, "Genome")
        if os.path.exists(first_file):
            return star_index_dir
    default_index_dir = outdir + f"/genome/{wildcards.genome}/index/star"
    logging.info(f"[get_star_index_ncRNAseq] using default star_index_dir: {default_index_dir}")
    return default_index_dir

rule star_align_ncRNAseq_single:
    """
    STAR alignment for single-end ncRNA-seq reads. # normal sample, not virus
    """
    input:
        fastq = outdir + "/ncRNAseq/cutadapt/{sample_id}_cutadapt2_trimmed.fq.gz",
        genome_index = get_star_index_ncRNAseq
    output:
        bam = outdir + "/ncRNAseq/star/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam",
        unmapped = outdir + "/ncRNAseq/star/{genome}/{sample_id}.Unmapped.out.mate1"
    log:
        outdir + "/log/ncRNAseq/STAR/{genome}/{sample_id}.log"
    threads: 14
    params:
        prefix = outdir + "/ncRNAseq/star/{genome}/{sample_id}.",
        STAR = config.get('Procedure', {}).get('STAR', "STAR") or "STAR"
    conda:
        config['conda']['smallRNAseq']
    shell:
        """
        {params.STAR} \
            --readFilesIn {input.fastq} \
            --runThreadN {threads} \
            --genomeDir {input.genome_index} \
            --genomeLoad LoadAndRemove \
            --readFilesCommand unpigz -c \
            --limitBAMsortRAM 20000000000 \
            --outFileNamePrefix {params.prefix} \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 99 \
            --outFilterMismatchNoverLmax 0.1 \
            --outFilterMatchNminOverLread 0.66 \
            --alignSJoverhangMin 999 \
            --alignSJDBoverhangMin 999 \
            > {log} 2>&1
        """

rule ncRNAseq_result:
    input:
        bam = outdir + "/ncRNAseq/star/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam"

