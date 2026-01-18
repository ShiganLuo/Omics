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
    logger.info(f"[get_star_index_ncRNAseq] called with wildcards: {wildcards}")
    star_index_dir = config.get('genome',{}).get(wildcards.genome,{}).get('star_index_dir') or None
    if star_index_dir:
        logger.info(f"[get_star_index_ncRNAseq] using provided star_index_dir: {star_index_dir}")
        first_file = os.path.join(star_index_dir, "Genome")
        if os.path.exists(first_file):
            return star_index_dir
    default_index_dir = outdir + f"/genome/{wildcards.genome}/index/star"
    logger.info(f"[get_star_index_ncRNAseq] using default star_index_dir: {default_index_dir}")
    return default_index_dir

rule star_align_ncRNAseq_single:
    """
    STAR alignment for single-end ncRNA-seq reads. # normal sample, not virus
    --outFilterMultimapNmax 99999 : allow multiple mapping for small RNA
    --outFilterMatchNminOverLread 0.66 : at least 2/3 of read length should match
    --alignSJoverhangMin 999 : disable spliced alignment
    --alignSJDBoverhangMin 999 : disable spliced alignment
    --outFilterMismatchNoverLmax 0.1 : at most 10% mismatches
    """
    input:
        fastq = outdir + "/ncRNAseq/cutadapt/{sample_id}_cutadapt2_trimmed.fq.gz",
        genome_index = get_star_index_ncRNAseq
    output:
        bam = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam",
        unmapped = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.Unmapped.out.mate1"
    log:
        outdir + "/log/ncRNAseq/STAR/{genome}/{sample_id}.log"
    threads: 14
    params:
        prefix = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.",
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
            --outFilterMultimapNmax 99999 \
            --outFilterMismatchNoverLmax 0.1 \
            --outFilterMatchNminOverLread 0.66 \
            --alignSJoverhangMin 999 \
            --alignSJDBoverhangMin 999 \
            > {log} 2>&1
        """
