import logging
logger = logging.getLogger("align_star")
rule star_index:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta'],
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    output:
        # 以索引目录下的核心文件作为标记
        index_file = outdir + "/genome/{genome}/index/star/Genome"
    log:
        outdir + "/log/genome/{genome}/star_index.log"
    threads: 12
    conda:
        config['conda']['run']
    params:
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR',
        # 索引目录路径
        index_dir = outdir + "/genome/{genome}/index/star"
    shell:
        """
        mkdir -p {params.index_dir}
        {params.STAR} --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.index_dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100 > {log} 2>&1
        """
def get_star_index(wildcards):
    logger.info(f"[get_star_index] called with wildcards: {wildcards}")
    star_index_dir = config.get('genome',{}).get(wildcards.genome,{}).get('star_index_dir') or None
    if star_index_dir:
        logger.info(f"[get_star_index] using provided star_index_dir: {star_index_dir}")
        first_file = os.path.join(star_index_dir, "Genome")
        if os.path.exists(first_file):
            return star_index_dir
    logger.info(f"[get_star_index] using default star_index_dir")
    return outdir + f"/genome/{wildcards.genome}/index/star"

rule star_align:
    input:
        fastq = get_alignment_input,
        genome_index = get_star_index
    output:
        outfile = outdir + "/Align/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam"
    log:
        outdir + "/log/Align/{sample_id}/{genome}/star_align.log"
    threads: 12
    conda:
        config['conda']['run']
    params:
        outPrefix = outdir + "/Align/{sample_id}/{genome}/{sample_id}.",
        input_params = lambda wildcards, input: \
            f"{input.fastq[0]} {input.fastq[1]}" if len(input.fastq) == 2 else f"{input.fastq[0]}",
        STAR = config.get('Procedure',{}).get('STAR') or 'STAR'
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