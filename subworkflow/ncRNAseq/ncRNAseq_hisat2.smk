
rule hisat2_index_ncRNAseq:
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output:
        index = expand(
            outdir + "/genome/{{genome}}/index/hista2/{{genome}}.{idx}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    threads: 8
    conda:
        config['conda']['run']
    params:
        prefix = lambda wildcards: outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}",
        HISAT2_BUILD = config.get('Procedure',{}).get('hisat2-build') or 'hisat2-build'
    log:
        outdir + "/log/genome/{genome}/hisat2_build.log"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        {params.HISAT2_BUILD} -p {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """

def get_hisat2_index_ncRNAseq(wildcards):
    logging.info(f"[get_hisat2_index] called with wildcards: {wildcards}")
    config_index_prefix = config.get('genome',{}).get(wildcards.genome,{}).get('hisat2_index_prefx') or None
    if config_index_prefix:
        first_file = f"{config_index_prefix}.1.ht2"
        if os.path.exists(first_file):
            return [f"{config_index_prefix}.{idx}.ht2" for idx in [1, 2, 3, 4, 5, 6, 7, 8]]
    return expand(
        outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}.{{idx}}.ht2",
        idx = [1, 2, 3, 4, 5, 6, 7, 8]
    )

rule hisat2_align_ncRNAseq_single:
    """
    HISAT2 alignment for single-end ncRNA-seq reads (STAR-equivalent).
    --no-spliced-alignment : disable spliced alignment
    -k 99999 : allow multiple mapping for small RNA
    --score-min L,0,-0.6 : at least 2/3 of read length should match
    """
    input:
        fastq = outdir + "/ncRNAseq/cutadapt/{sample_id}_cutadapt2_trimmed.fq.gz",
        genome_index = get_hisat2_index_ncRNAseq
    output:
        bam = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam",
        unmapped = outdir + "/ncRNAseq/bam/{genome}/{sample_id}.unmapped.fq.gz"
    log:
        outdir + "/log/ncRNAseq/HISAT2/{genome}/{sample_id}.log"
    threads: 12
    params:
        HISAT2 = config.get('Procedure',{}).get('hisat2') or 'hisat2',
        SAMTOOLS = config.get('Procedure',{}).get('samtools') or 'samtools',
    conda:
        config['conda']['smallRNAseq']
    shell:
        """
        {params.HISAT2} \
            -x {input.genome_index} \
            -U {input.fastq} \
            -p {threads} \
            --no-spliced-alignment \
            -k 99999 \
            --score-min L,0,-0.6 \
            --un-gz {output.unmapped} \
            2> {log} \
        | {params.SAMTOOLS} sort -@ {threads} -o {output.bam}
        """


