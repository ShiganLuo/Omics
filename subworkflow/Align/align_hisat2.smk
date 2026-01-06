rule hisat2_index:
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

def get_hisat2_index(wildcards):
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

rule hisat2_align:
    input:
        fastq = get_alignment_input,
        index = get_hisat2_index
    output:
        outfile = outdir + "/Align/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        outdir + "/log/Align/{sample_id}/{genome}/hisat2_align.log"
    threads: 12
    conda:
        config['conda']['run']
    params:
        HISAT2 = config.get('Procedure',{}).get('hisat2') or 'hisat2',
        SAMTOOLS = config.get('Procedure',{}).get('samtools') or 'samtools',
        index_prefix = lambda wildcards, input: input.index[0].rsplit('.', 2)[0],
        input_params = lambda wildcards, input: \
            f"-1 {input.fastq[0]} -2 {input.fastq[1]}" if len(input.fastq) == 2 else f"-U {input.fastq[0]}"
    shell:
        """
        mkdir -p $(dirname {output.outfile})
        {params.HISAT2} -x {params.index_prefix} \
            {params.input_params} \
            -p {threads} 2> {log} | \
        {params.SAMTOOLS} sort -@ {threads} -@ {threads} -o {output.outfile}
        """
