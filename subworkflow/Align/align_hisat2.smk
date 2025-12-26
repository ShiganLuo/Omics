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

rule hisat2_align:
    input:
        fastq = get_alignment_input,
        index = lambda wildcards: expand(
            outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}.{{idx}}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
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
        index_prefix = lambda wildcards, input: input.index[0].replace(".1.ht2", ""),
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
