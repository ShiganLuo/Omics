BW2_IDX_SUFFIX = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

rule bowtie2_index_small:
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta']
    output:
        index = expand(
            outdir + "/genome/{{genome}}/index/{{genome}}.{ext}",
            ext = BW2_IDX_SUFFIX
        )
    log:
        outdir + "/log/Align/bowtie2/{genome}/bowtie2_index.log"
    params:
        bowtie2_build = config.get('tools', {}).get('bowtie2-build') or 'bowtie2-build',
        prefix = outdir + "/genome/{genome}/index/{genome}"
    conda:
        config['conda']['run']
    threads: 8
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        {params.bowtie2_build} --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """

rule bowtie2_align:
    input:
        fastqs = get_alignment_input,
        index = expand(
            outdir + "/genome/{{genome}}/index/{{genome}}.{ext}",
            ext = BW2_IDX_SUFFIX
        )
    output:
        bam = outdir + "/Align/bam/{genome}/{sample_id}.bam"
    log:
        outdir + "/log/Align/bowtie2/{genome}/{sample_id}/bowtie2_align.log"
    params:
        bowtie2 = config.get('tools', {}).get('bowtie2') or 'bowtie2',
        samtools = config.get('tools', {}).get('samtools') or 'samtools',
        prefix = lambda wildcards, input: input.index[0].replace(".1.bt2", ""),
        args = lambda wildcards, input: \
            f"-1 {input.fastqs[0]} -2 {input.fastqs[1]}" if len(input.fastqs) == 2 else f"-U {input.fastqs[0]}"
    conda:
        config['conda']['run']
    threads: 10
    shell:
        """
        mkdir -p $(dirname {output.bam})
        {params.bowtie2} -x {params.prefix} {params.args} -N 1 -L 30 --threads {threads} 2> {log} | \
        {params.samtools} sort -@ {threads} -o {output.bam} - >> {log} 2>&1
        """