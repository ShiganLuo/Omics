rule TEtranscript_prepare_hisat2:
    input:
        fastq = get_alignment_input,
        index = lambda wildcards: expand(
            outdir + f"/genome/{wildcards.genome}/index/hista2/{wildcards.genome}.{{idx}}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    output:
        outfile = temp(outdir + "/TEtranscripts/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
    log:
        outdir + "/log/TEtranscripts/{sample_id}/{genome}/hisat2_align.log"
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
            --no-mixed \
            --no-discordant \
            -k 100 \
            --score-min L,0,-0.2 \
            -p {threads} 2> {log} | \
        {params.SAMTOOLS} sort -@ {threads} -@ {threads} -o {output.outfile}
        """
