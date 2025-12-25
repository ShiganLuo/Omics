rule hisat2_index:
    input:
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output:
        index = lambda wildcards: expand(
            outdir + "/genome/{wildcards.genome}/index/hista2/GRCm39.primary_assembly.genome.{idx}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    threads: 8
    params:
        prefix = lambda wildcards: config["genome"][wildcards.genome]["hisat2_index_prefix"],
        HISAT2_BUILD = config["Procedure"]["hisat2-build"]
    log:
        lambda wildcards: f"results/log/genome/{wildcards.genome}/hisat2_build.log"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        {params.HISAT2_BUILD} -p {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """

rule hisat2_align:
    input:
        fastq = get_alignment_input,
        index = lambda wildcards: expand(
            outdir + "/genome/{wildcards.genome}/index/hista2/GRCm39.primary_assembly.genome.{idx}.ht2",
            idx = [1, 2, 3, 4, 5, 6, 7, 8]
        )
    output:
        outfile = outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        outdir + "/log/Align/{sample_id}/{genome}/hisat2_align.log"
    threads: 12
    params:
        HISAT2 = config['Procedure']['hisat2'],
        SAMTOOLS = config['Procedure']['samtools'],
        input_params = lambda wildcards, input: \
            f"-1 {input[0]} -2 {input[1]}" if len(input) == 3 else f"-U {input[0]}"
    shell:
        """
        mkdir -p $(dirname {output.outfile})
        {params.HISAT2} -x {input.index} \
            {params.input_params} \
            -p {threads} 2> {log} | \
        {params.SAMTOOLS} sort -@ {threads} -o {output.outfile}
        """

