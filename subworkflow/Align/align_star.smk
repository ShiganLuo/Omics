rule star_align:
    input:
        fastq = get_alignment_input,
        genome_index = lambda wildcards: config['genome'][wildcards.genome]['genome_index']
    output:
        outfile = outdir + "/Align/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        outdir + "/log/Align/{sample_id}/{genome}/star_align.log"
    threads: 12
    params:
        outPrefix = outdir + "/Align/{sample_id}/{genome}/{sample_id}",
        # 动态判断输入参数,加上genome_index，如果三个参数，即为双端测序，两个参数即为单端测序
        input_params = lambda wildcards, input: \
            f"{input[0]} {input[1]}" if len(input) == 3 else f"{input[0]}",
        STAR = config['Procedure']['STAR']
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