def get_star_index(wildcards):
    logger.info(f"[get_star_index] called with wildcards: {wildcards}")
    star_index_dir = config.get('genome',{}).get(wildcards.genome,{}).get('star_index_dir') or None
    if star_index_dir:
        first_file = os.path.join(star_index_dir, "Genome")
        if os.path.exists(first_file):
            return star_index_dir

    return outdir + f"/genome/{wildcards.genome}/index/star"

rule TEtranscript_prepare_star:
    input:
        get_alignment_input,
        genome_index = get_star_index
    output:
        outfile = temp(outdir + "/TEtranscripts/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam")
    log:
        log = outdir + "/log/TEtranscripts/{genome}/{sample_id}/star_align.log",
        STAR_log = outdir + "/log/TEtranscripts/{genome}/{sample_id}/{sample_id}Log.out",
        STAR_progress = outdir + "/log/TEtranscripts/{genome}/{sample_id}/{sample_id}Log.progress.out",
        STAR_final = outdir + "/log/TEtranscripts/{genome}/{sample_id}/{sample_id}Log.final.out"
    threads: 12
    params:
        outPrefix = outdir + "/counts/bam/{genome}/{sample_id}",
        STAR = config['Procedure']['STAR'],
        # 动态判断输入参数,加上genome_index，如果三个参数，即为双端测序，两个参数即为单端测序
        input_params = lambda wildcards, input: \
            f"{input[0]} {input[1]}" if len(input) == 3 else f"{input[0]}"
    shell:
        """
        {params.STAR} --runThreadN {threads} \
            --genomeDir {input.genome_index} \
            --readFilesCommand zcat \
            --readFilesIn {params.input_params} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.outPrefix} \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 100  > {log.log} 2>&1
        cp {params.outPrefix}Log.out {log.STAR_log}
        cp {params.outPrefix}Log.progress.out {log.STAR_progress}
        cp {params.outPrefix}Log.final.out {log.STAR_final}
        """