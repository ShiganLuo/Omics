rule TEtranscript_prepare:
    input:
        get_alignment_input,
        genome_index = config['STAR']['human']['genome_index']
    output:
        outfile = outdir + "/counts/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        log=outdir+"/log/human/{sample_id}/TEtranscript_prepare.log"
    threads: 15
    params:
        outPrefix = outdir + "/counts/{sample_id}/human/{sample_id}",
        STAR = config["STAR"]["procedure"],
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
        """
rule TEcount:
    input:
        bam = outdir + "/counts/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/counts/{sample_id}/human/{sample_id}TEcount.cntTable"
    params:
        project = "{sample_id}TEcount",
        outdir = outdir + "/counts/{sample_id}/human",
        TE_gtf = config['TEtranscripts']['human']['TE_gtf'],
        gtf = config['TEtranscripts']['human']['gtf']
    log:
        log=outdir+"/log/human/{sample_id}/TEtranscripts.log"
    conda:
        config['conda']['TE']
    shell:
        """
        TEcount --sortByPos --format BAM --mode multi \
        -b {input.bam} --GTF {params.gtf} --TE {params.TE_gtf} \
        --project {params.project} --outdir {params.outdir} \
        > {log.log} 2>&1
        """

rule combine_TEcount:
    input:
        fileList = expand(outdir + "/counts/{sample_id}/human/{sample_id}TEcount.cntTable",sample_id=samples)
    output:
        outfile = outdir + "/counts/humanTEcount.cntTable"
    conda:
        config['conda']['RNA-SNP']
    params:
        combineTE = "scripts/combineTE.py",
        indir = outdir + "/counts"
    log:
        log = outdir + "/log/human/combine_TEcount.log"
    shell:
        """
        python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile} > {log.log} 2>&1
        """
rule TElocal:
    input:
        bam = outdir + "/counts/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/counts/{sample_id}/human/{sample_id}TElocal.cntTable"
    params:
        project = "{sample_id}TElocal",
        TE = config['TElocal']['human']['TEind'],
        GTF = config['TElocal']['human']['gtf'],
        procedure = "/opt/TElocal/TElocal"
    log:
        log = outdir+"/log/human/{sample_id}/TElocal.log"
    conda:
        config['conda']['TElocal']
    shell:
        """
        which python
        {params.procedure} --sortByPos -b {input.bam} \
        --GTF {params.GTF} --TE {params.TE} \
        --project {params.project} > {log.log} 2>&1
        mv {params.project}.cntTable {output.project}
        """

rule combine_TElocal:
    input:
        fileList = expand(outdir + "/counts/{sample_id}/human/{sample_id}TElocal.cntTable",sample_id=samples)
    output:
        outfile = outdir + "/counts/humanTElocal.cntTable"
    conda:
        config['conda']['RNA-SNP']
    params:
        combineTE = "scripts/combineTE.py",
        indir = outdir + "/counts"
    log:
        log = outdir + "/log/human/combine_TElocal.log"
    shell:
        """
        python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile} > {log.log} 2>&1
        """