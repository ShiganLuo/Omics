rule stringTie:
    input:
        bam = outdir + "/2pass/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        gtf = outdir + "/2pass/{sample_id}/human/{sample_id}.gtf"
    conda:
        config['conda']['RNA-SNP']
    log:
        log = outdir + "/log/human/{sample_id}/stringTie.log"
    params:
        gtf = config['stringTie']['human']['gtf'] #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
    threads: 5
    shell:
        """
        stringtie -o {output.gtf} {input.bam} -G {params.gtf} -p {threads} > {log.log} 2>&1
        """
rule stringTieMerge:
    input:
        gtf = expand(outdir + "/2pass/{sample_id}/human/{sample_id}.gtf",sample_id=samples)
    output:
        outfile = outdir + "/2pass/human.gtf"
    conda:
        config['conda']['RNA-SNP']
    log:
        log = outdir + "/log/human/stringTieMerge.log"
    params:
        gtf = config['stringTie']['human']['gtf'] #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
    shell:
        """
        stringtie --merge {input.gtf} -o {output.outfile} -G {params.gtf} > {log.log} 2>&1
        """

rule TEcountStringTie:
    input:
        bam = outdir + "/counts/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam",
        gtf = outdir + "/2pass/human.gtf"
    output:
        project = outdir + "/counts/{sample_id}/human/{sample_id}TEcountStringTie.cntTable"
    params:
        project = "{sample_id}TEcountStringTie",
        outdir = outdir + "/counts/{sample_id}/human",
        TE_gtf = config['TEtranscripts']['human']['TE_gtf'],
    threads:5 #防止过多并行运行爆内存
    log:
        log=outdir+"/log/human/{sample_id}/TEtranscriptsStringTie.log"
    conda:
        config['conda']['TE']
    shell:
        """
        TEcount --sortByPos --format BAM --mode multi \
        -b {input.bam} --GTF {input.gtf} --TE {params.TE_gtf} \
        --project {params.project} --outdir {params.outdir} \
        > {log.log} 2>&1
        """

rule combine_TEStringtie:
    input:
        fileList = expand(outdir + "/counts/{sample_id}/human/{sample_id}TEcountStringTie.cntTable",sample_id=samples)
    output:
        outfile = outdir + "/counts/humanTEcountStringTie.cntTable"
    conda:
        config['conda']['RNA-SNP']
    params:
        combineTE = "scripts/combineTE.py",
        indir = outdir + "/counts"
    log:
        log = outdir + "/log/human/combine_TEcountStringTie.log"
    shell:
        """
        python {params.combineTE} -p TEcountStringTie -i {params.indir} -o {output.outfile} > {log.log} 2>&1
        """

rule getStringtieBed:
    input:
        gtf = outdir + "/2pass/human.gtf",
        infile = outdir + "/counts/humanTEcountStringTie.cntTable"
    output:
        genefile = outdir + "/2pass/human_STG.bed",
        TEfile = outdir + "/2pass/human_TE.bed"
    log:
        log = outdir + "/log/human/getStringtieBed.log"
    conda:
        config['conda']['RNA-SNP']
    threads:2
    params:
        script = "scripts/SNP/getBed.py",
        TE_gtf = config['TEtranscripts']['human']['TE_gtf']
    shell:
        """
        python {params.script} \
            --mode StringTie \
            --input {input.infile} \
            --output {output.genefile} \
            --output {output.TEfile} \
            --Gtf {input.gtf} \
            --TEGtf {params.TE_gtf} > {log.log} 2>&1
        """
rule StgTEOverlap:
    input:
        genefile = outdir + "/2pass/human_STG.bed",
        TEfile = outdir + "/2pass/human_TE.bed"
    output:
        outfile = outdir + "/2pass/human_StgTEOverlap.bed"
    log:
        log = outdir + "/log/human_StgTEOverlap.log"
    conda:
        config['conda']['RNA-SNP']
    threads:2
    shell:
        """
        bedtools intersect -a {input.genefile} -b {input.TEfile} -wa -wb > {output.outfile} 2>{log.log}
        """
