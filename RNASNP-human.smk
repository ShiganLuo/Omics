shell.prefix("set -x; set -e;")
configfile: "config/RNASNP.yaml"
indir = config.get('indir', '../data')
outdir = config.get('outdir', '../output')
# 动态获取所有样本（包含单端和双端）
paired_samples = glob_wildcards(indir + "/fq/{sample_id}_1.fq.gz").sample_id
Allsamples = glob_wildcards(indir + "/fq/{sample_id}.fq.gz").sample_id
single_samples = [sample for sample in Allsamples if re.match(r"^SRR\d+$", sample)]
samples = paired_samples + single_samples
print(samples)
# genomes
genomes1=['mouse','human']
genomes=['human']

def request():
    ##对于能通过依赖关系寻找的中间文件不需要重复定义，否则执行次数会过多
    output = dict()
    ### star_align
    output['star_align'] =  [expand(outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam", sample_id=samples, genome=genomes1)]
    output['xenoInput'] = [
        outdir + "/xenofilterR/xenofilterR_input.csv"
    ]
    output['VarientCalling'] = [expand(outdir + "/Split/vcf/human/{sample_id}.vcf.gz",sample_id=samples)]
    output['filter'] = [expand(outdir + "/filter/vcf/human/{sample_id}.vcf.gz",sample_id=samples)]
    output['count_Gene'] = [expand(outdir + "/counts/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam",sample_id=samples)]
    output['count_TE'] = [expand(outdir + "/counts/{sample_id}/human/{sample_id}TElocal.cntTable",sample_id=samples)]
    output['combineTEcount'] = [expand(outdir + "/counts/humanTEcount.cntTable")]
    output['combineTElocal'] = [expand(outdir + "/counts/humanTElocal.cntTable")]
    output['stringTieOverlap'] = [expand(outdir + "/2pass/{genome}_StgTEOverlap.bed",genome=genomes)]
    output['annovar'] = expand(outdir + "/annovar/human/{sample_id}/{sample_id}.GRCh38_multianno.csv",sample_id=samples)
    # output['gatk_index'] = []
    return list(output.values())
# print(request())
rule all:
    input:
        request()
        
rule trimming_Paired:
    input:
        fastq1=indir+"/fq/{sample_id}_1.fq.gz",
        fastq2=indir+"/fq/{sample_id}_2.fq.gz"
    output:
        fastq1=outdir+"/cutadapt/{sample_id}_1.fq.gz",
        fastq2=outdir+"/cutadapt/{sample_id}_2.fq.gz",
        report1=outdir+"/log/{sample_id}/trimming_statistics_1.txt",
        report2=outdir+"/log/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir=outdir+"/cutadapt",
        quality=30,
        trim_galore="/opt/TrimGalore-0.6.10/trim_galore"
    threads: 6
    log:
        log=outdir+"/log/{sample_id}/trimming.txt"
    shell:
        """
        # trim_galore can automatically judge the fq quality scoring system,it's no need to add such as --phred33 --phred64
        {params.trim_galore} --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fq.gz_trimming_report.txt {output.report2}
        """
rule trimming_Single:
    input:
        fastq = indir+"/fq/{sample_id}.fq.gz"
    output:
        fastq = outdir+"/cutadapt/{sample_id}Single.fq.gz",
        report = outdir+"/log/{sample_id}/trimming_statistics.txt"
    params:
        outdir=outdir+"/cutadapt",
        quality=30,
        trim_galore="/opt/TrimGalore-0.6.10/trim_galore"
    threads: 6
    log:
        log=outdir+"/log/{sample_id}/trimming.txt"
    shell:
        """
        {params.trim_galore} --phred33  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_trimmed.fq.gz {output.fastq}
        mv {params.outdir}/{wildcards.sample_id}.fq.gz_trimming_report.txt {output.report}
        """

def get_alignment_input(wildcards):
    """动态判断输入文件类型"""
    # 构造可能的输入路径
    paired_r1 = f"{outdir}/cutadapt/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{outdir}/cutadapt/{wildcards.sample_id}_2.fq.gz"
    single = f"{outdir}/cutadapt/{wildcards.sample_id}Single.fq.gz"
    
    # 检查文件实际存在情况
    if wildcards.sample_id in paired_samples:
        print(f"双端测序：{[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        print(f"单端测序：{[single]}")
        return [single]
    else:
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample_id}\n"
            f"Checked paths:\n- {paired_r1}\n- {paired_r2}\n- {single}"
        )

rule star_align:
    input:
        fastq = get_alignment_input,
        genome_index = lambda wildcards: config['STAR'][wildcards.genome]['genome_index']
    output:
        outfile = outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        log = outdir + "/log/{sample_id}/{genome}/star_align.log"
    threads: 25
    params:
        outPrefix = outdir + "/2pass/{sample_id}/{genome}/{sample_id}",
        STAR = config["STAR"]["procedure"],
        # 动态判断输入参数,加上genome_index，如果三个参数，即为双端测序，两个参数即为单端测序
        input_params = lambda wildcards, input: \
            f"{input[0]} {input[1]}" if len(input) == 3 else f"{input[0]}"
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
            --outFileNamePrefix {params.outPrefix} > {log.log} 2>&1
        """

# rule star_index:
#     input:
#         genome_fa = lambda wildcards: config['STAR'][wildcards.genome]['genome_gtf'] ,
#         genome_gtf = lambda wildcards: config['STAR'][wildcards.genome]['genome_gtf'] 
#     output:
#         genome_index = directory("/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/star")
#     log:
#         log=outdir+"/log/{genome}/star_index.log"
#     params:
#         STAR = config["STAR"]["procedure"],
#         read_length = config["STAR"]["read_length"],
#         index_core = config["STAR"]["index_core"]
#     shell:
#         """
#         {params.STAR} --runMode genomeGenerate \
#             --runThreadN {params.index_core} \
#             --genomeDir {output.genome_index} \
#             --genomeFastaFiles {input.genome_fa} \
#             --sjdbGTFfile {input.genome_gtf} \
#             --sjdbOverhang {params.read_length} > {log.log} 2>&1
        # """

# 生成Xenofilter输入CSV,第一列为比对到需要的基因组bam，第二列比对到怀疑有小鼠基因组污染的bam
rule generate_xenofilter_input:
    input:
        expand(outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam", 
               sample_id=samples, genome=genomes)
    output:
        csvIn = outdir + "/xenofilterR/xenofilterR_input.csv",
        # csvRe = outdir +"/xenofilterR/xenofilterR_reName.csv"
    run:
        import csv
        with open(output.csvIn, 'w', newline='') as f:
            writer = csv.writer(f)
            for sample in samples:
                row = [
                    f"{outdir}/2pass/{sample}/human/{sample}Aligned.sortedByCoord.out.bam",
                    f"{outdir}/2pass/{sample}/mouse/{sample}Aligned.sortedByCoord.out.bam"
                ]
                writer.writerow(row)
        # 选择重命名文件
        # with open(output.csvRe, 'w', newline='') as f:
        #     writer = csv.writer(f)
        #     for sample in samples:
        #         row = [f"{sample}_xenofilterR"]
        #         writer.writerow(row)


# XenofilterR处理规则
rule XenofilterR:
    input:
        csvIn = outdir + "/xenofilterR/xenofilterR_input.csv",
    output:
        # expand(outdir + "/xenofilterR/Filtered_bams/{sample_id}_Filtered.bam",sample_id=samples),
        # expand(outdir + "/xenofilterR/Filtered_bams/{sample_id}_Filtered.bam.bai",sample_id=samples)
        outdir = directory(outdir + "/xenofilterR/bam") #XenofilteR设计不合理，没办法
    log:
        log = outdir + "/log/human/XenofilterR.log"
    threads: 6
    params:
        script = "scripts/XenofilteR.r",
        threshold=8
    shell:
        """
        /usr/bin/Rscript {params.script} \
            --inputFile {input.csvIn} \
            --outputDir {output.outdir} \
            --MM {params.threshold} \
            --workers {threads} > {log.log} 2>&1
        """
#####################SNP########################
rule  addReadsGroup:
    input:
        outdir = outdir + "/xenofilterR/bam"
    output:
        bam = temp(outdir + "/RG/human/{sample_id}.bam"),
        bai = temp(outdir + "/RG/human/{sample_id}.bam.bai")
    log:
        log = outdir + "/log/human/{sample_id}/addReadsGroup.log"
    threads:16
    conda:
        config['conda']['cfDNA_base']
    params:
        bam = outdir + "/xenofilterR/bam/Filtered_bams/{sample_id}Aligned.sortedByCoord.out_Filtered.bam",
        id="{sample_id}",
        java="--java-options -Xmx15G",
        RGLB=config["addReadsGroup"]["RGLB"],
        RGPL=config["addReadsGroup"]["RGPL"],
        RGPU=config["addReadsGroup"]["RGPU"]
    shell:
        """
        echo "sample_id: {wildcards.sample_id}" >> {log.log}
        gatk AddOrReplaceReadGroups {params.java} \
            --INPUT {params.bam} --OUTPUT {output.bam} \
            -SO coordinate --RGLB {params.RGLB} --RGPL {params.RGPL} --RGPU {params.RGPU} --RGSM {params.id} > {log.log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log.log} 2>&1
        """
rule MarkDuplicates:
    input:
        bam = outdir + "/RG/human/{sample_id}.bam"
    output:
        bam=temp(outdir + "/bam-sorted-Markdup/human/{sample_id}.bam"),
        bai=temp(outdir + "/bam-sorted-Markdup/human/{sample_id}.bai"),
        metrics=temp(outdir + "/log/human/{sample_id}/Markdup-metrics.txt")
    log:
        log=outdir+"/bam-sorted-Markdup/log/human/{sample_id}/MarkDuplicates.log"
    threads: 16
    conda:
        config['conda']['cfDNA_base']
    params:
        javaOptions="-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
    shell:
        """
        gatk --java-options "{params.javaOptions}" MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam}   \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --METRICS_FILE {output.metrics} > {log.log} 2>&1
        """
# rule gatk_index:
#     input:
#         genome=config['genome']['human']
#     output:
#         outdict="/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38.primary_assembly.genome.dict"
#     log:
#         log=outdir+"/log/gatk_index.log"
#     conda:
#         config['conda']['cfDNA_base']
#     params:
#     shell:
#         """
#         gatk CreateSequenceDictionary -R {input.genome} -O {output.outdict} > {log.log} 2>&1
#         samtools faidx {input.genome} >> {log.log} 2>&1
#         """

rule SplitNCigarReads:
    input:
        genome = config['genome']['human'],
        bam = outdir + "/bam-sorted-Markdup/human/{sample_id}.bam",
        indict = config['genome_dict']['human']
    output:
        bam = outdir + "/Split/bam/human/{sample_id}.bam"
    params:
        javaOptions="-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
    log:
        log = outdir + "/log/human/{sample_id}/SplitNCigarReads.log"
    conda:
        config['conda']['cfDNA_base']
    shell:
        """
        gatk --java-options "{params.javaOptions}" SplitNCigarReads \
        -R {input.genome} \
        -I {input.bam} \
        -O {output.bam} > {log.log} 2>&1
      """
rule VarientCalling:
    input:
        genome = config['genome']['human'],
        bam = outdir + "/Split/bam/human/{sample_id}.bam"
    output:
        vcf = outdir + "/Split/vcf/human/{sample_id}.vcf.gz"
    conda:
        config['conda']['cfDNA_base']
    log:
        log=outdir+"/log/human/{sample_id}/VarientCalling.log"
    params:
        javaOptions="-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
    shell:
        """
        gatk --java-options "{params.javaOptions}" \
		HaplotypeCaller \
		-R {input.genome} \
		-I {input.bam} \
		-O {output.vcf} \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20 > {log.log} 2>&1
        """
rule vcf_filter:
    input:
        genome = config['genome']['human'],
        vcf = outdir + "/Split/vcf/human/{sample_id}.vcf.gz"
    output:
        vcf = outdir + "/filter/vcf/human/{sample_id}.vcf.gz"
    log:
        log=outdir+"/log/human/{sample_id}/vcf_filter.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        javaOptions="-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        vcf = outdir+"/filter/vcf/human/{sample_id}.vcf",
    shell:
        """
        gatk --java-options "{params.javaOptions}" VariantFiltration \
        --R {input.genome} \
        --V {input.vcf} \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        -O {params.vcf} > {log.log} 2>&1 
        bgzip {params.vcf} >> {log.log} 2>&1 
        """
####################################TEtranscripts############################################
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
        exon_gtf = config['TEtranscripts']['human']['exon_gtf']
    log:
        log=outdir+"/log/human/{sample_id}/TEtranscripts.log"
    conda:
        config['conda']['TE']
    shell:
        """
        TEcount --sortByPos --format BAM --mode multi \
        -b {input.bam} --GTF {params.exon_gtf} --TE {params.TE_gtf} \
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
        GTF = config['TElocal']['human']['exon_gtf'],
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
###################转录组gtf#######################
rule stringTie:
    input:
        bam = outdir + "/2pass/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        gtf = outdir + "/2pass/{sample_id}/human/{sample_id}.gtf"
    conda:
        config['conda']['RNA-SNP']
    log:
        log = outdir + "/log/human/{sample_id}/stringTie.log"
    shell:
        """
        stringtie -o {output.gtf} {input.bam} > {log.log} 2>&1
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
    shell:
        """
        stringtie --merge {input.gtf} -o {output.outfile} > {log.log} 2>&1
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
            --exonGtf {input.gtf} \
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


#####################annovar########################
rule commonExpression:
    input:
        infile = outdir + "/counts/humanTEcount.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.cntTable"
    log:
        log = outdir + "/log/human/commonExpression.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = "scripts/SNP/commonExpression.py"
    shell:
        """
            python {params.script} --input {input.infile} --output {output.outfile} > {log.log} 2>&1
        """
rule getBed:
    input:
        infile = outdir + "/counts/humanTEcountCommon.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.bed"
    log:
        log = outdir + "/log/human/getBed.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = "scripts/SNP/getBed.py",
        exon_gtf = config['getBed']['human']['exon_gtf'],
        TE_gtf = config['getBed']['human']['TE_gtf']
    shell:
        """
            python {params.script} \
                --input {input.infile} \
                --output {output.outfile} \
                --exonGtf {params.exon_gtf} \
                --TEGtf {params.TE_gtf} > {log.log} 2>&1
        """
rule vcfIntersectBed:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}.vcf.gz",
        bed = outdir + "/counts/humanTEcountCommon.bed"
    output:
        outfile = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    log:
        log = outdir + "/log/human/{sample_id}/vcfIntersectBed.log"
    conda:
        config['conda']['RNA-SNP']
    threads:4 #防止同时执行太多，爆内存
    shell:
        """
            bedtools intersect -a {input.vcf} -b {input.bed} -wa -wb > {output.outfile} 2>{log.log}
        """
rule annovar_convert:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    output:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_convert.log"
    params:
        convert = "/opt/annovar/convert2annovar.pl"
    threads:4 #防止同时执行太多，爆内存
    shell:
        """
        /usr/bin/perl {params.convert} \
        -format vcf4 \
        -withfreq {input.vcf} > {output.avinput} 2>{log.log}
        """
rule annovar_table:
    input:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    output:
        outfile = outdir + "/annovar/human/{sample_id}/{sample_id}.GRCh38_multianno.csv"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_table.log"
    params:
        db = config['annovar']['human']['db'],
        buildver = config['annovar']['human']['buildver'],
        # annotate = "/opt/annovar/annotate_variation.pl",
        table = "/opt/annovar/table_annovar.pl",
        out = outdir + "/annovar/human/{sample_id}/{sample_id}"
    threads:4 #防止同时执行太多，爆内存
    shell:
        """
        /usr/bin/perl {params.table} \
        {input.avinput} {params.db} \
        -buildver {params.buildver} \
        -out {params.out} \
        -remove -protocol refGene \
        -operation g \
        -nastring . \
        -csvout > {log.log} 2>&1
        """