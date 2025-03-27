shell.prefix("set -x; set -e;")
configfile: "config/RNASNP.yaml"
indir = config.get('indir', '../data')  # 如果没有传递 indir，则使用 'data' 作为默认值
outdir = config.get('outdir', '../output')  # 如果没有传递 outdir，则使用 'output' 作为默认值
samples = glob_wildcards(indir + "/fq/{sample_id}_1.fq.gz").sample_id
# genomes = ['human','mouse']
genomes = ['human','mouse']
print("数据输入路径:", indir,"\t数据输出路径:", outdir)   
print(samples)
print(config['genome'])
def request():
    ##对于能通过依赖关系寻找的中间文件不需要重复定义，否则执行次数会过多
    output = dict()
    ### star_align
    # output['star_align'] =  expand(outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam", sample_id=samples, genome=genomes)
    ### xenofilter
    output['xenofilter'] = [outdir + "/xenofilterR/bam"]
    output['VarientCalling'] = [expand(outdir+"/Split/vcf/{sample_id}.vcf.gz",sample_id=samples)]
    output['filter'] = [expand(outdir+"/filter/vcf/{sample_id}.vcf.gz",sample_id=samples)]
    # output['gatk_index'] = []
    return list(output.values())
# print(request())
rule all:
    input:
        request()
        
rule trimming:
    input:
        fastq1=indir+"/fq/{sample_id}_1.fq.gz",
        fastq2=indir+"/fq/{sample_id}_2.fq.gz"
    output:
        fastq1=outdir+"/cutadapt/{sample_id}_1.fq.gz",
        fastq2=outdir+"/cutadapt/{sample_id}_2.fq.gz",
        report1=outdir+"/log/human/{sample_id}/trimming_statistics_1.txt",
        report2=outdir+"/log/human/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir=outdir+"/cutadapt",
        quality=30,
        trim_galore="/opt/TrimGalore-0.6.10/trim_galore"
    threads: 8
    log:
        log=outdir+"/log/human/{sample_id}/trimming.txt"
    shell:
        """
        {params.trim_galore} --phred33 --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fq.gz_trimming_report.txt {output.report2}
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
# 修改后的比对规则（支持多基因组）
rule star_align:
    input:
        fastq1 = outdir + "/cutadapt/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/cutadapt/{sample_id}_2.fq.gz",
        genome_index = lambda wildcards: config['STAR'][wildcards.genome]['genome_index'] ,
    output:
        outfile = outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        log=outdir+"/log/human/{sample_id}/{genome}/star_align.log"
    threads:25
    params:
        outPrefix = outdir + "/2pass/{sample_id}/{genome}/{sample_id}",
        STAR = config["STAR"]["procedure"]
    shell:
        """
        mkdir -p $(dirname {params.outPrefix})
        {params.STAR} --runThreadN {threads} \
            --genomeDir {input.genome_index} \
            --twopassMode Basic \
            --readFilesCommand zcat \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NM \
            --outFileNamePrefix {params.outPrefix} > {log.log} 2>&1
        """
# 生成Xenofilter输入CSV,第一列为比对到需要的基因组bam，第二列比对到怀疑有小鼠基因组污染的bam
rule generate_xenofilter_input:
    input:
        expand(outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam", 
               sample_id=samples, genome=genomes)
    output:
        csvIn = outdir + "/xenofilterR/xenofilterR_input.csv",
        csvRe = outdir +"/xenofilterR/xenofilterR_reName.csv"
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
        with open(output.csvRe, 'w', newline='') as f:
            writer = csv.writer(f)
            for sample in samples:
                row = [f"{sample}_xenofilterR"]
                writer.writerow(row)


# XenofilterR处理规则
rule XenofilterR:
    input:
        csvIn = outdir + "/xenofilterR/xenofilterR_input.csv",
        csvRe = outdir +"/xenofilterR/xenofilterR_reName.csv"
    output:
        outdir = directory(outdir + "/xenofilterR/bam")
    log:
        log = outdir + "/log/human/XenofilterR.log"
    threads: 3
    params:
        script = "scripts/XenofilteR.r",
        threshold=8
    shell:
        """
        #-p不仅递归创建还不会警告
        mkdir -p {output.outdir}
        /usr/bin/Rscript {params.script} \
            --inputFile {input.csvIn} \
            --renameFile {input.csvRe} \
            --outputDir {output.outdir} \
            --MM {params.threshold} \
            --workers {threads} > {log.log} 2>&1
        """
rule  addReadsGroup:
    input:
        outdir = outdir + "/xenofilterR/bam",
    output:
        bam = temp(outdir + "/RG/{sample_id}.bam"),
        bai = temp(outdir + "/RG/{sample_id}.bam.bai")
    log:
        log = outdir + "/log/human/{sample_id}/addReadsGroup.log"
    threads:16
    conda:
        config['conda']['cfDNA_base']
    params:
        bam = outdir + "/xenofilterR/bam/Filtered_bams/{sample_id}_xenofilterR_Filtered.bam",#放入params躲避检查，隐式依赖
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
        bam=outdir+"/RG/{sample_id}.bam"
    output:
        bam=temp(outdir+"/bam-sorted-Markdup/{sample_id}.bam"),
        bai=temp(outdir+"/bam-sorted-Markdup/{sample_id}.bai"),
        metrics=temp(outdir+"/log/{sample_id}/Markdup-metrics.txt")
    log:
        log=outdir+"/bam-sorted-Markdup/log/{sample_id}/MarkDuplicates.log"
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
        genome=config['genome']['human'],
        bam=outdir+"/bam-sorted-Markdup/{sample_id}.bam",
        indict=config['genome_dict']['human']
    output:
        bam = outdir+"/Split/bam/{sample_id}.bam"
    params:
        javaOptions="-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
    log:
        log=outdir+"/log/human/{sample_id}/SplitNCigarReads.log"
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
        genome=config['genome']['human'],
        bam=outdir+"/Split/bam/{sample_id}.bam"
    output:
        vcf=outdir+"/Split/vcf/{sample_id}.vcf.gz"
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
        genome=config['genome'],
        vcf=outdir+"/Split/vcf/{sample_id}.vcf.gz",
    output:
        vcf=outdir+"/filter/vcf/{sample_id}.vcf.gz",
    log:
        log=outdir+"/log/human/{sample_id}/vcf_filter.log"
    conda:
        config['conda']['cfDNA_base']
    params:
        javaOptions="-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
        vcf = outdir+"/filter/vcf/{sample_id}.vcf",
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