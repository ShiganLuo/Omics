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