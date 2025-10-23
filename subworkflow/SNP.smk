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
