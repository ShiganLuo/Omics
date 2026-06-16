from snakemake.logging import logger
import re
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"


rule SplitNCigarReads:
    input:
        bam = indir + "/bam-sorted-Markdup/{sample_id}.bam",
    output:
        bam = outdir + "/Split/{sample_id}/{sample_id}.bam"
    params:
        javaOptions = "-Xms20g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config["Procedure"]["gatk"],
        genome = lambda wildcards: config['genome']['fasta'],
        indict = lambda wildcards: config['genome']['dict_index'],
        fai = lambda wildcards: config['genome']['fai_index']
    threads: 8 
    log:
        logdir + "/{sample_id}/SplitNCigarReads.log"
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" SplitNCigarReads \
        -R {params.genome} \
        -I {input.bam} \
        -O {output.bam} > {log} 2>&1
      """

rule VarientCalling:
    input:
        bam = outdir + "/Split/{sample_id}/{sample_id}.bam"
    output:
        vcf = outdir + "/vcf/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/VarientCalling.log"
    params:
        javaOptions = "-Xms20g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config["Procedure"]["gatk"],
        genome = lambda wildcards: config['genome']['fasta'],
        fai = lambda wildcards: config['genome']['fai_index']
    threads: 8
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" \
		HaplotypeCaller \
		-R {params.genome} \
		-I {input.bam} \
		-O {output.vcf} \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20 > {log} 2>&1
        """

rule vcf_filter:
    input:
        vcf = outdir + "/vcf/{sample_id}/{sample_id}.vcf.gz"
    output:
        vcf = outdir + "/vcf-filtered/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/vcf_filter.log"
    threads: 8
    params:
        javaOptions = "-Xms20g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        vcf = outdir + "/SNP/vcf/filter/{genome}/{sample_id}.vcf",
        gatk = config["Procedure"]["gatk"],
        bgzip = config["Procedure"]["bgzip"],
        genome = lambda wildcards: config['genome']['fasta'],
        fai = lambda wildcards: config['genome']['fai_index']
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" VariantFiltration \
        --R {params.genome} \
        --V {input.vcf} \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        -O {params.vcf} > {log} 2>&1 
        {params.bgzip} {params.vcf} >> {log} 2>&1 
        """
