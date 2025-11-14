SNAKEFILE_FULL_PATH_SNP = workflow.snakefile
SNAKEFILE_DIR_SNP = os.path.dirname(SNAKEFILE_FULL_PATH_SNP)
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_SNP ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
SNPYaml = get_yaml_path("SNP")
configfile: SNPYaml
logging.info(f"Include SNP config: {SNPYaml}")

rule addReadsGroup:
    input:
        outdir = outdir + "/xenofilterR/bam"
    output:
        bam = temp(outdir + "/SNP/RG/{genome}/{sample_id}.bam"),
        bai = temp(outdir + "/SNP/RG/{genome}/{sample_id}.bam.bai")
    log:
        outdir + "/log/SNP/{genome}/{sample_id}/addReadsGroup.log"
    threads:16
    params:
        bam = lambda wildcards: \
            f"{outdir}/xenofilterR/bam/Filtered_bams/{wildcards.sample_id}Aligned.sortedByCoord.out_Filtered.bam" \
            if wildcards.genome == "Homo_sapiens" else \
            f"{outdir}/2pass/{wildcards.sample_id}/{wildcards.genome}/{wildcards.sample_id}Aligned.sortedByCoord.out.bam",
        id = "{sample_id}",
        javaOptions = "--java-options -Xmx15G",
        RGLB = config["addReadsGroup"]["RGLB"],
        RGPL = config["addReadsGroup"]["RGPL"],
        RGPU = config["addReadsGroup"]["RGPU"],
        gatk = config["Procedure"]["gatk"],
        samtools = config["Procedure"]["samtools"]
    shell:
        """
        echo "sample_id: {wildcards.sample_id}" > {log}
        {params.gatk} AddOrReplaceReadGroups {params.javaOptions} \
            --INPUT {params.bam} --OUTPUT {output.bam} \
            -SO coordinate --RGLB {params.RGLB} --RGPL {params.RGPL} --RGPU {params.RGPU} --RGSM {params.id} >> {log} 2>&1
        {params.samtools} index -@ {threads} {output.bam} >> {log} 2>&1
        """

rule MarkDuplicates:
    input:
        bam = outdir + "/SNP/RG/{genome}/{sample_id}.bam"
    output:
        bam = temp(outdir + "/SNP/bam-sorted-Markdup/{genome}/{sample_id}.bam"),
        bai = temp(outdir + "/SNP/bam-sorted-Markdup/{genome}{sample_id}.bai"),
        metrics = temp(outdir + "/SNP/bam-sorted-Markdup/{genome}/{sample_id}_Markdup-metrics.txt")
    log:
        outdir + "/log/SNP/{genome}/{sample_id}/{sample_id}_MarkDuplicates.log"
    threads: 16
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config["Procedure"]["gatk"]
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam}   \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --METRICS_FILE {output.metrics} > {log} 2>&1
        """
# rule gatk_index:
#     input:
#         genome=config['genome']['human']
#     output:
#         outdict="/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38.primary_assembly.genome.dict"
#     log:
#         outdir+"/log/gatk_index.log"
#     params:
#         gatk = config["Procedure"]["gatk"],
#         samtools = config["Procedure"]["samtools"]
#     shell:
#         """
#         {params.gatk} CreateSequenceDictionary -R {input.genome} -O {output.outdict} > {log} 2>&1
#         {params.samtools} faidx {input.genome} >> {log} 2>&1
#         """

rule SplitNCigarReads:
    input:
        genome = config['genome']['human']['fasta'],
        bam = outdir + "/SNP/bam-sorted-Markdup/{genome}/{sample_id}.bam",
        indict = config['genome']['human']['dict'],
        fai = config['genome']['human']['fai']
    output:
        bam = outdir + "/SNP/Split/{genome}/{sample_id}.bam"
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config["Procedure"]["gatk"]
    log:
        outdir + "/log/SNP/{genome}/{sample_id}/{sample_id}_SplitNCigarReads.log"
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" SplitNCigarReads \
        -R {input.genome} \
        -I {input.bam} \
        -O {output.bam} > {log} 2>&1
      """

rule VarientCalling:
    input:
        genome = config['genome']['human']['fasta'],
        bam = outdir + "/SNP/Split/{genome}/{sample_id}.bam",
        fai = config['genome']['human']['fai']
    output:
        vcf = outdir + "/SNP/vcf/origin/{genome}/{sample_id}.vcf.gz"
    log:
        outdir + "/log/SNP/{genome}/{sample_id}/{sample_id}_VarientCalling.log"
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config["Procedure"]["gatk"]
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" \
		HaplotypeCaller \
		-R {input.genome} \
		-I {input.bam} \
		-O {output.vcf} \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20 > {log} 2>&1
        """

rule vcf_filter:
    input:
        genome = config['genome']['human']['fasta'],
        vcf = outdir + "/SNP/vcf/origin/{genome}/{sample_id}.vcf.gz",
        fai = config['genome']['human']['fai']
    output:
        vcf = outdir + "/SNP/vcf/filter/{genome}/{sample_id}.vcf.gz"
    log:
        outdir + "/log/SNP/{genome}/{sample_id}/{sample_id}_vcf_filter.log"
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        vcf = outdir + "/SNP/vcf/filter/{genome}/{sample_id}.vcf",
        gatk = config["Procedure"]["gatk"],
        bgzip = config["Procedure"]["bgzip"],
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" VariantFiltration \
        --R {input.genome} \
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
