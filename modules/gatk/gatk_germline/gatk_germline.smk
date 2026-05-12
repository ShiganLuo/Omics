outdir = config.get("outdir") or "output"
fasta = config.get("fasta")
fai_index = config.get("fai_index")
dict_index = config.get("dict_index")
################################################################
# SNV (option1: GATK germline mode using Haplotypecaller)
################################################################
def get_input_HaplotypeCaller(wildcards):
    """
    This function determines the appropriate input BAM file for the HaplotypeCaller step based on the presence of known_sites and interval parameters.
    If both known_sites and interval are provided, it assumes that BQSR has been performed and uses the BQSR-corrected BAM file. Otherwise, it falls back to using the MarkDuplicates BAM file. 
    It also checks for the presence of the corresponding index files and includes them as input if they exist.
    """
    in_dcit = {}
    if known_sites and interval:
        print(f"Using known_sites: {known_sites} and interval: {interval}")
        include: "bqsr.smk"
        in_dcit["bam"] = f"{outdir}/{wildcards.genome}/gatk/bqsr/{wildcards.sample_id}.sorted.markdup.BQSR.bam"
    else:
        print("No known_sites or interval specified in config, proceeding without them.")
        in_dcit["bam"] = f"{outdir}/{wildcards.genome}/gatk/bam-sorted-Markdup/{wildcards.sample_id}.bam"
    
    if fai_index and dict_index:
        in_dcit["fai"] = fai_index
        in_dcit["dict"] = dict_index
    else:
        in_dcit["fai"] = f"{outdir}/{wildcards.genome}/gatk/index/{wildcards.genome}.fai"
        in_dcit["dict"] = f"{outdir}/{wildcards.genome}/gatk/index/{wildcards.genome}.dict"
    return in_dcit

rule HaplotypeCaller:
    input:
        bam = get_input_HaplotypeCaller
    output:
        vcf = outdir + "/gatk/vcf/{genome}/{sample_id}.vcf.gz"
    log:
        outdir + "/log/gatk/{genome}/{sample_id}/haplotypeCaller.log"
    conda:
        "../gatk.yaml"
    params:
        ref = fasta,
        javaOptions =  "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk"
    threads: 10
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" HaplotypeCaller \
            -R {params.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            > {log} 2>&1
        """

rule filterHaplotypeCallerVcf:
    input:
        vcf = outdir + "/gatk/vcf/{genome}/{sample_id}.vcf.gz"
    output:
        vcf = outdir+"/gatk/vcf-filtered/{genome}/{sample_id}.vcf.gz"
    log:
        outdir + "/log/gatk/{genome}/{sample_id}/haplotypeCaller-filtering.log"
    conda:
        "../gatk.yaml"
    threads: 10
    params:
        ref = fasta,
        javaOptions = "--java-options -Xmx35G", # "--java-options -Xmx15G",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk"
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" VariantFiltration \
            -R {params.ref} -V {input.vcf} -O {output.vcf} \
            -window 35 -cluster 3 \
            --filter-name FS20 -filter "FS > 20.0" \
            --filter-name QD2 -filter "QD < 2.0" \
            --filter-name DP10 -filter "DP < 10.0" \
            --filter-name QUAL20 -filter "QUAL < 20.0" \
            > {log} 2>&1
        """

rule gatk_haplotypeCaller_result:
    input:
        vcf = outdir + "/gatk/vcf-filtered/{genome}/{sample_id}.vcf.gz"
