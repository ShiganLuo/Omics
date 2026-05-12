outdir = config.get("outdir") or "output"
fasta = config.get("fasta")
fai_index = config.get("genome", {}).get("fai_index")
dict_index = config.get("genome", {}).get("dict_index")
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
        in_dcit["bam"] = f"{indir}/bqsr/{wildcards.sample_id}.sorted.markdup.BQSR.bam"
    else:
        print("No known_sites or interval specified in config, proceeding without them.")
        in_dcit["bam"] = f"{indir}/bam-sorted-Markdup/{wildcards.sample_id}.bam"

    if fai_index and dict_index:
        in_dcit["fai"] = fai_index
        in_dcit["dict"] = dict_index
    else:
        in_dcit["fai"] = f"{outdir}/index/genome.fa.fai"
        in_dcit["dict"] = f"{outdir}/index/genome.dict"
    return in_dcit

rule HaplotypeCaller:
    input:
        bam = get_input_HaplotypeCaller
    output:
        vcf = outdir + "/vcf/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/haplotypeCaller.log"
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
        vcf = outdir + "/vcf/{sample_id}/{sample_id}.vcf.gz"
    output:
        vcf = outdir+"/vcf-filtered/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/haplotypeCaller-filtered.log"
    conda:
        "../gatk.yaml"
    threads: 10
    params:
        ref = fasta,
        javaOptions = "--java-options -Xmx35G", # "--java-options -Xmx15G",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        FS_threshold = config.get("Params", {}).get("gatk", {}).get("FS_threshold") or 20.0,
        QD_threshold = config.get("Params", {}).get("gatk", {}).get("QD_threshold") or 2.0,
        DP_threshold = config.get("Params", {}).get("gatk", {}).get("DP_threshold") or 10.0,
        QUAL_threshold = config.get("Params", {}).get("gatk", {}).get("QUAL_threshold") or 20.0
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" VariantFiltration \
            -R {params.ref} -V {input.vcf} -O {output.vcf} \
            -window 35 -cluster 3 \
            --filter-name FS -filter "FS > {params.FS_threshold}" \
            --filter-name QD -filter "QD < {params.QD_threshold}" \
            --filter-name DP -filter "DP < {params.DP_threshold}" \
            --filter-name QUAL -filter "QUAL < {params.QUAL_threshold}" \
            > {log} 2>&1
        """

rule gatk_haplotypeCaller_result:
    input:
        vcf = outdir + "/vcf-filtered/{sample_id}/{sample_id}.vcf.gz"
