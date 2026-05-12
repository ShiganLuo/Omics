outdir = config.get("outdir") or "output"
fai_index = config.get("fai_index")
dict_index = config.get("dict_index")
fasta = config.get("fasta")
known_sites = config.get("known_sites")
interval = config.get("interval")
################################################################
# SNV (option2: GATK somatic mode using Mutect2)
################################################################
def get_input_somaticMutect2(wildcards):
    """
    This function determines the appropriate input BAM file for the HaplotypeCaller step based on the presence of known_sites and interval parameters.
    If both known_sites and interval are provided, it assumes that BQSR has been performed and uses the BQSR-corrected BAM file. Otherwise, it falls back to using the MarkDuplicates BAM file. 
    It also checks for the presence of the corresponding index files and includes them as input if they exist.
    """
    in_dict = {}

    if known_sites and interval:
        print(f"Using known_sites: {known_sites} and interval: {interval}")
        include: "bqsr.smk"
        in_dict["normal_bam"] = f"{outdir}/gatk/bqsr/{wildcards.genome}/{wildcards.normal_sample_id}.sorted.markdup.BQSR.bam"
        in_dict["experimental_bam"] = f"{outdir}/gatk/bqsr/{wildcards.genome}/{wildcards.experimental_sample_id}.sorted.markdup.BQSR.bam"
    else:
        print("No known_sites or interval specified in config, proceeding without them.")
        in_dict["normal_bam"] = f"{outdir}/gatk/bam-sorted-Markdup/{wildcards.genome}/{wildcards.normal_sample_id}.bam"
        in_dict["experimental_bam"] = f"{outdir}/gatk/bam-sorted-Markdup/{wildcards.genome}/{wildcards.experimental_sample_id}.bam"

    if fai_index and dict_index:
        in_dict["fai"] = fai_index
        in_dict["dict"] = dict_index
        in_dict["fasta"] = fasta
    else:
        in_dict["fai"] = f"{outdir}/gatk/index/{wildcards.genome}/{wildcards.genome}.fa.fai"
        in_dict["dict"] = f"{outdir}/gatk/index/{wildcards.genome}/{wildcards.genome}.dict"
        in_dict["fasta"] = f"{outdir}/gatk/index/{wildcards.genome}/{wildcards.genome}.fa"
    return in_dict


rule somaticMutect2:
    input:
       unpack(get_input_somaticMutect2)
    output:
        vcf = outdir + "/gatk/mutect2-vcf/{genome}/{normal_sample_id}_{experimental_sample_id}.vcf.gz"
    log:
        outdir + "/log/gatk/{genome}/{normal_sample_id}_{experimental_sample_id}/mutect2.log"
    conda:
        "../gatk.yaml"
    params:
        javaOptions =  "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        parameters = config.get("mutect2_parameters") or ""
    threads: 10
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" Mutect2 \
            -R {input.fasta} \
            -I {input.normal_bam} \
            -I {input.experimental_bam} \
            -normal {wildcards.normal_sample_id} \
            -O {output.vcf} \
            {params.parameters} \
            --native-pair-hmm-threads {threads} > {log} 2>&1
        """

