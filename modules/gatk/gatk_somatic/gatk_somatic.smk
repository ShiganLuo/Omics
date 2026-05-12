from snakemake.logging import logger
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fai_index = config.get("genome", {}).get("fai_index")
dict_index = config.get("genome", {}).get("dict_index")
fasta = config.get("genome", {}).get("fasta")
known_sites = config.get("genome", {}).get("known_sites")
interval = config.get("genome", {}).get("interval")

def get_input_somaticMutect2(wildcards):
    """
    This function determines the appropriate input BAM file for the HaplotypeCaller step based on the presence of known_sites and interval parameters.
    If both known_sites and interval are provided, it assumes that BQSR has been performed and uses the BQSR-corrected BAM file. Otherwise, it falls back to using the MarkDuplicates BAM file. 
    It also checks for the presence of the corresponding index files and includes them as input if they exist.
    """
    in_dict = {}

    if known_sites and interval:
        logger.info(f"Using known_sites: {known_sites} and interval: {interval}")
        include: "../gatk_bqsr/gatk_bqsr.smk"
        in_dict["normal_bam"] = f"{indir}/bqsr/{wildcards.normal_sample_id}.sorted.markdup.BQSR.bam"
        in_dict["experimental_bam"] = f"{indir}/bqsr/{wildcards.experimental_sample_id}.sorted.markdup.BQSR.bam"
    else:
        logger.info("No known_sites or interval specified in config, proceeding without them.")
        in_dict["normal_bam"] = f"{indir}/bam-sorted-Markdup/{wildcards.normal_sample_id}.bam"
        in_dict["experimental_bam"] = f"{indir}/bam-sorted-Markdup/{wildcards.experimental_sample_id}.bam"

    if fai_index and dict_index:
        logger.info(f"Using provided fai_index: {fai_index} and dict_index: {dict_index}")
        in_dict["fai"] = fai_index
        in_dict["dict"] = dict_index
        in_dict["fasta"] = fasta
    else:
        logger.info("No fai_index or dict_index specified in config, using rule to generate them.")
        in_dict["fai"] = f"{outdir}/index/genome.fa.fai"
        in_dict["dict"] = f"{outdir}/index/genome.dict"
        in_dict["fasta"] = f"{outdir}/index/genome.fa"
    return in_dict


rule somaticMutect2:
    input:
       unpack(get_input_somaticMutect2)
    output:
        vcf = outdir + "/mutect2-vcf/{normal_sample_id}_{experimental_sample_id}.vcf.gz"
    log:
        logdir + "/all/gatk/{normal_sample_id}_{experimental_sample_id}/mutect2.log"
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
