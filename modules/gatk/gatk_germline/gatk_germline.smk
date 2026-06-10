from snakemake.logging import logger
import time
import os
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fasta = config.get("genome", {}).get("fasta")
fai_index = config.get("genome", {}).get("fai_index")
dict_index = config.get("genome", {}).get("dict_index")
interval = config.get("reference", {}).get("interval")
known_sites = config.get("reference", {}).get("known_sites", [])

def get_input_HaplotypeCaller(wildcards):
    """
    This function determines the appropriate input BAM file for the HaplotypeCaller step based on the presence of known_sites and interval parameters.
    If both known_sites and interval are provided, it assumes that BQSR has been performed and uses the BQSR-corrected BAM file. Otherwise, it falls back to using the MarkDuplicates BAM file. 
    It also checks for the presence of the corresponding index files and includes them as input if they exist.
    """
    in_dict = {}
    if known_sites and interval:
        logger.info(f"Using known_sites: {known_sites} and interval: {interval}")
        include: "bqsr.smk"
        in_dict["bam"] = f"{indir}/bqsr/wildcards.sample_id/{wildcards.sample_id}.sorted.markdup.BQSR.bam"
    else:
        logger.info("No known_sites or interval specified in config, proceeding without them.")
        in_dict["bam"] = f"{indir}/bam-sorted-Markdup/{wildcards.sample_id}/{wildcards.sample_id}.bam"

    if fai_index and dict_index and fasta:
        logger.info(f"Using provided fai_index: {fai_index} , dict_index: {dict_index}")
        in_dict["fai"] = fai_index
        in_dict["dict"] = dict_index
        in_dict["fasta"] = fasta
    else:
        logger.info("No fai_index or dict_index specified in config, using rule to generate them.")
        in_dict["fai"] = f"{indir}/index/genome.fa.fai"
        in_dict["dict"] = f"{indir}/index/genome.dict"
        in_dict["fasta"] = f"{indir}/index/genome.fa"
    return in_dict

rule HaplotypeCaller:
    input:
        unpack(get_input_HaplotypeCaller)
    output:
        vcf = outdir + "/vcf/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/haplotypeCaller.log"
    conda:
        "../gatk.yaml"
    params:
        javaOptions =  "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk"
    threads: 10
    run:
        current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        logger.info(f"Start HaplotypeCaller for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir,f"HaplotypeCaller_{current_time}.sh")
        cmd = [
            params.gatk, "--java-options", params.javaOptions, "HaplotypeCaller",
            "-R", input.fasta,
            "-I", input.bam,
            "-O", output.vcf
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")

def get_input_filterHaplotypeCallerVcf(wildcards):
    in_dict = {}
    in_dict["vcf"] = f"{outdir}/vcf/{wildcards.sample_id}/{wildcards.sample_id}.vcf.gz"
    if fai_index and dict_index and fasta:
        logger.info(f"Using provided fai_index: {fai_index} and dict_index: {dict_index}")
        in_dict["fai"] = fai_index
        in_dict["dict"] = dict_index
        in_dict["fasta"] = fasta
    else:
        logger.info("No fai_index or dict_index specified in config, using rule to generate them.")
        in_dict["fai"] = f"{indir}/index/genome.fa.fai"
        in_dict["dict"] = f"{indir}/index/genome.dict"
        in_dict["fasta"] = f"{indir}/index/genome.fa"
    return in_dict
rule filterHaplotypeCallerVcf:
    input:
        unpack(get_input_filterHaplotypeCallerVcf)
    output:
        vcf = outdir + "/vcf-filtered/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/haplotypeCaller-filtered.log"
    conda:
        "../gatk.yaml"
    threads: 10
    params:
        javaOptions = "-Xmx35G", # "--java-options -Xmx15G",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        FS_threshold = config.get("Params", {}).get("gatk", {}).get("FS_threshold") or 20.0,
        QD_threshold = config.get("Params", {}).get("gatk", {}).get("QD_threshold") or 2.0,
        DP_threshold = config.get("Params", {}).get("gatk", {}).get("DP_threshold") or 10.0,
        QUAL_threshold = config.get("Params", {}).get("gatk", {}).get("QUAL_threshold") or 20.0
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" VariantFiltration \
            -R {input.fasta} -V {input.vcf} -O {output.vcf} \
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
