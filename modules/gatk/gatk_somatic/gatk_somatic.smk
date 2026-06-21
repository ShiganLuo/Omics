from snakemake.logging import logger
include: "../../common/common.smk"
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
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bqsr.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bqsr.bam"
    else:
        logger.info("No known_sites or interval specified in config, proceeding without them.")
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bam"

    if fai_index and dict_index:
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


rule somaticMutect2:
    input:
       unpack(get_input_somaticMutect2)
    output:
        vcf = outdir + "/mutect2-vcf/{normal_sample_id}_vs_{experimental_sample_id}/{normal_sample_id}_vs_{experimental_sample_id}.vcf.gz"
    log:
        logdir + "/all/gatk/{normal_sample_id}_{experimental_sample_id}/mutect2.log"
    conda:
        "../gatk.yaml"
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        parameters = config.get("mutect2_parameters") or ""
    threads: 10
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="somaticMutect2", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} at {current_time}")
            script = os.path.join(outdir,f"mutect2-vcf/{wildcards.normal_sample_id}_vs_{wildcards.experimental_sample_id}/mutect2_{current_time}.sh")
            cmd = [
                params.gatk, "Mutect2",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-I", input.normal_bam,
                "-I", input.experimental_bam,
                "-normal", wildcards.normal_sample_id,
                "-O", output.vcf,
                "--native-pair-hmm-threads", str(threads)
            ] + params.parameters.split()
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"Error during somaticMutect2 execution: {str(e)}\n")
            raise f"Error occurred while running somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id}: {e}, you can check the log file {log[0]} for more details."


