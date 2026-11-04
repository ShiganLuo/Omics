from snakemake.logging import logger
include: "../../common/common.smk"

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fai_index = config.get("genome", {}).get("fai_index")
dict_index = config.get("genome", {}).get("dict_index")
fasta = config.get("genome", {}).get("fasta")
targets = config.get("genome", {}).get("targets")
known_sites = config.get("genome", {}).get("known_sites")
interval = config.get("genome", {}).get("interval")
# Germline suppression resources for Mutect2/FilterMutectCalls (both optional):
#   genome.germline_resource - population allele-frequency VCF (e.g. gnomAD),
#                       assembly-tied so it lives in genome.references.<name>
#   Params.gatk.Mutect2.panel_of_normals - sites observed in a SELF-BUILT
#                       normal panel (cohort/lab property, not a genome
#                       property); recurrent artifacts get hard-filtered
germline_resource = config.get("genome", {}).get("germline_resource")
panel_of_normals = config.get("Params", {}).get("gatk", {}).get("Mutect2", {}).get("panel_of_normals")

def get_input_somaticMutect2(wildcards):
    """
    This function determines the appropriate input BAM file for the HaplotypeCaller step based on the presence of known_sites and interval parameters.
    If both known_sites and interval are provided, it assumes that BQSR has been performed and uses the BQSR-corrected BAM file. Otherwise, it falls back to using the MarkDuplicates BAM file. 
    It also checks for the presence of the corresponding index files and includes them as input if they exist.
    """
    in_dict = {}

    if known_sites and interval:
        logger.debug(f"Using known_sites: {known_sites} and interval: {interval}")
        include: "../gatk_bqsr/gatk_bqsr.smk"
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bqsr.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bqsr.bam"
    else:
        logger.debug("No known_sites or interval specified in config, proceeding without them.")
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bam"

    if fai_index and dict_index:
        logger.debug(f"Using provided fai_index: {fai_index} and dict_index: {dict_index}")
        in_dict["fai"] = fai_index
        in_dict["dict"] = dict_index
        in_dict["fasta"] = fasta
    else:
        logger.debug("No fai_index or dict_index specified in config, using rule to generate them.")
        in_dict["fai"] = f"{indir}/index/genome.fa.fai"
        in_dict["dict"] = f"{indir}/index/genome.dict"
        in_dict["fasta"] = f"{indir}/index/genome.fa"
    return in_dict

rule somaticMutect2:
    input:
       unpack(get_input_somaticMutect2),
       targets = lambda wildcards: [targets] if targets else [],
       germline_resource = lambda wildcards: [germline_resource] if germline_resource else [],
       panel_of_normals = lambda wildcards: [panel_of_normals] if panel_of_normals else [],
    output:
        vcf = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/{normal_sample_id}_vs_{experimental_sample_id}.vcf.gz"
    log:
        logdir + "/{normal_sample_id}_{experimental_sample_id}_mutect2.log"
    conda:
        "../gatk.yaml"
    container:
        sif("../gatk.yaml", "somaticMutect2")
    params:
        javaOptions =  config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("tmp_dir") or None
    threads: rule_res("somaticMutect2").get("threads", 10)
    resources:
        mem_mb = rule_res("somaticMutect2").get("mem_mb", 32000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="somaticMutect2", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"somaticMutect2_{wildcards.normal_sample_id}_{current_time}.sh")
            sample_outdir = os.path.dirname(output.vcf)
            plain_script = os.path.join(sample_outdir, "somaticMutect2.sh")
            cmd = [
                params.gatk, "Mutect2",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-I", input.normal_bam,
                "-I", input.experimental_bam,
                "-normal", wildcards.normal_sample_id,
                "-O", output.vcf,
                "--native-pair-hmm-threads", str(threads)
            ]
            if input.targets:
                cmd.extend(["-L", str(input.targets[0])])
            # Germline suppression: population AF resource + normal panel.
            # Both flow through to FilterMutectCalls via the raw VCF header.
            if input.germline_resource:
                cmd.extend(["--germline-resource", str(input.germline_resource[0])])
            if input.panel_of_normals:
                cmd.extend(["--pon", str(input.panel_of_normals[0])])
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'somaticMutect2 completed successfully for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during somaticMutect2 execution: {str(e)}")
            raise RuntimeError(f"Error occurred while running somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id}: {e}, you can check the log file {log_path} for more details.")

