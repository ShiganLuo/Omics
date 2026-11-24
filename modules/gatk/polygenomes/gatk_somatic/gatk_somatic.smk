include: "../../../common/common.smk"

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
# 参数精确到 rule:javaOptions/tmp-dir 按 rule 独立配置(Params.gatk.<rule_name>
# .javaOptions/.tmp-dir),无 module 级共享。

def get_input_somaticMutect2(wildcards):
    """
    This function determines the appropriate input BAM file for the HaplotypeCaller step based on the presence of known_sites and interval parameters.
    If both known_sites and interval are provided, it assumes that BQSR has been performed and uses the BQSR-corrected BAM file. Otherwise, it falls back to using the MarkDuplicates BAM file.
    It also checks for the presence of the corresponding index files and includes them as input if they exist.
    """
    in_dict = {}
    known_sites = config.get("genome", {}).get(wildcards.genome, {}).get("known_sites")
    interval = config.get("genome", {}).get(wildcards.genome, {}).get("interval")

    if known_sites and interval:
        logger.debug(f"Using known_sites: {known_sites} and interval: {interval}")
        include: "../gatk_bqsr/gatk_bqsr.smk"
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bqsr.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bqsr.bam"
    else:
        logger.debug("No known_sites or interval specified in config, proceeding without them.")
        in_dict["normal_bam"] = f"{indir}/{wildcards.normal_sample_id}/{wildcards.normal_sample_id}.sorted_markdup.bam"
        in_dict["experimental_bam"] = f"{indir}/{wildcards.experimental_sample_id}/{wildcards.experimental_sample_id}.sorted_markdup.bam"

    fai_index = config.get("genome", {}).get(wildcards.genome, {}).get("fai_index")
    dict_index = config.get("genome", {}).get(wildcards.genome, {}).get("dict_index")
    fasta = config.get("genome", {}).get(wildcards.genome, {}).get("fasta")
    if fai_index and dict_index:
        logger.debug(f"Using provided fai_index: {fai_index} and dict_index: {dict_index}")
        in_dict["fai"] = fai_index
        in_dict["dict"] = dict_index
        in_dict["fasta"] = fasta
    else:
        logger.debug("No fai_index or dict_index specified in config, using rule to generate them.")
        in_dict["fai"] = f"{indir}/index/{wildcards.genome}/genome.fa.fai"
        in_dict["dict"] = f"{indir}/index/{wildcards.genome}/genome.dict"
        in_dict["fasta"] = f"{indir}/index/{wildcards.genome}/genome.fa"
    return in_dict

rule somaticMutect2:
    input:
       unpack(get_input_somaticMutect2)
    output:
        vcf = outdir + "/{genome}/mutect2-vcf/{normal_sample_id}_vs_{experimental_sample_id}/{normal_sample_id}_vs_{experimental_sample_id}.vcf.gz"
    log:
        logdir + "/{genome}/all/gatk/{normal_sample_id}_{experimental_sample_id}/mutect2.log"
    conda:
        "../../gatk.yaml"
    container:
        sif("../../gatk.yaml", "somaticMutect2")
    params:
        javaOptions = config.get("Params", {}).get("gatk", {}).get("somaticMutect2", {}).get("javaOptions", "-Xmx30g"),
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("somaticMutect2", {}).get("tmp-dir", None),
        parameters = config.get("mutect2_parameters") or ""
    threads: rule_res("somaticMutect2").get("threads", 10)
    benchmark: logdir + "/{genome}/all/gatk/{normal_sample_id}_{experimental_sample_id}/mutect2.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "somaticMutect2", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 32000),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger(logger_name="somaticMutect2", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} genome {wildcards.genome} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"mutect2_{wildcards.normal_sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"{wildcards.genome}/mutect2-vcf/{wildcards.normal_sample_id}_vs_{wildcards.experimental_sample_id}/mutect2.sh")
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
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} on genome {wildcards.genome} successfully completed!"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during somaticMutect2 execution: {str(e)}")
            raise RuntimeError(f"Error occurred while running somaticMutect2 for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.")
