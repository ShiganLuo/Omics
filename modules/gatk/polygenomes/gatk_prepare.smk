from snakemake.logging import logger
include: "../../common/common.smk"

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
logdir_combine = config.get("logdir_combine") or "log"
input_bam_substring = config.get("input_bam_substring") or ""

def get_input_for_gatk_index(wildcards):
    logger.debug(f"gatk_index called with {wildcards}")
    in_dict = {}
    fasta = config.get("genome", {}).get('references', {}).get(wildcards.genome, {}).get("fasta")
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"Fasta file for genome {wildcards.genome} not found in config or does not exist: {fasta}")
    return fasta
rule gatk_index:
    input:
        fasta = get_input_for_gatk_index
    output:
        fai_index = outdir + "/index/{genome}/{genome}.fa.fai",
        dict_index = outdir + "/index/{genome}/{genome}.dict",
        fasta_link = outdir + "/index/{genome}/{genome}.fa"
    log:
        logdir_combine + "/index/{genome}/gatk_index.log"
    threads: rule_res("gatk_index").get("threads", 4)
    resources:
        mem_mb=rule_res("gatk_index").get("mem_mb", 8000),
    conda: "../gatk.yaml"
    container:
        sif("../gatk.yaml")
    params:
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="gatk_index", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start gatk_index for genome {wildcards.genome} at {current_time}")
            index_outdir = os.path.join(outdir, "index", wildcards.genome)
            os.makedirs(index_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"gatk_index_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(index_outdir, f"gatk_index.sh")
            cmd1 = ["ln", "-s", input.fasta, output.fasta_link]
            cmd2 = [
                params.gatk, "CreateSequenceDictionary",
                "--java-options", params.javaOptions,
                "-R", input.fasta,
                "-O", output.dict_index
            ]
            if params.tmp_dir:
                cmd2.extend(["--tmp-dir", params.tmp_dir])
            cmd3 = [
                "samtools", "faidx",
                input.fasta,
                "-o", output.fai_index
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(" ".join(cmd3) + "\n")
                f.write(f'echo "gatk_index for genome {wildcards.genome} successfully completed!"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred while running gatk_index for genome {wildcards.genome}: {e}")
            raise RuntimeError(f"Error occurred while running gatk_index for genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.")

def get_input_for_addReadsGroup(wildcards):
    logger.debug(f"addReadsGroup called with {wildcards}")
    if input_bam_substring != "":
        bam = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}." + input_bam_substring + ".bam"
    else:
        bam = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    return bam

rule addReadsGroup:
    input:
        bam = get_input_for_addReadsGroup
    output:
        bam = temp(outdir + "/{genome}/{sample_id}/{sample_id}.addReadsGroup.bam"),
        bai = temp(outdir + "/{genome}/{sample_id}/{sample_id}.addReadsGroup.bai")
    log:
        logdir + "/{sample_id}/{genome}/addReadsGroup.log"
    conda: "../gatk.yaml"
    container:
        sif("../gatk.yaml")
    threads: rule_res("addReadsGroup").get("threads", 8)
    resources:
        mem_mb=rule_res("addReadsGroup").get("mem_mb", 32000),
    params:
        id = "{sample_id}",
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None,
        RGLB = config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGLB") or "lib1",
        RGPL = config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGPL") or "illumina",
        RGPU = config.get("Params", {}).get("gatk", {}).get("addReadsGroup", {}).get("RGPU") or "unit1",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="addReadsGroup", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start addReadsGroup for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.genome, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"addReadsGroup_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"addReadsGroup.sh")
            cmd1 = [
                params.gatk, "AddOrReplaceReadGroups", "--java-options", f'"{params.javaOptions}"',
                "--INPUT", input.bam,
                "--OUTPUT", output.bam,
                "-SO", "coordinate",
                "--RGLB", params.RGLB,
                "--RGPL", params.RGPL,
                "--RGPU", params.RGPU,
                "--RGSM", params.id
            ]
            if params.tmp_dir:
                cmd1.extend(["--tmp-dir", params.tmp_dir])
            cmd2 = [
                params.samtools, "index",
                "-@", str(threads),
                output.bam,
                "-o", output.bai
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
                f.write(f'echo "addReadsGroup for sample {wildcards.sample_id} on genome {wildcards.genome} successfully completed!"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred while running addReadsGroup for sample {wildcards.sample_id} genome {wildcards.genome}: {e}.")
            raise RuntimeError(f"Error occurred while running addReadsGroup for sample {wildcards.sample_id} genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.")

rule MarkDuplicates:
    input:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.addReadsGroup.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.addReadsGroup.bai"
    output:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.sorted_markdup.bai",
        metrics = outdir + "/{genome}/{sample_id}/{sample_id}.Markdup-metrics.txt"
    log:
        logdir + "/{sample_id}/{genome}/MarkDuplicates.log"
    conda: "../gatk.yaml"
    container:
        sif("../gatk.yaml")
    threads: rule_res("MarkDuplicates").get("threads", 8)
    resources:
        mem_mb=rule_res("MarkDuplicates").get("mem_mb", 32000),
    params:
        javaOptions = config.get("Params", {}).get("gatk", {}).get("javaOptions") or "-Xmx30g",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = config.get("Params", {}).get("gatk", {}).get("tmp-dir") or None
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="MarkDuplicates", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start MarkDuplicates for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.genome, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"MarkDuplicates_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"MarkDuplicates.sh")
            cmd = [
                params.gatk, "MarkDuplicates", "--java-options", f'"{params.javaOptions}"',
                "--INPUT", input.bam,
                "--OUTPUT", output.bam,
                "--CREATE_INDEX", "true",
                "--VALIDATION_STRINGENCY", "SILENT",
                "--METRICS_FILE", output.metrics
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "MarkDuplicates for sample {wildcards.sample_id} on genome {wildcards.genome} successfully completed!"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred while running MarkDuplicates for sample {wildcards.sample_id} genome {wildcards.genome}: {e}.\n")
            raise RuntimeError(f"Error occurred while running MarkDuplicates for sample {wildcards.sample_id} genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.")

rule gatk_prepare_result:
    input:
        bam = outdir + "/{genome}/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = outdir + "/{genome}/{sample_id}/{sample_id}.sorted_markdup.bai",
        metrics = outdir + "/{genome}/{sample_id}/{sample_id}.Markdup-metrics.txt"
