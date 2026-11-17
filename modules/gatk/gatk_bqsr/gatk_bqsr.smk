include: "../../common/common.smk"

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fasta = config.get("reference", {}).get("fasta")
interval = config.get("reference", {}).get("interval")
known_sites = config.get("reference", {}).get("known_sites", [])
# 参数精确到 rule:javaOptions/tmp-dir 按 rule 独立配置(Params.gatk.<rule_name>
# .javaOptions/.tmp-dir),无 module 级共享;rule_param 见 common.smk。
_gatk_params = config.get("Params", {}).get("gatk", {})

rule BaseRecalibrator:
    input:
        bam = indir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = indir + "/{sample_id}/{sample_id}.sorted_markdup.bai",
        ref = fasta,
        interval = interval
    output:
        table = outdir + "/bqsr/{sample_id}.recal_data.table"
    log:
        logdir + "/{sample_id}/gatk_bqsr.log"
    threads: rule_res("BaseRecalibrator").get("threads", 8)
    benchmark: logdir + "/{sample_id}/gatk_bqsr.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "BaseRecalibrator", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../gatk.yaml"
    container:
        sif("../gatk.yaml", "BaseRecalibrator")
    params:
        javaOptions = rule_param(_gatk_params, "BaseRecalibrator", "javaOptions", "-Xmx30g"),
        tmp_dir = rule_param(_gatk_params, "BaseRecalibrator", "tmp-dir", None),
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        known_sites = " ".join([f"--known-sites {ks}" for ks in known_sites])
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="BaseRecalibrator", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start BaseRecalibrator for sample {wildcards.sample_id} at {current_time}")
            os.makedirs(os.path.dirname(str(output.table)), exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"BaseRecalibrator_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/BaseRecalibrator.sh")
            os.makedirs(os.path.dirname(plain_script), exist_ok=True)
            cmd = [
                params.gatk, "BaseRecalibrator",
                "--java-options", params.javaOptions,
                "-R", input.ref,
                "-I", input.bam,
                params.known_sites
            ]
            if input.interval:
                cmd.extend(["-L", input.interval])
            cmd.extend([
                "-O", output.table
            ])
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'BaseRecalibrator completed successfully for sample {wildcards.sample_id}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during BaseRecalibrator execution: {str(e)}")
            raise RuntimeError(f"Error occurred while running BaseRecalibrator for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

rule ApplyBQSR:
    input:
        bam = indir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        table = outdir + "/bqsr/{sample_id}.recal_data.table",
        ref = fasta
    output:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam"
    log:
        logdir + "/{sample_id}/ApplyBQSR.log"
    conda:
        "../gatk.yaml"
    container:
        sif("../gatk.yaml", "ApplyBQSR")
    threads:
        8
    benchmark: logdir + "/{sample_id}/ApplyBQSR.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "ApplyBQSR", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 16000),
    params:
        javaOptions = rule_param(_gatk_params, "ApplyBQSR", "javaOptions", "-Xmx30g"),
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        tmp_dir = rule_param(_gatk_params, "ApplyBQSR", "tmp-dir", None),
        interval = interval
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger(logger_name="ApplyBQSR", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ApplyBQSR for sample {wildcards.sample_id} at {current_time}")
            os.makedirs(os.path.dirname(str(output.bam)), exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"ApplyBQSR_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/ApplyBQSR.sh")
            os.makedirs(os.path.dirname(plain_script), exist_ok=True)
            cmd = [
                params.gatk, "ApplyBQSR",
                "--java-options", params.javaOptions,
                "-R", input.ref,
                "-I", input.bam,
                "-bqsr", input.table,
                "-O", output.bam
            ]
            if params.tmp_dir:
                cmd.extend(["--tmp-dir", params.tmp_dir])
            if params.interval:
                cmd.extend(["-L", params.interval])
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'ApplyBQSR completed successfully for sample {wildcards.sample_id}'\n")
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during ApplyBQSR execution: {str(e)}")
            raise RuntimeError(f"Error occurred while running ApplyBQSR for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

rule gatk_bqsr_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam.bai"

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "gatk_bqsr_result", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
