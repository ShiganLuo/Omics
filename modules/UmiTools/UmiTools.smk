include: "../common/common.smk"

outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")

# 参数精确到 rule:Params.umi_tools.<rule_name>.<key>(config 全量路径直读;
# 无 module 级共享回退,每 rule 各自独立键)。

rule umi_tools_dedup_for_hisat2:
    input:
        bam = indir + "/{sample_id}.bam"
    output:
        bam = temp(outdir + "/{sample_id}.dedup.bam"),
        bai = temp(outdir + "/{sample_id}.dedup.bam.bai"),
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        method = config.get("Params", {}).get("umi_tools", {}).get("umi_tools_dedup_for_hisat2", {}).get("method") or "unique",
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    threads: rule_res("umi_tools_dedup_for_hisat2").get("threads", 2)
    benchmark: logdir + "/{sample_id}/umi_tools_dedup_run.benchmark.txt",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "umi_tools_dedup_for_hisat2", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "UmiTools.yaml"
    container:
        sif("UmiTools.yaml", "umi_tools_dedup_for_hisat2")
    log:
        log = logdir + "/{sample_id}/umi_tools_dedup_run.txt"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("umi_tools_dedup_for_hisat2", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start umi_tools_dedup_for_hisat2 for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"umi_tools_dedup_for_hisat2_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"umi_tools_dedup_for_hisat2.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.umi_tools} dedup \\\n")
                f.write(f"    --method={params.method} \\\n")
                f.write(f"    -I {input.bam} \\\n")
                f.write(f"    -S {output.bam} \\\n")
                f.write(f"{params.samtools} index {output.bam}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"umi_tools_dedup_for_hisat2 for sample {wildcards.sample_id} completed")
        except Exception as e:
            rule_logger.error(f"Error occurred: {e}")
            raise RuntimeError(f"Error occurred: {e}, you can check the log file {log_path} for more details.")

rule umi_tools_dedup_for_star:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = temp(outdir + "/{sample_id}.dedup.bam"),
        bai = temp(outdir + "/{sample_id}.dedup.bam.bai")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        method = config.get("Params", {}).get("umi_tools", {}).get("umi_tools_dedup_for_star", {}).get("method") or "unique",
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    threads: rule_res("umi_tools_dedup_for_star").get("threads", 2)
    benchmark: logdir + "/{sample_id}/umi_tools_dedup_run.benchmark.txt",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "umi_tools_dedup_for_star", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "UmiTools.yaml"
    container:
        sif("UmiTools.yaml", "umi_tools_dedup_for_star")
    log:
        log = logdir + "/{sample_id}/umi_tools_dedup_run.txt"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("umi_tools_dedup_for_star", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start umi_tools_dedup_for_star for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"umi_tools_dedup_for_star_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"umi_tools_dedup_for_star.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.umi_tools} dedup --method={params.method} \\\n")
                f.write(f"    -I {input.bam} -S {output.bam}\n")
                f.write(f"{params.samtools} index {output.bam}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"umi_tools_dedup_for_star for sample {wildcards.sample_id} completed")
        except Exception as e:
            rule_logger.error(f"Error occurred: {e}")
            raise RuntimeError(f"Error occurred: {e}, you can check the log file {log_path} for more details.")
