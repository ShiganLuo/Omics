include: "../common/common.smk"

outdir = config.get("outdir", "output")
indir = config.get("indir", "input")
logdir = config.get("logdir", "logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
mode = config.get("mode") or None

# 参数精确到 rule:Params.cutadapt.<rule_name>.<key>(config 全量路径直读;
# 无 module 级共享回退,trimming_Paired / trimming_Single 各自独立键)。

def get_cutadapt_input_paired(wildcards):
    """Return paired-end FASTQ inputs, with optional UMI suffix."""
    if mode == "UMI":
        return [
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.umi.fq.gz",
            f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.umi.fq.gz",
        ]
    return [
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz",
        f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz",
    ]

def get_cutadapt_input_single(wildcards):
    """Return single-end FASTQ input, with optional UMI suffix."""
    if mode == "UMI":
        return f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.umi.single.fq.gz"
    return f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

rule trimming_Paired:
    input:
        get_cutadapt_input_paired
    output:
        fastq1 = outdir + "/{sample_id}/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{sample_id}/{sample_id}_2.fq.gz",
        report  = outdir + "/{sample_id}/{sample_id}.cutadapt_report.txt",
    log:
        logdir + "/{sample_id}/cutadapt_paired.log"
    threads: rule_res("trimming_Paired").get("threads", 6)
    benchmark: logdir + "/{sample_id}/cutadapt_paired.benchmark.tsv",
    resources:
        mem_mb=rule_res("trimming_Paired").get("mem_mb", 8000),
    params:
        cutadapt     = config.get("Procedure", {}).get("cutadapt") or "cutadapt",
        quality      = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("quality") or 20,
        adapter_r1   = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("adapter_r1", None),
        adapter_r2   = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("adapter_r2", None),
        minimum_length = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("minimum_length", None),
        maximum_length = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("maximum_length", None),
        match_read_wildcards = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("match_read_wildcards", False),
        cut = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("cut", None),
        trimmed_only = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("trimmed_only", False),
        report = config.get("Params", {}).get("cutadapt", {}).get("trimming_Paired", {}).get("report") or "full"
    conda:
        "cutadapt.yaml"
    container:
        sif("cutadapt.yaml", "trimming_Paired")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        logger = setup_logger(logger_name="cutadapt_paired", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start cutadapt paired trimming for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"cutadapt_paired_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/cutadapt_paired.sh")
            cmd = [
                params.cutadapt,
                "-j", str(threads),
                "-q", str(params.quality),
                "--report", str(params.report),
                "-o", output.fastq1,
                "-p", output.fastq2,
            ]
            if params.adapter_r1:
                for adapter in params.adapter_r1:
                    cmd += ["-a", adapter]
            if params.adapter_r2:
                for adapter in params.adapter_r2:
                    cmd += ["-A", adapter]
            if params.minimum_length:
                cmd += ["-m", str(params.minimum_length)]
            if params.maximum_length:
                cmd += ["-M", str(params.maximum_length)]
            if params.match_read_wildcards:
                cmd += ["--match-read-wildcards"]
            if params.cut:
                cmd += ["--cut", str(params.cut)]
            if params.trimmed_only:
                cmd += ["--trimmed-only"]
            cmd += [str(input[0]), str(input[1])]
            cmd_str = " ".join(cmd)
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(cmd_str + "\n")
                f.write(f"echo 'cutadapt paired trimming finished for {wildcards.sample_id}'\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} > {output.report} 2>> {log_path}")
            logger.info(f"Finished cutadapt paired trimming for sample {wildcards.sample_id}")
        except Exception as e:
            raise

rule trimming_Single:
    input:
        get_cutadapt_input_single
    output:
        fastq  = outdir + "/{sample_id}/{sample_id}.single.fq.gz",
        report = outdir + "/{sample_id}/{sample_id}.single.cutadapt_report.txt",
    log:
        logdir + "/{sample_id}/cutadapt_single.log"
    threads: rule_res("trimming_Single").get("threads", 6)
    benchmark: logdir + "/{sample_id}/cutadapt_single.benchmark.tsv",
    resources:
        mem_mb=rule_res("trimming_Single").get("mem_mb", 4000),
    params:
        cutadapt       = config.get("Procedure", {}).get("cutadapt") or "cutadapt",
        quality        = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("quality") or 20,
        adapter_r1     = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("adapter_r1", None),
        minimum_length = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("minimum_length", None),
        maximum_length = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("maximum_length", None),
        match_read_wildcards = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("match_read_wildcards", False),
        cut = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("cut", None),
        trimmed_only = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("trimmed_only", False),
        report = config.get("Params", {}).get("cutadapt", {}).get("trimming_Single", {}).get("report") or "full"
    conda:
        "cutadapt.yaml"
    container:
        sif("cutadapt.yaml", "trimming_Single")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger(logger_name="cutadapt_single", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cutadapt single trimming for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"cutadapt_single_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/cutadapt_single.sh")
            cmd = [
                params.cutadapt,
                "-j", str(threads),
                "-q", str(params.quality),
                "--report", str(params.report),
                "-o", output.fastq,
            ]
            if params.adapter_r1:
                for adapter in params.adapter_r1:
                    cmd += ["-a", adapter]
            if params.minimum_length:
                cmd += ["-m", str(params.minimum_length)]
            if params.maximum_length:
                cmd += ["-M", str(params.maximum_length)]
            if params.match_read_wildcards:
                cmd += ["--match-read-wildcards"]
            if params.cut:
                cmd += ["--cut", str(params.cut)]
            if params.trimmed_only:
                cmd += ["--trimmed-only"]
            cmd += [str(input[0])]
            cmd_str = " ".join(cmd)
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(cmd_str + "\n")
                f.write(f"echo 'cutadapt single trimming finished for {wildcards.sample_id}'\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} > {output.report} 2>> {log_path}")
        except Exception as e:
            rule_logger.error()
            raise
