include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", "log")
indir = config.get("indir", "input")
bam_suffix = config.get("bam_suffix", "")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

# ① parse 期读入 + 规范化（str 收成单元素 list；合法值由 schema 声明）
# 参数精确到 rule:formats 挂 plotFragmentSize 名下
# (Params.fragment_size.plotFragmentSize.formats,单 plot rule 消费;
# 无 module 级共享回退),config 全量路径读取 + str 收成单元素 list。
# 每 rule 一个 parse 期变量;fragment_size_result 聚合 plotFragmentSize
# 输出文件名,展开同一 rule 键。
_formats_plotFragmentSize = config.get("Params", {}).get("fragment_size", {}).get("plotFragmentSize", {}).get("formats") or ["png"]
if isinstance(_formats_plotFragmentSize, str):
    _formats_plotFragmentSize = [_formats_plotFragmentSize]
_formats_fragment_size_result = config.get("Params", {}).get("fragment_size", {}).get("plotFragmentSize", {}).get("formats") or ["png"]
if isinstance(_formats_fragment_size_result, str):
    _formats_fragment_size_result = [_formats_fragment_size_result]

rule samtools_stats:
    """Get BAM statistics including fragment size distribution."""
    input:
        bam = indir + "/{sample_id}/{sample_id}" + bam_suffix + ".bam"
    output:
        stats = outdir + "/stats/{sample_id}/{sample_id}_stats.txt"
    log:
        logdir + "/{sample_id}/samtools_stats.log"
    benchmark: logdir + "/{sample_id}/samtools_stats.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "samtools_stats", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "fragment_size.yaml"
    container:
        sif("fragment_size.yaml", "samtools_stats")
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("samtools_stats", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start samtools_stats for sample {wildcards.sample_id} at {current_time}")
            outdir_sample = os.path.dirname(str(output.stats))
            os.makedirs(outdir_sample, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"samtools_stats_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir_sample, f"samtools_stats.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.samtools} stats {input.bam} > {output.stats}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} > {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during samtools_stats for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule getFragmentSize:
    """Extract and summarize fragment size distribution from all samples."""
    input:
        stats = expand(outdir + "/stats/{sample_id}/{sample_id}_stats.txt", sample_id=samples)
    output:
        hist = outdir + "/fragment/FragmentSize.txt"
    log:
        logdir_combine + "/fragment/getFragmentSize.log"
    benchmark: logdir_combine + "/fragment/getFragmentSize.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "getFragmentSize", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/fragment_size/bin/getFragmentSize.py")
    conda:
        "fragment_size.yaml"
    container:
        sif("fragment_size.yaml", "getFragmentSize")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("getFragmentSize", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start getFragmentSize at {current_time}")
            outdir_fragment = os.path.dirname(str(output.hist))
            os.makedirs(outdir_fragment, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"getFragmentSize_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir_fragment, f"getFragmentSize.sh")
            stats_str = " ".join(input.stats)
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"python {params.script} --input {stats_str} --out {output.hist}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} > {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during getFragmentSize: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule plotFragmentSize:
    """Plot fragment size distribution."""
    input:
        hist = outdir + "/fragment/FragmentSize.txt"
    output:
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        plots = [outdir + "/fragment/FragmentSize." + fmt for fmt in _formats_plotFragmentSize],
    log:
        logdir_combine + "/fragment/plotFragmentSize.log"
    benchmark: logdir_combine + "/fragment/plotFragmentSize.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "plotFragmentSize", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/fragment_size/bin/plotFragmentSize.py"),
        formats = _formats_plotFragmentSize,   # ③ 转发
        out_base = outdir + "/fragment/FragmentSize",
    conda:
        "fragment_size.yaml"
    container:
        sif("fragment_size.yaml", "plotFragmentSize")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("plotFragmentSize", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start plotFragmentSize at {current_time}")
            outdir_fragment = os.path.dirname(str(output.plots[0]))
            script = os.path.join(RECEIPT_DIR, f"plotFragmentSize_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir_fragment, f"plotFragmentSize.sh")
            fmt_args = " ".join(f"--format {fmt}" for fmt in params.formats)
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"python {params.script} --input {input.hist} --output {params.out_base} {fmt_args}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} > {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error during plotFragmentSize: {e}")
            raise RuntimeError(f"Error occurred during plotFragmentSize: {e}, you can check the log file {log_path} for more details.")

rule fragment_size_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        hist = outdir + "/fragment/FragmentSize.txt",
        plots = [outdir + "/fragment/FragmentSize." + fmt for fmt in _formats_fragment_size_result]
