include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
Procedure = config.get("Procedure", {})
gtf = config.get("gtf", "")
bam_suffix = config.get("bam_suffix", ".dedup")
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
treated_ip_samples = config.get("treated_ip_samples", [])
treated_input_samples = config.get("treated_input_samples", [])


def get_input_for_call_exomePeak(wildcards):
    """声明峰调用全部输入文件;缺文件在建图期即失败。"""
    return {
        "ip_bams": [indir + f"/{sample_id}/{sample_id}" + bam_suffix + ".bam" for sample_id in ip_samples],
        "input_bams": [indir + f"/{sample_id}/{sample_id}" + bam_suffix + ".bam" for sample_id in input_samples],
        "gtf": gtf,
    }


def get_input_for_diff_exomePeak(wildcards):
    """声明差异分析全部输入文件;缺文件在建图期即失败。"""
    files = get_input_for_call_exomePeak(wildcards)
    files["treated_ip_bams"] = [indir + f"/{sample_id}/{sample_id}" + bam_suffix + ".bam" for sample_id in treated_ip_samples]
    files["treated_input_bams"] = [indir + f"/{sample_id}/{sample_id}" + bam_suffix + ".bam" for sample_id in treated_input_samples]
    return files


rule diff_exomePeak:
    """差异 exomePeak: 处理 vs 未处理 IP/Input 的差异甲基化分析,产出 gene name 注释表。"""
    input:
        unpack(get_input_for_diff_exomePeak),
    output:
        diff_peak_bed = outdir + "/diff_peaks_gene_names.bed",
        diff_peak_xls = outdir + "/diff_peaks_gene_names.xls",
        sig_diff_bed = outdir + "/sig_diff_peak_gene_names.bed",
        sig_diff_xls = outdir + "/sig_diff_peak_gene_names.xls",
        con_sig_diff_bed = outdir + "/con_sig_diff_peak_gene_names.bed",
        con_sig_diff_xls = outdir + "/con_sig_diff_peak_gene_names.xls"
    log:
        logdir + "/exomePeak/diff_exomePeak.log"
    benchmark: logdir + "/exomePeak/diff_exomePeak.benchmark.tsv",
    threads: rule_res("diff_exomePeak").get("threads", 1)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "diff_exomePeak", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        exomePeak_script = os.path.join(ROOT_DIR, "modules", "exomePeak", "bin", "exomePeak.r"),
        rscript = Procedure.get("Rscript") or "Rscript",
    conda:
        "exomePeak.yaml"
    container:
        sif("exomePeak.yaml", "diff_exomePeak")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("diff_exomePeak", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start exomePeak diff at {current_time}")
            sample_outdir = os.path.dirname(str(output.diff_peak_bed))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = receipt_script(f"diff_exomePeak_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "diff_exomePeak.sh")

            cmd = [
                params.rscript, params.exomePeak_script,
                "--gtf", str(input.gtf),
                "--ip_bams", *[str(x) for x in input.ip_bams],
                "--input_bams", *[str(x) for x in input.input_bams],
                "--treated_ip_bams", *[str(x) for x in input.treated_ip_bams],
                "--treated_input_bams", *[str(x) for x in input.treated_input_bams],
                "--outdir", sample_outdir,
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "exomePeak diff was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during exomePeak diff: {e}\n")
            raise RuntimeError(f"Error occurred during exomePeak diff: {e}, you can check the log file {log_path} for more details.\n")


rule call_exomePeak:
    """峰调用 exomePeak: 单条件 IP/Input 的甲基化峰调用,产出 gene name 注释表。"""
    input:
        unpack(get_input_for_call_exomePeak),
    output:
        all_peak_bed = outdir + "/all_peaks_gene_names.bed",
        all_peak_xls = outdir + "/all_peaks_gene_names.xls",
        con_peaks_bed = outdir + "/con_peaks_gene_names.bed",
        con_peaks_xls = outdir + "/con_peaks_gene_names.xls"
    log:
        logdir + "/exomePeak/call_exomePeak.log"
    benchmark: logdir + "/exomePeak/call_exomePeak.benchmark.tsv",
    threads: rule_res("call_exomePeak").get("threads", 12)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "call_exomePeak", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        exomePeak_script = os.path.join(ROOT_DIR, "modules", "exomePeak", "bin", "exomePeak.r"),
        rscript = Procedure.get("Rscript") or "Rscript",
    conda:
        "exomePeak.yaml"
    container:
        sif("exomePeak.yaml", "call_exomePeak")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("call_exomePeak", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start exomePeak call at {current_time}")
            sample_outdir = os.path.dirname(str(output.all_peak_bed))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = receipt_script(f"call_exomePeak_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "call_exomePeak.sh")

            cmd = [
                params.rscript, params.exomePeak_script,
                "--gtf", str(input.gtf),
                "--ip_bams", *[str(x) for x in input.ip_bams],
                "--input_bams", *[str(x) for x in input.input_bams],
                "--outdir", sample_outdir,
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "exomePeak call was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during exomePeak call: {e}\n")
            raise RuntimeError(f"Error occurred during exomePeak call: {e}, you can check the log file {log_path} for more details.\n")


rule exomePeak_result:
    input:
        call_exomePeak = rules.call_exomePeak.output,
        diff_exomePeak = rules.diff_exomePeak.output,
