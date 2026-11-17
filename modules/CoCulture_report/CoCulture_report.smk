include: "../common/common.smk"

indir = config.get("indir", "data")
outdir = config.get("outdir", "results")
logdir = config.get("logdir", "logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
genome1_samples = config.get("genome1_samples", [])
genome2_samples = config.get("genome2_samples", [])
genome1, genome2 = config.get("genome_pairs", ["genome1", "genome2"])
# formats(模块 md #9 / subworkflow.md formats 规范):matplotlib savefig 格式
# 集,一种格式一个输出文件;str 收成单元素列表(格式合法性由 schema 声明)
_formats = config.get("formats") or ["png"]
if isinstance(_formats, str):
    _formats = [_formats]
# 出图基名(不含扩展名;每种格式一个文件)。bin/HRT.py 只为样本列表非空的
# 基因组出图,所以 output 的格式展开以列表为条件(module 级变量——output
# 指令无法调用函数,见 decoil.smk 注释)
_plot_bases = []
if genome1_samples:
    _plot_bases.append(genome1 + "_2_" + genome2)
if genome2_samples:
    _plot_bases.append(genome2 + "_2_" + genome1)
# need for test
rule HRT:
    input:
        genome1_count = indir + "/" + genome1 + "/all_TEcount_name.tsv",
        genome2_count = indir + "/" + genome2 + "/all_TEcount_name.tsv"
    output:
        report = outdir + "/CoCulture_report/CoCulture_report.html",
        # ② 按格式展开:每种格式一个输出文件(base 不带扩展名)
        plots = [outdir + "/CoCulture_report/" + base + "." + fmt
                 for base in _plot_bases for fmt in _formats],
    log:
        logdir + "/CoCulture_report/CoCulture_report.log"
    benchmark: logdir + "/CoCulture_report/CoCulture_report.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "HRT", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        genome1 = genome1,
        genome2 = genome2,
        genome1_samples = genome1_samples,
        genome2_samples = genome2_samples,
        formats = _formats,   # ③ 转发
        HTR_script = ROOT_DIR + "/modules/CoCulture_report/HRT.py"
    conda:
        "CoCulture_report.yaml"
    container:
        sif("CoCulture_report.yaml", "HRT")
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        logger = setup_logger(logger_name="HRT_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start HRT run at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"HRT_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"HRT.sh")
            cmd = ["python", params.HTR_script,
                   "--genome1_count", input.genome1_count,
                   "--genome2_count", input.genome2_count,
                   "--genome1", params.genome1,
                   "--genome2", params.genome2,
                   "--genome1_samples", ' '.join(params.genome1_samples),
                   "--genome2_samples", ' '.join(params.genome2_samples),
                   "--output", output.report]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")
