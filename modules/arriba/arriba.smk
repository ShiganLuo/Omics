include: "../common/common.smk"

indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
samples = config.get("samples",[])
bam_substring = config.get("bam_substring") or ""
# ① parse 期读入 + 规范化(str 收成单元素 list;合法值由 schema 声明)。
# output 指令无法调用函数,格式列表必须是模块级变量。
_formats = config.get("Params", {}).get("arriba", {}).get("formats") or ["png"]
if isinstance(_formats, str):
    _formats = [_formats]
# figure base names (no extension) — expanded per format below
FUSION_FIG_BASES = [
    "fig1_per_sample_counts",
    "fig2_fusion_type_distribution",
    "fig3_type_by_sample",
    "fig4_reading_frame",
    "fig5_recurrent_heatmap",
    "fig6_inframe_support",
]

def get_input_for_arriba(wildcards):
    in_dict = {}
    if bam_substring != "" :
        in_dict["bam"] = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.{bam_substring}.bam"
    else:
        in_dict["bam"] = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    in_dict["fasta"] = config.get('genome',{}).get('fasta')
    in_dict["gtf"] = config.get('genome',{}).get('gtf')
    in_dict["blacklist"] = config.get('Params',{}).get('arriba',{}).get('blacklist')
    in_dict["known_fusions"] = config.get('Params',{}).get('arriba',{}).get('known_fusions')
    return in_dict
rule arriba:
    input:
        unpack(get_input_for_arriba)
    output:
        passed_fusion_tsv = outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv",
        discarded_fusion_tsv = outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv",
    log:
        logdir + "/{sample_id}/arriba.log"
    threads: rule_res("arriba").get("threads", 4)
    benchmark: logdir + "/{sample_id}/arriba.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "arriba", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        arriba = config.get('Procedure',{}).get('arriba') or 'arriba',
        t = config.get('Params',{}).get('arriba',{}).get('t') or None,
        d = config.get('Params',{}).get('arriba',{}).get('d') or None,
        E = config.get('Params',{}).get('arriba',{}).get('E') or 0.3,
        p = config.get('Params',{}).get('arriba',{}).get('p') or None
    conda:
        "arriba.yaml"
    container:
        sif("arriba.yaml", "arriba")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("arriba", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = os.path.join(RECEIPT_DIR, f"{wildcards.sample_id}_arriba_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = f"{outdir}/{wildcards.sample_id}/{wildcards.sample_id}_arriba.sh"
            if not os.path.exists(input.blacklist):
                rule_logger.error(f"Blacklist file not found: {input.blacklist}")
                raise ValueError(f"Blacklist file not found: {input.blacklist}")
            if not os.path.exists(input.known_fusions):
                rule_logger.error(f"Known fusions file not found: {input.known_fusions}")
                raise ValueError(f"Known fusions file not found: {input.known_fusions}")
            cmd = [
                params.arriba,
                "-x", input.bam,
                "-o", output.passed_fusion_tsv,
                "-O", output.discarded_fusion_tsv,
                "-a", input.fasta,
                "-g", input.gtf,
                "-b", input.blacklist,
                "-k", input.known_fusions,
                "-E", str(params.E)
            ]
            if params.t:
                cmd += ["-t", params.t]
            if params.d:
                cmd += ["-d", params.d]
            if params.p:
                cmd += ["-p", params.p]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Arriba fusion detection for sample {wildcards.sample_id} successfully completed!"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Arriba fusion detection failed for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")
            raise RuntimeError(f"Arriba fusion detection failed for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

rule arriba_report:
    input:
        passed_fusions = expand(outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv", sample_id=samples),
        discarded_fusions = expand(outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv", sample_id=samples)
    output:
        report = outdir + "/arriba_report/arriba_fusion_report.html",
        fusion_summary = outdir + "/arriba_report/per_sample_summary.tsv",
        recurrent_fusions = outdir + "/arriba_report/recurrent_fusions.tsv",
        high_medium_fusions = outdir + "/arriba_report/high_medium_confidence_fusions.tsv",
        inframe_fusions = outdir + "/arriba_report/inframe_fusions.tsv",
        # ② 按格式展开:每种格式一个输出文件(base 不带扩展名)
        fusion_figs = expand(outdir + "/arriba_report/figures/{base}.{fmt}",
                             base=FUSION_FIG_BASES, fmt=_formats),
    log:
        logdir + "/arriba/arriba_report.log"
    benchmark: logdir + "/arriba/arriba_report.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "arriba_report", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        summary_script = os.path.join(ROOT_DIR, "modules/arriba/bin/summarize_arriba_fusions.py"),
        formats = _formats,   # ③ 转发
    conda:
        "arriba.yaml"
    container:
        sif("arriba.yaml", "arriba_report")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("arriba_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = os.path.join(RECEIPT_DIR, f"arriba_report_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = f"{outdir}/arriba_report.sh"
            cmd = [
                "python", params.summary_script,
                "-p", ",".join(input.passed_fusions),
                "-d", ",".join(input.discarded_fusions),
                "-o", outdir + "/arriba_report"
            ]
            # script's -f/--format is nargs="+" — single occurrence, all values
            cmd += ["--format"] + [str(fmt) for fmt in params.formats]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Arriba report generation successfully completed!"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Arriba report generation failed: {e}\n")
            raise RuntimeError(f"Arriba report generation failed: {e}\n")
rule arriba_result:
    input:
        passed_fusion_tsv = outdir + "/{sample_id}/{sample_id}_passed_fusions.tsv",
        discarded_fusion_tsv = outdir + "/{sample_id}/{sample_id}_discarded_fusions.tsv"
