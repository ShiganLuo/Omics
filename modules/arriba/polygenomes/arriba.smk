include: "../../common/common.smk"

indir = config.get("indir","data/fastq")
outdir = config.get("outdir","output")
logdir = config.get("logdir","logs")
logdir_combine = config.get("logdir_combine","logs")
ROOT_DIR = config.get("ROOT_DIR", ".")
genome_samples = config.get("genome_samples") or {}
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
        in_dict["bam"] = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.{bam_substring}.bam"
    else:
        in_dict["bam"] = indir + f"/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}.bam"
    fasta = config.get("genome", {}).get('references', {}).get(wildcards.genome, {}).get("fasta")
    if not fasta or not os.path.exists(fasta):
        raise ValueError(f"Fasta file for genome {wildcards.genome} not found in config or does not exist: {fasta}")
    gtf = config.get("genome", {}).get('references', {}).get(wildcards.genome, {}).get("gtf")
    if not gtf or not os.path.exists(gtf):
        raise ValueError(f"GTF file for genome {wildcards.genome} not found in config or does not exist: {gtf}")
    blacklist = config.get('Params',{}).get('arriba',{}).get('blacklist')
    if not blacklist or not os.path.exists(blacklist):
        raise ValueError(f"Blacklist file for genome {wildcards.genome} not found in config or does not exist: {blacklist}")
    known_fusions = config.get('Params',{}).get('arriba',{}).get('known_fusions')
    if not known_fusions or not os.path.exists(known_fusions):
        raise ValueError(f"Known fusions file for genome {wildcards.genome} not found in config or does not exist: {known_fusions}")
    in_dict["fasta"] = fasta
    in_dict["gtf"] = gtf
    in_dict["blacklist"] = blacklist
    in_dict["known_fusions"] = known_fusions
    return in_dict
rule arriba:
    input:
        unpack(get_input_for_arriba)
    output:
        passed_fusion_tsv = outdir + "/{genome}/{sample_id}/{sample_id}_passed_fusions.tsv",
        discarded_fusion_tsv = outdir + "/{genome}/{sample_id}/{sample_id}_discarded_fusions.tsv",
    log:
        logdir + "/{sample_id}/{genome}/arriba.log"
    threads: rule_res("arriba").get("threads", 4)
    benchmark: logdir + "/{sample_id}/{genome}/arriba.benchmark.tsv",
    resources:
        mem_mb=rule_res("arriba").get("mem_mb", 8000),
    params:
        arriba = config.get('Procedure',{}).get('arriba') or 'arriba',
        t = config.get('Params',{}).get('arriba',{}).get('t') or None,
        d = config.get('Params',{}).get('arriba',{}).get('d') or None,
        E = config.get('Params',{}).get('arriba',{}).get('E') or 0.3,
        p = config.get('Params',{}).get('arriba',{}).get('p') or None
    conda:
        "../arriba.yaml"
    container:
        sif("../arriba.yaml", "arriba")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("arriba", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            script = receipt_script(f"{wildcards.sample_id}_arriba_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = f"{outdir}/{wildcards.genome}/{wildcards.sample_id}/{wildcards.sample_id}_arriba.sh"
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
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Arriba fusion detection for sample {wildcards.sample_id} on genome {wildcards.genome} successfully completed!"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Arriba fusion detection failed for sample {wildcards.sample_id} genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.")
            raise RuntimeError(f"Arriba fusion detection failed for sample {wildcards.sample_id} genome {wildcards.genome}: {e}, you can check the log file {log_path} for more details.")

def get_input_for_arriba_report(wildcards):
    logger.debug(f"arriba_report called with {wildcards}")
    in_dict = {}
    passed_fusions = []
    discarded_fusions = []
    for sample_id in genome_samples.get(wildcards.genome, []):
        passed_fusions.append(outdir + f"/{wildcards.genome}/{sample_id}/{sample_id}_passed_fusions.tsv")
        discarded_fusions.append(outdir + f"/{wildcards.genome}/{sample_id}/{sample_id}_discarded_fusions.tsv")
        in_dict["passed_fusions"] = passed_fusions
        in_dict["discarded_fusions"] = discarded_fusions
    if len(passed_fusions) == 0:
        raise ValueError(f"No passed fusion files found for genome {wildcards.genome}. genomes_samples: {genome_samples}")
    if len(discarded_fusions) == 0:
        raise ValueError(f"No discarded fusion files found for genome {wildcards.genome}. genomes_samples: {genome_samples}")
    return in_dict
rule arriba_report:
    input:
        unpack(get_input_for_arriba_report)
    output:
        report = outdir + "/{genome}/arriba_report/arriba_fusion_report.html",
        fusion_summary = outdir + "/{genome}/arriba_report/per_sample_summary.tsv",
        recurrent_fusions = outdir + "/{genome}/arriba_report/recurrent_fusions.tsv",
        high_medium_fusions = outdir + "/{genome}/arriba_report/high_medium_confidence_fusions.tsv",
        inframe_fusions = outdir + "/{genome}/arriba_report/inframe_fusions.tsv",
        # ② 按格式展开:每种格式一个输出文件(base 不带扩展名)
        fusion_figs = expand(outdir + "/{genome}/arriba_report/figures/{base}.{fmt}",
                             base=FUSION_FIG_BASES, fmt=_formats, genome="{genome}"),
    log:
        logdir_combine + "/arriba/{genome}/arriba_report.log"
    benchmark: logdir_combine + "/arriba/{genome}/arriba_report.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "arriba_report", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        summary_script = os.path.join(ROOT_DIR, "modules/arriba/bin/summarize_arriba_fusions.py"),
        formats = _formats,   # ③ 转发
    conda:
        "../arriba.yaml"
    container:
        sif("../arriba.yaml", "arriba_report")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("arriba_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
            sample_outdir = os.path.dirname(output.report)
            script = receipt_script(f"arriba_report_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = f"{sample_outdir}/arriba_report.sh"
            cmd = [
                "python", params.summary_script,
                "-p", ",".join(input.passed_fusions),
                "-d", ",".join(input.discarded_fusions),
                "-o", outdir + f"/{wildcards.genome}/arriba_report"
            ]
            # script's -f/--format is nargs="+" — single occurrence, all values
            cmd += ["--format"] + [str(fmt) for fmt in params.formats]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Arriba report generation for genome {wildcards.genome} successfully completed!"\n')
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
        passed_fusion_tsv = outdir + "/{genome}/{sample_id}/{sample_id}_passed_fusions.tsv",
        discarded_fusion_tsv = outdir + "/{genome}/{sample_id}/{sample_id}_discarded_fusions.tsv"
