include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/diff")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
# 参数精确到 rule:formats 按 rule 独立配置(Params.fiberseq_diff.<rule_name>.formats,
# 无 module 级共享回退),每 rule 一个 parse 期变量(config 全量路径读取 + str 规范化)。
_formats_run = config.get("Params", {}).get("fiberseq_diff", {}).get("fiberseq_diff_run", {}).get("formats") or ["pdf"]
if isinstance(_formats_run, str):
    _formats_run = [_formats_run]

# upstream contract: matrices from fiberseq_analysis/quant, annotations from
# fiberseq_analysis/annotate
QUANT_DIR = config.get("quant_dir") or indir + "/quant"
ANNOT_DIR = config.get("annot_dir") or indir + "/annotate"
# <- Params.fiberseq_diff.fiberseq_diff_run.annot_sample(input 指令无法调用
# 函数,每 rule 一个 parse 期变量;默认取首个 sample)
ANNOT_SAMPLE = (config.get("Params", {}).get("fiberseq_diff", {}).get("fiberseq_diff_run", {}).get("annot_sample", None)
                or (samples[0] if samples else ""))

rule fiberseq_diff_run:
    """Per-region differential accessibility for one contrast."""
    input:
        matrix = QUANT_DIR + "/matrix/percent_accessible_matrix.tsv",
        annot = ANNOT_DIR + "/" + ANNOT_SAMPLE + "/" + ANNOT_SAMPLE + ".peaks_annotated.tsv",
    output:
        diff = outdir + "/{contrast}/diff_{contrast}.tsv",
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        plots = [outdir + "/{contrast}/volcano_{contrast}." + fmt
                 for fmt in _formats_run],
        class_summary = outdir + "/{contrast}/class_summary_{contrast}.tsv",
        gene_table = outdir + "/{contrast}/TEcount/{contrast}.TEcount_Gene.name.tsv",
    log:
        logdir + "/group/fiberseq_diff_run_{contrast}.log"
    threads: rule_res("fiberseq_diff_run").get("threads", 4)
    benchmark: logdir + "/group/fiberseq_diff_run_{contrast}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_diff_run", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_diff_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/diff_accessibility.py",
        alpha = config.get("Params", {}).get("fiberseq_diff", {}).get("fiberseq_diff_run", {}).get("alpha", 0.05),
        group_spec = config.get("Params", {}).get("fiberseq_diff", {}).get("fiberseq_diff_run", {}).get("group_spec", "") or "",
        formats = _formats_run,   # ③ 转发
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_diff_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start diff analysis for {wildcards.contrast} at {current_time}")
            contrast_outdir = os.path.dirname(str(output.diff))
            os.makedirs(contrast_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_diff_run_{wildcards.contrast}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(contrast_outdir, f"fiberseq_diff_run_{wildcards.contrast}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.matrix),
                "--group-spec", params.group_spec,
                "-a", str(input.annot),
                "-C", str(wildcards.contrast),
                "-o", contrast_outdir,
                "-p", str(params.alpha),
                "-l", log_path,
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_diff_run for {wildcards.contrast} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_diff_run for {wildcards.contrast}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_diff_run for {wildcards.contrast}: {e}\n")

rule fiberseq_diff_result:
    input:
        matrix = QUANT_DIR + "/matrix/percent_accessible_matrix.tsv",
