include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/footprint")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
# 参数精确到 rule:formats 按 rule 独立配置(Params.fiberseq_footprint.<rule_name>.formats,
# 无 module 级共享回退),每 rule 一个 parse 期变量(config 全量路径读取 + str 规范化)。
_formats_run = config.get("Params", {}).get("fiberseq_footprint", {}).get("fiberseq_footprint_run", {}).get("formats") or ["pdf"]
if isinstance(_formats_run, str):
    _formats_run = [_formats_run]
_formats_result = config.get("Params", {}).get("fiberseq_footprint", {}).get("fiberseq_footprint_result", {}).get("formats") or ["pdf"]
if isinstance(_formats_result, str):
    _formats_result = [_formats_result]

# upstream contract: motif sites come from fiberseq_analysis/motif
MOTIF_DIR = config.get("motif_dir") or indir + "/motif"

rule fiberseq_footprint_run:
    """Single-fiber TF footprinting (3-category) around motif sites."""
    input:
        sites = MOTIF_DIR + "/{sample_id}/{sample_id}.motif_sites.bed",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        categories = outdir + "/{sample_id}/{sample_id}.footprint_categories.tsv",
        profiles = outdir + "/{sample_id}/{sample_id}.footprint_profiles.tsv",
        scores = outdir + "/{sample_id}/{sample_id}.footprint_scores.tsv",
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        plots = [outdir + "/{sample_id}/{sample_id}.footprint_profile." + fmt
                 for fmt in _formats_run],
    log:
        logdir + "/{sample_id}/fiberseq_footprint_run.log"
    threads: rule_res("fiberseq_footprint").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_footprint_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_footprint", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_footprint_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/footprint_profile.py",
        flank = config.get("Params", {}).get("fiberseq_footprint", {}).get("fiberseq_footprint_run", {}).get("flank", 1000),
        formats = _formats_run,   # ③ 转发
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_footprint_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start footprint for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.categories))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_footprint_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_footprint_run_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.sites),
                "-a", str(input.m6a),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-f", str(params.flank),
                "-l", log_path,
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_footprint_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_footprint_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_footprint_run for {wildcards.sample_id}: {e}\n")

rule fiberseq_footprint_result:
    input:
        categories = expand(outdir + "/{sample_id}/{sample_id}.footprint_categories.tsv",
                            sample_id=samples),
        scores = expand(outdir + "/{sample_id}/{sample_id}.footprint_scores.tsv",
                        sample_id=samples),
        plots = expand(outdir + "/{sample_id}/{sample_id}.footprint_profile.{fmt}",
                       sample_id=samples, fmt=_formats_result),
