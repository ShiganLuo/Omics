include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/nucpos")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
# 参数精确到 rule:formats 按 rule 独立配置(Params.fiberseq_nucpos.<rule_name>.formats,
# 无 module 级共享回退),每 rule 一个 parse 期变量(config 全量路径读取 + str 规范化)。
_formats_run = config.get("Params", {}).get("fiberseq_nucpos", {}).get("fiberseq_nucpos_run", {}).get("formats") or ["pdf"]
if isinstance(_formats_run, str):
    _formats_run = [_formats_run]
_formats_result = config.get("Params", {}).get("fiberseq_nucpos", {}).get("fiberseq_nucpos_result", {}).get("formats") or ["pdf"]
if isinstance(_formats_result, str):
    _formats_result = [_formats_result]

rule fiberseq_nucpos_run:
    """Nucleosome positioning/phasing around TSS anchors."""
    input:
        gtf = lambda wildcards: _require(genome_cfg.get("gtf"), "genome.gtf", "fiberseq_nucpos"),
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        per_anchor = outdir + "/{sample_id}/{sample_id}.nuc_offset_per_anchor.tsv",
        summary = outdir + "/{sample_id}/{sample_id}.nuc_offset_summary.tsv",
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        plots = [outdir + "/{sample_id}/{sample_id}.nucpos." + fmt
                 for fmt in _formats_run],
    log:
        logdir + "/{sample_id}/fiberseq_nucpos_run.log"
    threads: rule_res("fiberseq_nucpos").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_nucpos_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_nucpos", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 6000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_nucpos_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/nucpos.py",
        upstream = config.get("Params", {}).get("fiberseq_nucpos", {}).get("fiberseq_nucpos_run", {}).get("upstream_nucs", 5),
        downstream = config.get("Params", {}).get("fiberseq_nucpos", {}).get("fiberseq_nucpos_run", {}).get("downstream_nucs", 5),
        min_fibers = config.get("Params", {}).get("fiberseq_nucpos", {}).get("fiberseq_nucpos_run", {}).get("min_fibers", 10),
        formats = _formats_run,   # ③ 转发
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_nucpos_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start nucpos for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.per_anchor))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"fiberseq_nucpos_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_nucpos_run_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-g", str(input.gtf),
                "--anchor", "tss",
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-5", str(params.upstream),
                "-3", str(params.downstream),
                "-i", str(params.min_fibers),
                "-l", log_path,
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_nucpos_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_nucpos_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_nucpos_run for {wildcards.sample_id}: {e}\n")

rule fiberseq_nucpos_result:
    input:
        summary = expand(outdir + "/{sample_id}/{sample_id}.nuc_offset_summary.tsv",
                         sample_id=samples),
        plots = expand(outdir + "/{sample_id}/{sample_id}.nucpos.{fmt}",
                       sample_id=samples, fmt=_formats_result),
