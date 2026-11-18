include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/coactuation")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
rule fiberseq_coactuation_run:
    """Same-fiber co-actuation of region pairs (Fisher + actual vs expected)."""
    input:
        peaks = indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        fire = indir + "/{sample_id}/{sample_id}.fire.bed.gz",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        pairs = outdir + "/{sample_id}/{sample_id}.coactuation.tsv",
        summary = outdir + "/{sample_id}/{sample_id}.coactuation_summary.tsv",
    log:
        logdir + "/{sample_id}/fiberseq_coactuation_run.log"
    threads: rule_res("fiberseq_coactuation").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_coactuation_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_coactuation", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_coactuation_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/coactuation.py",
        # 参数精确到 rule:Params.fiberseq_coactuation.fiberseq_coactuation_run.<key>(config 全量路径读取)
        max_distance = config.get("Params", {}).get("fiberseq_coactuation", {}).get("fiberseq_coactuation_run", {}).get("max_distance", 50000),
        min_fibers = config.get("Params", {}).get("fiberseq_coactuation", {}).get("fiberseq_coactuation_run", {}).get("min_fibers", 10),
        max_pairs = config.get("Params", {}).get("fiberseq_coactuation", {}).get("fiberseq_coactuation_run", {}).get("max_pairs", 200000),
        alpha = config.get("Params", {}).get("fiberseq_coactuation", {}).get("fiberseq_coactuation_run", {}).get("alpha", 0.05),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_coactuation_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start co-actuation for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.pairs))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_coactuation_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_coactuation_run_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-p", str(input.peaks),
                "-f", str(input.fire),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-d", str(params.max_distance),
                "-i", str(params.min_fibers),
                "-m", str(params.max_pairs),
                "-a", str(params.alpha),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_coactuation_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_coactuation_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_coactuation_run for {wildcards.sample_id}: {e}\n")

rule fiberseq_coactuation_result:
    input:
        pairs = expand(outdir + "/{sample_id}/{sample_id}.coactuation.tsv",
                       sample_id=samples),
