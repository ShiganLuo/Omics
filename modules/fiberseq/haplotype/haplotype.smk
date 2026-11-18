include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/haplotype")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})

# upstream contract: per-hap quant tables from fiberseq_analysis/quant
QUANT_DIR = config.get("quant_dir") or indir + "/quant"

rule fiberseq_haplotype_diff:
    """Haplotype-specific accessibility differences (Fisher per region)."""
    input:
        hap1 = QUANT_DIR + "/{sample_id}/{sample_id}.hap1.region_quant.tsv",
        hap2 = QUANT_DIR + "/{sample_id}/{sample_id}.hap2.region_quant.tsv",
    output:
        diff = outdir + "/{sample_id}/{sample_id}.haplotype_diff.tsv",
        summary = outdir + "/{sample_id}/{sample_id}.haplotype_summary.tsv",
    log:
        logdir + "/{sample_id}/fiberseq_haplotype_diff.log"
    threads: rule_res("fiberseq_haplotype_diff").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_haplotype_diff.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_haplotype_diff", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_haplotype_diff")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/haplotype_diff.py",
        # 参数精确到 rule:Params.fiberseq_haplotype.fiberseq_haplotype_diff.<key>(config 全量路径读取)
        min_fibers = config.get("Params", {}).get("fiberseq_haplotype", {}).get("fiberseq_haplotype_diff", {}).get("min_fibers", 10),
        alpha = config.get("Params", {}).get("fiberseq_haplotype", {}).get("fiberseq_haplotype_diff", {}).get("alpha", 0.05),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_haplotype_diff", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start haplotype diff for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.diff))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_haplotype_diff_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_haplotype_diff_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-1", str(input.hap1),
                "-2", str(input.hap2),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-i", str(params.min_fibers),
                "-a", str(params.alpha),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_haplotype_diff for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_haplotype_diff for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_haplotype_diff for {wildcards.sample_id}: {e}\n")

rule fiberseq_haplotype_result:
    input:
        diff = expand(outdir + "/{sample_id}/{sample_id}.haplotype_diff.tsv",
                      sample_id=samples),
