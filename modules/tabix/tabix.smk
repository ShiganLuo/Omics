include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule tabix_bgzip:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf"
    output:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
    log:
        logdir + "/{sample_id}/bgzip.log"
    benchmark: logdir + "/{sample_id}/bgzip.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "tabix_bgzip", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 8000),
    params:
        bgzip = config.get("Procedure", {}).get("bgzip") or "bgzip"
    conda:
        "tabix.yaml"
    container:
        sif("tabix.yaml", "tabix_bgzip")
    run:
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        logger.info(f"Start bgzip for sample {wildcards.sample_id} at {current_time}")
        script = receipt_script(f"bgzip_{wildcards.sample_id}_{current_time}.sh")
        os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
        plain_script = os.path.join(outdir, f"bgzip.sh")
        cmd1 = [
            params.bgzip, input.vcf,
            "-o", output.vcf_gz
        ]
        cmd2 = [
            params.tabix, "-p", "vcf", output.vcf_gz
        ]
        with open(script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(" ".join(cmd1) + "\n")
            f.write(" ".join(cmd2) + "\n")
            f.write(f'echo "script completed successfully"\n')
        shutil.copyfile(script, plain_script)
        if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
            script_only_touch(output)
        else:
            shell("bash {plain_script} >> {log} 2>&1")
rule tabix_result:
    input:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
