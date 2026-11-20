include: "../common/common.smk"

indir = config.get("indir","input")
outdir = config.get("outdir", "output")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule bam_flagstat:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"
    log:
        logdir + "/{sample_id}/flagstat.log"
    benchmark: logdir + "/{sample_id}/flagstat.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "bam_flagstat", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "samtools.yaml"
    container:
        sif("samtools.yaml", "bam_flagstat")
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("bam_flagstat", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start samtools flagstat for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.flagstat))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"bam_flagstat_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"bam_flagstat.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.samtools} flagstat {input.bam} > {output.flagstat}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during samtools flagstat for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule samtools_result:
    input:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"
