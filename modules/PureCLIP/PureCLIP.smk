include: "../common/common.smk"

outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")
rule pureclip:
    input:
        bam = indir + "/{sample_id}.dedup.bam",
        bai = indir + "/{sample_id}.dedup.bam.bai",
        fasta = config.get('genome',{}).get('fasta')
    output:
        sites = outdir + "/{sample_id}.pureclip.sites.bed",
        region = outdir + "/{sample_id}.pureclip.region.bed"
    params:
        pureclip = config.get('Procedure',{}).get('pureclip') or 'pureclip',
        ld = '--ld' if config.get('Params',{}).get('pureclip',{}).get('ld', True) else '',
        outdir = outdir
    threads: rule_res("pureclip").get("threads", 8)
    benchmark: logdir + "/{sample_id}/pureclip_run.benchmark.txt",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "pureclip", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "PureCLIP.yaml"
    container:
        sif("PureCLIP.yaml", "pureclip")
    log:
        log = logdir + "/{sample_id}/pureclip_run.txt"
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("pureclip", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start pureclip for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.sites))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"pureclip_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"pureclip_{wildcards.sample_id}.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.pureclip} {params.ld} -nt {threads} \\\n")
                f.write(f"    -i {input.bam} \\\n")
                f.write(f"    -bai {input.bai} \\\n")
                f.write(f"    -g {input.fasta} \\\n")
                f.write(f"    -o {output.sites} \\\n")
                f.write(f"    -or {output.region} \\\n")
                f.write(f"    >> {log} 2>&1\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during pureclip for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")
