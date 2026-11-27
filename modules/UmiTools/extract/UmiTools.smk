include: "../../common/common.smk"

outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")


rule UmiTools_extract_single:
    input:
        fastq = indir + "/{sample_id}.single.fq.gz"
    output:
        fastq = maybe_temp(outdir + "/{sample_id}.umi.single.fq.gz", "UmiTools")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        extract_method = config.get('Procedure',{}).get('extract_method') or 'string',
        bc_pattern = config.get("Params", {}).get("umi_tools", {}).get("UmiTools_extract_single", {}).get("bc_pattern") or "NNNXXXXNN"
    log:
        logdir + "/{sample_id}/umi_tools_extract_single_run.txt"
    conda:
        "../UmiTools.yaml"
    container:
        sif("../UmiTools.yaml", "UmiTools_extract_single")
    threads: rule_res("UmiTools_extract_single").get("threads", 2)
    benchmark: logdir + "/{sample_id}/umi_tools_extract_single_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "UmiTools_extract_single", input_size=input_size_named(input, 'fastq'),
            attempt=attempt).get("mem_mb", 8000),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("UmiTools_extract_single", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start UmiTools_extract_single for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"UmiTools_extract_single_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"UmiTools_extract_single.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.umi_tools} extract \\\n")
                f.write(f"    --extract-method={params.extract_method} \\\n")
                f.write(f"    --bc-pattern={params.bc_pattern} \\\n")
                f.write(f"    -I {input.fastq} \\\n")
                f.write(f"    -S {output.fastq} \\\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"UmiTools_extract_single for sample {wildcards.sample_id} completed")
        except Exception as e:
            rule_logger.error(f"Error occurred: {e}")
            raise RuntimeError(f"Error occurred: {e}, you can check the log file {log_path} for more details.")

rule UmiTools_extract_paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        fastq1 = maybe_temp(outdir + "/{sample_id}_1.umi.fq.gz", "UmiTools"),
        fastq2 = maybe_temp(outdir + "/{sample_id}_2.umi.fq.gz", "UmiTools")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        extract_method = config.get('Procedure',{}).get('extract_method') or 'string',
        bc_pattern = config.get("Params", {}).get("umi_tools", {}).get("UmiTools_extract_paired", {}).get("bc_pattern") or "NNNXXXXNN",
        bc_pattern2 = config.get("Params", {}).get("umi_tools", {}).get("UmiTools_extract_paired", {}).get("bc_pattern2") or "NNNXXXXNN"
    log:
        logdir + "/{sample_id}/umi_tools_extract_paired_run.txt"
    threads: rule_res("UmiTools_extract_paired").get("threads", 2)
    benchmark: logdir + "/{sample_id}/umi_tools_extract_paired_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "UmiTools_extract_paired", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../UmiTools.yaml"
    container:
        sif("../UmiTools.yaml", "UmiTools_extract_paired")
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("UmiTools_extract_paired", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start UmiTools_extract_paired for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"UmiTools_extract_paired_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"UmiTools_extract_paired.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.umi_tools} extract \\\n")
                f.write(f"    --extract-method={params.extract_method} \\\n")
                f.write(f"    --bc-pattern={params.bc_pattern} \\\n")
                f.write(f"    -S {output.fastq1} \\\n")
                f.write(f"    --bc-pattern2={params.bc_pattern2} \\\n")
                f.write(f"    --read2-in={input.fastq2} \\\n")
                f.write(f"    --read2-out={output.fastq2} \\\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"UmiTools_extract_paired for sample {wildcards.sample_id} completed")
        except Exception as e:
            rule_logger.error(f"Error occurred: {e}")
            raise RuntimeError(f"Error occurred: {e}, you can check the log file {log_path} for more details.")
