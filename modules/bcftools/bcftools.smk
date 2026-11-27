include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule bcftools_sort:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf"
    output:
        vcf = outdir + "/sort/{sample_id}/{sample_id}.sorted.vcf.gz"
    log:
        logdir + "/{sample_id}/bcftools_sort.log"
    benchmark: logdir + "/{sample_id}/bcftools_sort.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "bcftools_sort", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "bcftools.yaml"
    container:
        sif("bcftools.yaml", "bcftools_sort")
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("bcftools_sort", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start bcftools sort for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.vcf))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"bcftools_sort_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"bcftools_sort.sh")
            cmd = [
                params.bcftools, "sort",
                "-O", "z",
                "-o", str(output.vcf),
                str(input.vcf),
            ]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'bcftools sort for {wildcards.sample_id} was completed successfully at {current_time}'")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"bcftools sort failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"bcftools sort failed for sample {wildcards.sample_id}: {e}")

rule bcftools_index_csi:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf.gz"
    output:
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
    log:
        logdir + "/{sample_id}/bcftools_index.log"
    benchmark: logdir + "/{sample_id}/bcftools_index.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "bcftools_index_csi", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "bcftools.yaml"
    container:
        sif("bcftools.yaml", "bcftools_index_csi")
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("bcftools_index", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start bcftools index for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.csi))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"bcftools_index_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"bcftools_index.sh")
            cmd = [
                params.bcftools, "index",
                str(input.vcf),
            ]
            with open(script, 'w') as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f"echo 'bcftools index for {wildcards.sample_id} was completed successfully at {current_time}'")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"bcftools index failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"bcftools index failed for sample {wildcards.sample_id}: {e}")

rule bcftools_result:
    input:
        vcf = outdir + "/sort/{sample_id}/{sample_id}.sorted.vcf.gz",
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
