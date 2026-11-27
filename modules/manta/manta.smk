include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
bam_suffix = config.get("bam_suffix", "")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

rule manta_config_germline:
    """Configure Manta SV caller."""
    input:
        bam = indir + "/{sample_id}/{sample_id}" + bam_suffix + ".bam",
        fasta = config.get("genome", {}).get("fasta")
    output:
        config = outdir + "/{sample_id}/runWorkflow.py.config.pickle",
        script = outdir + "/{sample_id}/runWorkflow.py"
    log:
        logdir + "/{sample_id}/manta_config.log"
    benchmark: logdir + "/{sample_id}/manta_config.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "manta_config", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    params:
        manta_outdir = outdir + "/{sample_id}"
    conda:
        "manta.yaml"
    container:
        sif("manta.yaml", "manta_config")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger(logger_name="manta_config", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start manta config for sample {wildcards.sample_id} at {current_time}")
            script = receipt_script(f"manta_config_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/manta_config.sh")
            cmd = [
                "configManta.py",
                "--bam", str(input.bam),
                "--referenceFasta", str(input.fasta),
                "--runDir", str(params.manta_outdir)
            ]
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"manta_config failed for sample {wildcards.sample_id} with error: {e}")
            raise RuntimeError(f"manta_config failed for sample {wildcards.sample_id} with error: {e}")

rule manta_run_germline:
    """Run Manta SV calling."""
    input:
        script = outdir + "/{sample_id}/runWorkflow.py",
        config = outdir + "/{sample_id}/runWorkflow.py.config.pickle"
    output:
        vcf = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz",
        tbi = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz.tbi"
    log:
        logdir + "/{sample_id}/manta_run.log"
    threads:
        config.get("Params", {}).get("manta", {}).get("manta_run_germline", {}).get("threads", 8)
    benchmark: logdir + "/{sample_id}/manta_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "manta_run_germline", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "manta.yaml"
    container:
        sif("manta.yaml", "manta_run")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger(logger_name="manta_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start manta run for sample {wildcards.sample_id} at {current_time}")
            script = receipt_script(f"manta_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/manta_run.sh")
            cmd = [
                "python", str(input.script),
                "-j", str(threads)
            ]
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"manta_run failed for sample {wildcards.sample_id} with error: {e}")
            raise RuntimeError(f"manta_run failed for sample {wildcards.sample_id} with error: {e}")

rule manta_config_somatic:
    """Configure Manta SV caller for tumor-normal (somatic) mode."""
    input:
        tumor_bam = indir + "/{experimental_sample_id}/{experimental_sample_id}" + bam_suffix + ".bam",
        normal_bam = indir + "/{normal_sample_id}/{normal_sample_id}" + bam_suffix + ".bam",
        fasta = config.get("genome", {}).get("fasta")
    output:
        config = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/runWorkflow.py.config.pickle",
        script = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/runWorkflow.py"
    log:
        logdir + "/{normal_sample_id}_vs_{experimental_sample_id}/manta_config_somatic.log"
    benchmark: logdir + "/{normal_sample_id}_vs_{experimental_sample_id}/manta_config_somatic.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "manta_config_somatic", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        manta_outdir = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}"
    conda:
        "manta.yaml"
    container:
        sif("manta.yaml", "manta_config_somatic")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger(logger_name="manta_config_somatic", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start manta somatic config for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} at {current_time}")
            script = receipt_script(f"manta_config_somatic_{wildcards.experimental_sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, f"{wildcards.normal_sample_id}_vs_{wildcards.experimental_sample_id}/manta_config_somatic.sh")
            cmd = [
                "configManta.py",
                "--tumorBam", str(input.tumor_bam),
                "--normalBam", str(input.normal_bam),
                "--referenceFasta", str(input.fasta),
                "--runDir", str(params.manta_outdir)
            ]
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"manta_config_somatic failed for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} with error: {e}")
            raise RuntimeError(f"manta_config_somatic failed for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} with error: {e}")


rule manta_run_somatic:
    """Run Manta SV calling for tumor-normal (somatic) mode."""
    input:
        script = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/runWorkflow.py",
        config = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/runWorkflow.py.config.pickle"
    output:
        vcf = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/results/variants/somaticSV.vcf.gz",
        tbi = outdir + "/{normal_sample_id}_vs_{experimental_sample_id}/results/variants/somaticSV.vcf.gz.tbi"
    log:
        logdir + "/manta/{normal_sample_id}_vs_{experimental_sample_id}/manta_run_somatic.log"
    threads:
        config.get("Params", {}).get("manta", {}).get("manta_run_somatic", {}).get("threads", 8)
    benchmark: logdir + "/manta/{normal_sample_id}_vs_{experimental_sample_id}/manta_run_somatic.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "manta_run_somatic", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "manta.yaml"
    container:
        sif("manta.yaml", "manta_run_somatic")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger(logger_name="manta_run_somatic", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start manta somatic run for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} at {current_time}")
            script = receipt_script(f"manta_run_somatic_{wildcards.experimental_sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, f"{wildcards.normal_sample_id}_vs_{wildcards.experimental_sample_id}/manta_run_somatic.sh")
            cmd = [
                "python", str(input.script),
                "-j", str(threads)
            ]
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"manta_run_somatic failed for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} with error: {e}")
            raise RuntimeError(f"manta_run_somatic failed for normal sample {wildcards.normal_sample_id} and experimental sample {wildcards.experimental_sample_id} with error: {e}")

rule manta_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        vcf = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz",
        tbi = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz.tbi"
