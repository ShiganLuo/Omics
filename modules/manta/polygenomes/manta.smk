include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
# 参数精确到 rule:threads 按 rule 独立配置(Params.manta.<rule_name>.threads),
# 无 module 级共享。

rule manta_config:
    """Configure Manta SV caller."""
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"]
    output:
        config = outdir + "/{genome}/{sample_id}/runWorkflow.py.config.pickle",
        script = outdir + "/{genome}/{sample_id}/runWorkflow.py"
    log:
        logdir + "/{genome}/{sample_id}/manta_config.log"
    benchmark: logdir + "/{genome}/{sample_id}/manta_config.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "manta_config", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    params:
        manta_outdir = lambda wildcards: outdir + f"/{wildcards.genome}/{wildcards.sample_id}"
    conda:
        "../manta.yaml"
    container:
        sif("../manta.yaml", "manta_config")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        logger = setup_logger(logger_name="manta_config", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start manta config for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"manta_config_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.sample_id}/manta_config.sh")
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
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise RuntimeError(f"Error: {e}")
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")

rule manta_run:
    """Run Manta SV calling."""
    input:
        script = outdir + "/{genome}/{sample_id}/runWorkflow.py",
        config = outdir + "/{genome}/{sample_id}/runWorkflow.py.config.pickle"
    output:
        vcf = outdir + "/{genome}/{sample_id}/results/variants/candidateSV.vcf.gz",
        tbi = outdir + "/{genome}/{sample_id}/results/variants/candidateSV.vcf.gz.tbi"
    log:
        logdir + "/{genome}/{sample_id}/manta_run.log"
    threads:
        config.get("Params", {}).get("manta", {}).get("manta_run", {}).get("threads", 8)
    benchmark: logdir + "/{genome}/{sample_id}/manta_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "manta_run", input_size=input_size_named(input, 'vcf'),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../manta.yaml"
    container:
        sif("../manta.yaml", "manta_run")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        logger = setup_logger(logger_name="manta_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start manta run for sample {wildcards.sample_id} genome {wildcards.genome} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"manta_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"{wildcards.genome}/{wildcards.sample_id}/manta_run.sh")
            cmd = [
                "python", str(input.script),
                "-j", str(threads)
            ]
            with open(script, 'w') as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(' '.join(cmd) + '\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise RuntimeError(f"Error: {e}")
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")

rule manta_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        vcf = outdir + "/{genome}/{sample_id}/results/variants/candidateSV.vcf.gz",
        tbi = outdir + "/{genome}/{sample_id}/results/variants/candidateSV.vcf.gz.tbi"
