include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

rule manta_config:
    """Configure Manta SV caller."""
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        fasta = config.get("genome", {}).get("fasta")
    output:
        config = outdir + "/{sample_id}/runWorkflow.py.config.pickle",
        script = outdir + "/{sample_id}/runWorkflow.py"
    log:
        logdir + "/{sample_id}/manta_config.log"
    params:
        manta_outdir = outdir + "/{sample_id}"
    conda:
        "manta.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="manta_config", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start manta config for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/manta_config_{current_time}.sh")
            cmd = [
                "configManta.py",
                "--bam", str(input.bam),
                "--referenceFasta", str(input.fasta),
                "--runDir", str(params.manta_outdir)
            ]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")

rule manta_run:
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
        config.get("Params", {}).get("manta", {}).get("threads") or 8
    conda:
        "manta.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger(logger_name="manta_run", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start manta run for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/manta_run_{current_time}.sh")
            cmd = [
                "python", str(input.script),
                "-j", str(threads)
            ]
            with open(script, 'w') as f:
                f.write(' '.join(cmd) + '\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, 'a') as f:
                f.write(f"Error: {e}\n")
            raise f"Error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Completed at {current_time}")

rule manta_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        vcf = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz",
        tbi = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz.tbi"
