include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "output/raw_fastq")

rule soapnuke_filter_paired:
    input:
        fastq1 = indir + "/{sample_id}/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}/{sample_id}_2.fq.gz"
    output:
        clean_fastq1 = outdir + "/{sample_id}/{sample_id}_1.fq.gz",
        clean_fastq2 = outdir + "/{sample_id}/{sample_id}_2.fq.gz"
    log:
        logdir + "/{sample_id}/soapnuke_filter_paired.log"
    params:
        min_len = 15,
        low_qual = 0.2,
        max_n = 0.05,
        qual_type = 2,
        SOAPnuke = config.get('Procedure', {}).get('SOAPnuke') or 'SOAPnuke',
        workdir = lambda wildcards: f"{outdir}/{wildcards.sample_id}",
        clean_fastq1 = lambda wildcards: f"{wildcards.sample_id}_1.fq.gz",
        clean_fastq2 = lambda wildcards: f"{wildcards.sample_id}_2.fq.gz",
    conda:
        "SOAPnuke.yaml"
    threads: 8
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("soapnuke_filter_paired", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start soapnuke_filter_paired for sample {wildcards.sample_id} at {current_time}")

            sample_outdir = params.workdir
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"soapnuke_filter_paired_{current_time}.sh")

            cmd = [
                params.SOAPnuke, "filter",
                "-l", str(params.min_len),
                "-q", str(params.low_qual),
                "-n", str(params.max_n),
                "-1", input.fastq1,
                "-2", input.fastq2,
                "-C", params.clean_fastq1,
                "-D", params.clean_fastq2,
                "-Q", str(params.qual_type),
                "-o", sample_outdir,
                "-T", str(threads),
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + "\n")
                f.write(f"mv {sample_outdir}/{params.clean_fastq1} {output.clean_fastq1}\n")
                f.write(f"mv {sample_outdir}/{params.clean_fastq2} {output.clean_fastq2}\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"soapnuke_filter_paired for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"soapnuke_filter_paired failed for sample {wildcards.sample_id}: {e}\n")
            raise e

rule soapnuke_filter_single:
    input:
        fastq = indir + "/{sample_id}/{sample_id}.single.fq.gz"
    output:
        clean_fastq = outdir + "/{sample_id}/{sample_id}.single.fq.gz"
    log:
        logdir + "/{sample_id}/soapnuke_filter_single.log"
    params:
        min_len = 15,
        low_qual = 0.2,
        max_n = 0.05,
        qual_type = 2,
        SOAPnuke = config.get('Procedure', {}).get('SOAPnuke') or 'SOAPnuke',
        workdir = lambda wildcards: f"{outdir}/{wildcards.sample_id}",
        clean_fastq = lambda wildcards: f"{wildcards.sample_id}.single.fq.gz"
    conda:
        "SOAPnuke.yaml"
    threads: 8
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("soapnuke_filter_single", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start soapnuke_filter_single for sample {wildcards.sample_id} at {current_time}")

            sample_outdir = params.workdir
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"soapnuke_filter_single_{current_time}.sh")

            cmd = [
                params.SOAPnuke, "filter",
                "-l", str(params.min_len),
                "-q", str(params.low_qual),
                "-n", str(params.max_n),
                "-Q", str(params.qual_type),
                "-1", input.fastq,
                "-C", params.clean_fastq,
                "-o", sample_outdir,
                "-T", str(threads),
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + "\n")
                f.write(f"mv {sample_outdir}/{params.clean_fastq} {output.clean_fastq}\n")
            shell(f"bash {script} >> {log_path} 2>&1")

            rule_logger.info(f"soapnuke_filter_single for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"soapnuke_filter_single failed for sample {wildcards.sample_id}: {e}\n")
            raise e

rule soapnuke_filter_paired_result:
    input:
        clean_fastq1 = outdir + "/{sample_id}/{sample_id}_1.fq.gz",
        clean_fastq2 = outdir + "/{sample_id}/{sample_id}_2.fq.gz",

rule soapnuke_filter_single_result:
    input:
        clean_fastq = outdir + "/{sample_id}/{sample_id}.single.fq.gz",
