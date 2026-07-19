include: "../common/common.smk"

from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")
# neet for test
rule trimmomatic_Paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        fastq1 = temp(outdir + "/{sample_id}/{sample_id}_1.fq.gz"),
        fastq2 = temp(outdir + "/{sample_id}/{sample_id}_2.fq.gz"),
        report = outdir + "/{sample_id}/trimmomatic_report.txt"
    params:
        outdir = outdir,
        trimmomatic = config.get('Procedure',{}).get('trimmomatic') or 'trimmomatic',
        adapter = config.get('Params',{}).get("trimmomatic", {}).get('adapter_pe'),
        cmd = lambda wildcards: (
            f"java -jar {config.get('Procedure',{}).get('trimmomatic')}"
            if str(config.get('Procedure',{}).get('trimmomatic')).endswith(".jar")
            else config.get('Procedure',{}).get('trimmomatic') or "trimmomatic"
        )
    threads: 6
    conda:
        "trimmomatic.yaml" if not str(config.get('Procedure',{}).get('trimmomatic')).endswith(".jar") else None
    log:
        log = logdir + "/{sample_id}/trimmomatic.txt"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("trimmomatic_Paired", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start trimmomatic_Paired for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(params.outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"trimmomatic_Paired_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.cmd} PE -threads {threads} \\\n")
                f.write(f"    {input.fastq1} {input.fastq2} \\\n")
                f.write(f"    -summary {output.report} \\\n")
                f.write(f"    {output.fastq1} {params.outdir}/{wildcards.sample_id}/{wildcards.sample_id}_1.unpaired.fq.gz \\\n")
                f.write(f"    {output.fastq2} {params.outdir}/{wildcards.sample_id}/{wildcards.sample_id}_2.unpaired.fq.gz \\\n")
                f.write(f"    ILLUMINACLIP:{params.adapter}:2:30:10 \\\n")
                f.write(f"    LEADING:3 TRAILING:3 \\\n")
                f.write(f"    SLIDINGWINDOW:4:15 \\\n")
                f.write(f"    MINLEN:80 \\\n")
                f.write(f"    2> {log_path}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"trimmomatic_Paired for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"trimmomatic_Paired failed for sample {wildcards.sample_id}: {e}\n")
            raise e

rule trimmomatic_Single:
    input:
        fastq = indir + "/{sample_id}.single.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}/{sample_id}.single.fq.gz"),
        report = outdir + "/{sample_id}/trimmomatic_report.txt"
    params:
        outdir = outdir,
        trimmomatic = config.get('Procedure',{}).get('trimmomatic') or 'trimmomatic',
        adapter = config.get('Params',{}).get("trimmomatic", {}).get('adapter_se'),
        cmd = lambda wildcards: (
            f"java -jar {config.get('Procedure',{}).get('trimmomatic')}"
            if str(config.get('Procedure',{}).get('trimmomatic')).endswith(".jar")
            else config.get('Procedure',{}).get('trimmomatic') or "trimmomatic"
        )
    threads: 6
    conda:
        "trimmomatic.yaml"
    log:
        log = logdir + "/{sample_id}/trimmomatic.txt"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("trimmomatic_Single", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start trimmomatic_Single for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(params.outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"trimmomatic_Single_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.cmd} SE -threads {threads} \\\n")
                f.write(f"    {input.fastq} {output.fastq} \\\n")
                f.write(f"    -summary {output.report} \\\n")
                f.write(f"    ILLUMINACLIP:{params.adapter}:2:30:10 \\\n")
                f.write(f"    LEADING:3 TRAILING:3 \\\n")
                f.write(f"    SLIDINGWINDOW:4:15 \\\n")
                f.write(f"    MINLEN:80 \\\n")
                f.write(f"    2> {log_path}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"trimmomatic_Single for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"trimmomatic_Single failed for sample {wildcards.sample_id}: {e}\n")
            raise e
