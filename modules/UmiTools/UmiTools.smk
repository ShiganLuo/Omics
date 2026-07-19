include: "../common/common.smk"

from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")

rule umi_tools_dedup_for_hisat2:
    input:
        bam = indir + "/{sample_id}.bam"
    output:
        bam = temp(outdir + "/{sample_id}.dedup.bam"),
        bai = temp(outdir + "/{sample_id}.dedup.bam.bai"),
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        method = config.get('Params',{}).get('umi_tools',{}).get('method') or 'unique',
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    threads: 2
    conda:
        "UmiTools.yaml"
    log:
        log = logdir + "/{sample_id}/umi_tools_dedup_run.txt"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("umi_tools_dedup_for_hisat2", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start umi_tools_dedup_for_hisat2 for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"umi_tools_dedup_for_hisat2_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.umi_tools} dedup \\\n")
                f.write(f"    --method={params.method} \\\n")
                f.write(f"    -I {input.bam} \\\n")
                f.write(f"    -S {output.bam} \\\n")
                f.write(f"    > {log_path} 2>&1\n")
                f.write(f"{params.samtools} index {output.bam}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"umi_tools_dedup_for_hisat2 for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"umi_tools_dedup_for_hisat2 failed for sample {wildcards.sample_id}: {e}\n")
            raise e

rule umi_tools_dedup_for_star:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = temp(outdir + "/{sample_id}.dedup.bam"),
        bai = temp(outdir + "/{sample_id}.dedup.bam.bai")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        method = config.get('Params',{}).get('umi_tools',{}).get('method') or 'unique',
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    threads: 2
    conda:
        "UmiTools.yaml"
    log:
        log = logdir + "/{sample_id}/umi_tools_dedup_run.txt"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("umi_tools_dedup_for_star", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start umi_tools_dedup_for_star for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"umi_tools_dedup_for_star_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.umi_tools} dedup --method={params.method} \\\n")
                f.write(f"    -I {input.bam} -S {output.bam} > {log_path} 2>&1\n")
                f.write(f"{params.samtools} index {output.bam} 2>>{log_path}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"umi_tools_dedup_for_star for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"umi_tools_dedup_for_star failed for sample {wildcards.sample_id}: {e}\n")
            raise e
