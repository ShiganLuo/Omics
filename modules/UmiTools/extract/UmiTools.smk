include: "../../common/common.smk"

from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")
rule UmiTools_extract_single:
    input:
        fastq = indir + "/{sample_id}.single.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}.umi.single.fq.gz")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        extract_method = config.get('Procedure',{}).get('extract_method') or 'string',
        bc_pattern = config.get('Params',{}).get('umi_tools',{}).get('bc_pattern') or 'NNNXXXXNN'
    log:
        logdir + "/{sample_id}/umi_tools_extract_single_run.txt"
    conda:
        "../UmiTools.yaml"
    threads: 2
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("UmiTools_extract_single", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start UmiTools_extract_single for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"UmiTools_extract_single_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.umi_tools} extract \\\n")
                f.write(f"    --extract-method={params.extract_method} \\\n")
                f.write(f"    --bc-pattern={params.bc_pattern} \\\n")
                f.write(f"    -I {input.fastq} \\\n")
                f.write(f"    -S {output.fastq} \\\n")
                f.write(f"    > {log_path} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"UmiTools_extract_single for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"UmiTools_extract_single failed for sample {wildcards.sample_id}: {e}\n")
            raise e

rule UmiTools_extract_paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        fastq1 = temp(outdir + "/{sample_id}_1.umi.fq.gz"),
        fastq2 = temp(outdir + "/{sample_id}_2.umi.fq.gz")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        extract_method = config.get('Procedure',{}).get('extract_method') or 'string',
        bc_pattern = config.get('Params',{}).get('umi_tools',{}).get('bc_pattern') or 'NNNXXXXNN',
        bc_pattern2 = config.get('Params',{}).get('umi_tools',{}).get('bc_pattern2') or 'NNNXXXXNN'
    log:
        logdir + "/{sample_id}/umi_tools_extract_paired_run.txt"
    threads: 2
    conda:
        "../UmiTools.yaml"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("UmiTools_extract_paired", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start UmiTools_extract_paired for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"UmiTools_extract_paired_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.umi_tools} extract \\\n")
                f.write(f"    --extract-method={params.extract_method} \\\n")
                f.write(f"    --bc-pattern={params.bc_pattern} \\\n")
                f.write(f"    -S {output.fastq1} \\\n")
                f.write(f"    --bc-pattern2={params.bc_pattern2} \\\n")
                f.write(f"    --read2-in={input.fastq2} \\\n")
                f.write(f"    --read2-out={output.fastq2} \\\n")
                f.write(f"    > {log_path} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"UmiTools_extract_paired for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"UmiTools_extract_paired failed for sample {wildcards.sample_id}: {e}\n")
            raise e
