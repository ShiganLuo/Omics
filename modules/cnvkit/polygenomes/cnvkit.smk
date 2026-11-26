include: "../../common/common.smk"


outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
control_samples = config.get("control_samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
# 参数精确到 rule:method 按 rule 独立配置(Params.cnvkit.<rule_name>.method),
# 无 module 级共享。

def get_cnvkit_bams(wildcards):
    """Get BAM files for CNVkit batch analysis."""
    treat_bams = [f"{indir}/{sid}/{sid}.bam" for sid in samples]
    control_bams = [f"{indir}/{sid}/{sid}.bam" for sid in control_samples] if control_samples else []
    return {"treat": treat_bams, "control": control_bams}

rule cnvkit_reference:
    """Build CNVkit reference from control samples."""
    input:
        control_bams = expand(indir + "/{sid}/{sid}.bam", sid=control_samples),
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"],
        access = lambda wildcards: config["genome"][wildcards.genome]["access"]
    output:
        ref = outdir + "/{genome}/reference/reference.cnn"
    log:
        logdir + "/{genome}/cnvkit_reference.log"
    threads: rule_res("cnvkit_reference").get("threads", 8)
    benchmark: logdir + "/{genome}/cnvkit_reference.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "cnvkit_reference", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "../cnvkit.yaml"
    container:
        sif("../cnvkit.yaml", "cnvkit_reference")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_reference", {}).get("method", "wgs"),
        control_bams = lambda wildcards: " ".join([f"{indir}/{sid}/{sid}.bam" for sid in control_samples])
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("cnvkit_reference", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit reference build for genome {wildcards.genome} at {current_time}")
            sample_outdir = os.path.dirname(str(output.ref))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"cnvkit_reference_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(sample_outdir, f"cnvkit_reference.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.cnvkit} batch \\\n")
                f.write(f"    {params.control_bams} \\\n")
                f.write(f"    --method {params.method} \\\n")
                f.write(f"    --fasta {input.fasta} \\\n")
                f.write(f"    --access {input.access} \\\n")
                f.write(f"    --processes {threads} \\\n")
                f.write(f"    --output-reference {output.ref} \\\n")
                f.write(f"    >> {log} 2>&1\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit reference build for genome {wildcards.genome}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule cnvkit_batch:
    """Run CNVkit batch analysis on all samples."""
    input:
        treat_bams = expand(indir + "/{sid}/{sid}.bam", sid=samples),
        ref = outdir + "/{genome}/reference/reference.cnn" if control_samples else [],
        fasta = lambda wildcards: config["genome"][wildcards.genome]["fasta"],
        access = lambda wildcards: config["genome"][wildcards.genome]["access"]
    output:
        cnr = expand(outdir + "/{{genome}}/cnv/{sample_id}.cnr", sample_id=samples),
        cns = expand(outdir + "/{{genome}}/cnv/{sample_id}.cns", sample_id=samples)
    log:
        logdir + "/{genome}/cnvkit_batch.log"
    threads: rule_res("cnvkit_batch").get("threads", 8)
    benchmark: logdir + "/{genome}/cnvkit_batch.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "cnvkit_batch", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../cnvkit.yaml"
    container:
        sif("../cnvkit.yaml", "cnvkit_batch")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_batch", {}).get("method", "wgs"),
        treat_bams = lambda wildcards: " ".join([f"{indir}/{sid}/{sid}.bam" for sid in samples]),
        batch_outdir = lambda wildcards: outdir + f"/{wildcards.genome}/cnv"
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("cnvkit_batch", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit batch for genome {wildcards.genome} at {current_time}")
            script = receipt_script(f"cnvkit_batch_{wildcards.genome}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(params.batch_outdir, f"cnvkit_batch.sh")
            os.makedirs(params.batch_outdir, exist_ok=True)
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"if [ -f {input.ref} ]; then\n")
                f.write(f"    {params.cnvkit} batch \\\n")
                f.write(f"        {params.treat_bams} \\\n")
                f.write(f"        --reference {input.ref} \\\n")
                f.write(f"        --processes {threads} \\\n")
                f.write(f"        --output-dir {params.batch_outdir} \\\n")
                f.write(f"        --scatter --diagram \\\n")
                f.write(f"        >> {log} 2>&1\n")
                f.write(f"else\n")
                f.write(f"    {params.cnvkit} batch \\\n")
                f.write(f"        {params.treat_bams} \\\n")
                f.write(f"        --method {params.method} \\\n")
                f.write(f"        --fasta {input.fasta} \\\n")
                f.write(f"        --access {input.access} \\\n")
                f.write(f"        --processes {threads} \\\n")
                f.write(f"        --output-dir {params.batch_outdir} \\\n")
                f.write(f"        --scatter --diagram \\\n")
                f.write(f"        >> {log} 2>&1\n")
                f.write(f"fi\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit batch for genome {wildcards.genome}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule cnvkit_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        cnr = expand(outdir + "/{{genome}}/cnv/{sample_id}.cnr", sample_id=samples),
        cns = expand(outdir + "/{{genome}}/cnv/{sample_id}.cns", sample_id=samples)
