include: "../common/common.smk"
from snakemake.logging import logger
import time
import os

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
control_samples = config.get("control_samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

def get_cnvkit_bams(wildcards):
    """Get BAM files for CNVkit batch analysis."""
    treat_bams = [f"{indir}/{sid}/{sid}.bam" for sid in samples]
    control_bams = [f"{indir}/{sid}/{sid}.bam" for sid in control_samples] if control_samples else []
    return {"treat": treat_bams, "control": control_bams}

rule cnvkit_reference:
    """Build CNVkit reference from control samples."""
    input:
        control_bams = expand(f"{indir}/{sid}/{sid}.bam", sid=control_samples),
        fasta = config.get("genome", {}).get("fasta"),
        access = config.get("genome", {}).get("access")
    output:
        ref = outdir + "/reference/reference.cnn"
    log:
        logdir + "/cnvkit_reference.log"
    threads:
        config.get("Params", {}).get("cnvkit", {}).get("processes") or 8
    conda:
        "cnvkit.yaml"
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs",
        control_bams = lambda wildcards: " ".join([f"{indir}/{sid}/{sid}.bam" for sid in control_samples])
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("cnvkit_reference", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit reference build at {current_time}")
            sample_outdir = os.path.dirname(str(output.ref))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"cnvkit_reference_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"{params.cnvkit} batch \\\n")
                f.write(f"    {params.control_bams} \\\n")
                f.write(f"    --method {params.method} \\\n")
                f.write(f"    --fasta {input.fasta} \\\n")
                f.write(f"    --access {input.access} \\\n")
                f.write(f"    --processes {threads} \\\n")
                f.write(f"    --output-reference {output.ref} \\\n")
                f.write(f"    > {log} 2>&1\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during cnvkit reference build: {e}\n")
            logger.error(f"Error occurred during cnvkit reference build: {e}")
            raise e

rule cnvkit_batch:
    """Run CNVkit batch analysis on all samples."""
    input:
        treat_bams = expand(f"{indir}/{sid}/{sid}.bam", sid=samples),
        ref = outdir + "/reference/reference.cnn" if control_samples else [],
        fasta = config.get("genome", {}).get("fasta"),
        access = config.get("genome", {}).get("access")
    output:
        cnr = expand(outdir + "/cnv/{sample_id}.cnr", sample_id=samples),
        cns = expand(outdir + "/cnv/{sample_id}.cns", sample_id=samples)
    log:
        logdir + "/cnvkit_batch.log"
    threads:
        config.get("Params", {}).get("cnvkit", {}).get("processes") or 8
    conda:
        "cnvkit.yaml"
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs",
        treat_bams = lambda wildcards: " ".join([f"{indir}/{sid}/{sid}.bam" for sid in samples]),
        outdir = outdir + "/cnv"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("cnvkit_batch", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit batch at {current_time}")
            script = os.path.join(params.outdir, f"cnvkit_batch_{current_time}.sh")
            os.makedirs(params.outdir, exist_ok=True)
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"if [ -f {input.ref} ]; then\n")
                f.write(f"    {params.cnvkit} batch \\\n")
                f.write(f"        {params.treat_bams} \\\n")
                f.write(f"        --reference {input.ref} \\\n")
                f.write(f"        --processes {threads} \\\n")
                f.write(f"        --output-dir {params.outdir} \\\n")
                f.write(f"        --scatter --diagram \\\n")
                f.write(f"        > {log} 2>&1\n")
                f.write(f"else\n")
                f.write(f"    {params.cnvkit} batch \\\n")
                f.write(f"        {params.treat_bams} \\\n")
                f.write(f"        --method {params.method} \\\n")
                f.write(f"        --fasta {input.fasta} \\\n")
                f.write(f"        --access {input.access} \\\n")
                f.write(f"        --processes {threads} \\\n")
                f.write(f"        --output-dir {params.outdir} \\\n")
                f.write(f"        --scatter --diagram \\\n")
                f.write(f"        > {log} 2>&1\n")
                f.write(f"fi\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during cnvkit batch: {e}\n")
            logger.error(f"Error occurred during cnvkit batch: {e}")
            raise e

rule cnvkit_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        cnr = expand(outdir + "/cnv/{sample_id}.cnr", sample_id=samples),
        cns = expand(outdir + "/cnv/{sample_id}.cns", sample_id=samples)
