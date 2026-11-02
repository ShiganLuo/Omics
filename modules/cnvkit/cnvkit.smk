include: "../common/common.smk"

from snakemake.logging import logger
import time
import os

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
bam_suffix = config.get("bam_suffix", "")
samples = config.get("samples", [])
control_samples = config.get("control_samples", [])
access = config.get("genome", {}).get("access")
ROOT_DIR = config.get("ROOT_DIR", ".")

def get_cnvkit_bams(wildcards):
    """Get BAM files for CNVkit batch analysis."""
    treat_bams = [f"{indir}/{sid}/{sid}{bam_suffix}.bam" for sid in samples]
    control_bams = [f"{indir}/{sid}/{sid}{bam_suffix}.bam" for sid in control_samples] if control_samples else []
    return {"treat": treat_bams, "control": control_bams}

rule cnvkit_reference:
    """Build CNVkit reference from control samples."""
    input:
        control_bams = expand(indir + "/{sid}/{sid}" + bam_suffix + ".bam", sid=control_samples),
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "cnvkit"),
        access = lambda wildcards: [access] if access else [],
    output:
        ref = outdir + "/reference/reference.cnn"
    log:
        logdir + "/cnvkit_reference.log"
    threads:
        config.get("Params", {}).get("cnvkit", {}).get("processes") or 8
    resources:
        mem_mb = rule_res("cnvkit_reference").get("mem_mb", 16000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_reference")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs",
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("cnvkit_reference", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit reference build at {current_time}")
            sample_outdir = os.path.dirname(str(output.ref))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"cnvkit_reference_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"cnvkit_reference.sh")
            cmd = [params.cnvkit, "batch"]
            cmd += [str(bam) for bam in input.control_bams]
            cmd += ["--method", params.method, "--fasta", str(input.fasta)]
            if input.access:
                cmd += ["--access", str(input.access[0])]
            cmd += ["--processes", str(threads), "--output-reference", str(output.ref)]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit reference build: {e}")
            raise RuntimeError(f"Error occurred during cnvkit reference build: {e}, you can check the log file {log_path} for more details.")

rule cnvkit_batch:
    """Run CNVkit batch analysis on all samples."""
    input:
        treat_bams = expand(indir + "/{sid}/{sid}" + bam_suffix + ".bam", sid=samples),
        ref = lambda wildcards: [outdir + "/reference/reference.cnn"] if control_samples else [],
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "cnvkit"),
        access = lambda wildcards: [access] if access else [],
    output:
        cnr = expand(outdir + "/cnv/{sample_id}.cnr", sample_id=samples),
        cns = expand(outdir + "/cnv/{sample_id}.cns", sample_id=samples)
    log:
        logdir + "/cnvkit_batch.log"
    threads:
        config.get("Params", {}).get("cnvkit", {}).get("processes") or 8
    resources:
        mem_mb = rule_res("cnvkit_batch").get("mem_mb", 32000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_batch")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs",
        outdir = outdir + "/cnv"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("cnvkit_batch", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit batch at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"cnvkit_batch_{current_time}.sh")
            plain_script = os.path.join(params.outdir, f"cnvkit_batch.sh")
            os.makedirs(params.outdir, exist_ok=True)
            cmd = [params.cnvkit, "batch"]
            cmd += [str(bam) for bam in input.treat_bams]
            if input.ref:
                cmd += ["--reference", str(input.ref[0])]
            else:
                cmd += ["--method", params.method, "--fasta", str(input.fasta)]
                if input.access:
                    cmd += ["--access", str(input.access[0])]
            cmd += ["--processes", str(threads), "--output-dir", params.outdir, "--scatter", "--diagram"]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit batch: {e}")
            raise RuntimeError(f"Error occurred during cnvkit batch: {e}, you can check the log file {log_path} for more details.")

rule cnvkit_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        cnr = expand(outdir + "/cnv/{sample_id}.cnr", sample_id=samples),
        cns = expand(outdir + "/cnv/{sample_id}.cns", sample_id=samples)
    resources:
        mem_mb = rule_res("cnvkit_result").get("mem_mb", 8000),
