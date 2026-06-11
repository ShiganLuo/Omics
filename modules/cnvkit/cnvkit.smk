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
        control_bams = lambda wildcards: [f"{indir}/{sid}/{sid}.bam" for sid in control_samples],
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
    shell:
        """
        {params.cnvkit} batch \
            {params.control_bams} \
            --method {params.method} \
            --fasta {input.fasta} \
            --access {input.access} \
            --processes {threads} \
            --output-reference {output.ref} \
            > {log} 2>&1
        """

rule cnvkit_batch:
    """Run CNVkit batch analysis on all samples."""
    input:
        treat_bams = lambda wildcards: [f"{indir}/{sid}/{sid}.bam" for sid in samples],
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
    shell:
        """
        if [ -f {input.ref} ]; then
            {params.cnvkit} batch \
                {params.treat_bams} \
                --reference {input.ref} \
                --processes {threads} \
                --output-dir {params.outdir} \
                --scatter --diagram \
                > {log} 2>&1
        else
            {params.cnvkit} batch \
                {params.treat_bams} \
                --method {params.method} \
                --fasta {input.fasta} \
                --access {input.access} \
                --processes {threads} \
                --output-dir {params.outdir} \
                --scatter --diagram \
                > {log} 2>&1
        fi
        """

rule cnvkit_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        cnr = expand(outdir + "/cnv/{sample_id}.cnr", sample_id=samples),
        cns = expand(outdir + "/cnv/{sample_id}.cns", sample_id=samples)
