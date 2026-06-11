from snakemake.logging import logger
import time
import os

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})

def get_macs3_input(wildcards):
    """
    Get treatment (IP) and optional control (Input) BAM files for MACS3.
    sample_ip_input_map: dict mapping IP sample_id -> input sample_id (or None)
    """
    sample = wildcards.sample_id
    genome = wildcards.genome
    
    bam_treatment = f"{indir}/{sample}/{sample}.bam"
    
    # Check if there's a matched input control
    input_sample = sample_ip_input_map.get(sample)
    if input_sample:
        bam_control = f"{indir}/{input_sample}/{input_sample}.bam"
        return {
            "bam_treatment": bam_treatment,
            "bam_control": bam_control
        }
    
    return {"bam_treatment": bam_treatment}

rule macs3_callpeak:
    """
    MACS3 peak calling for ChIP-seq/DIP-seq data.
    Supports both with-control and without-control modes.
    """
    input:
        unpack(get_macs3_input)
    output:
        peak = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.narrowPeak",
        xls = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.xls"
    log:
        logdir + "/macs3/{genome}/{sample_id}.log"
    threads: 4
    conda:
        "macs3.yaml"
    params:
        macs3 = config.get("Procedure", {}).get("macs3") or "macs3",
        outdir = lambda wildcards: f"{outdir}/{wildcards.genome}/{wildcards.sample_id}",
        name = lambda wildcards: wildcards.sample_id,
        bw = config.get("Params", {}).get("macs3", {}).get("bw") or 200,
        pvalue = config.get("Params", {}).get("macs3", {}).get("pvalue") or "1e-5",
        genome_size = config.get("Params", {}).get("macs3", {}).get("genome_size") or "mm"
    shell:
        """
        {params.macs3} callpeak \
            -t {input.bam_treatment} \
            {('-c ' + input.bam_control) if 'bam_control' in input else ''} \
            --bw {params.bw} \
            -p {params.pvalue} \
            -g {params.genome_size} \
            --outdir {params.outdir} \
            --name {params.name} \
            --seed 2346 \
            &> {log}
        """

rule macs3_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        peak = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.narrowPeak",
        xls = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.xls"
