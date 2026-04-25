import os
indir = config.get("indir", "output/raw_fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
log_suffix = config.get("log_suffix", "txt")
# need test
def get_fastqc_input(wildcards):
    pe1 = os.path.join(outdir, f"{wildcards.sample_id}_1.fq.gz")
    pe2 = os.path.join(outdir, f"{wildcards.sample_id}_2.fq.gz")
    se  = os.path.join(outdir, f"{wildcards.sample_id}.fq.gz")

    if os.path.exists(pe1) and os.path.exists(pe2):
        return [str(pe1), str(pe2)]
    else:
        return [str(se)]


rule fastqc:
    input:
        get_fastqc_input
    output:
        directory(outdir + "/{sample_id}")
    params:
        fastqc = config.get("Procedure", {}).get("fastqc") or "fastqc"
    threads: 6
    conda:
        "fastqc.yaml"
    log:
        logdir + "/{sample_id}/fastqc." + log_suffix
    shell:
        """
        {params.fastqc} \
            --threads {threads} \
            -o {output} \
            {input} \
            > {log} 2>&1
        """

rule fastqc_result:
    input:
        directory(outdir + "/fastqc/{sample_id}")
