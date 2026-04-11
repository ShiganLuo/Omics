import os

def get_fastqc_input(wildcards):
    pe1 = os.path.join(outdir, "cutadapt", f"{wildcards.sample_id}_1.fq.gz")
    pe2 = os.path.join(outdir, "cutadapt", f"{wildcards.sample_id}_2.fq.gz")
    se  = os.path.join(outdir, "cutadapt", f"{wildcards.sample_id}.fq.gz")

    if os.path.exists(pe1) and os.path.exists(pe2):
        return [str(pe1), str(pe2)]
    else:
        return [str(se)]


rule fastqc:
    input:
        get_fastqc_input
    output:
        directory(outdir + "/fastqc/{sample_id}")
    params:
        fastqc = config.get("Procedure", {}).get("fastqc", "fastqc")
    threads: 6
    conda:
        "fastqc.yaml"
    log:
        outdir + "/log/qc/{sample_id}/fastqc.txt"
    shell:
        """
        mkdir -p {output}
        {params.fastqc} \
            --threads {threads} \
            -o {output} \
            {input} \
            > {log} 2>&1
        """

rule fastqc_result:
    input:
        directory(outdir + "/fastqc/{sample_id}")