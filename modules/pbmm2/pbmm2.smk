from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")

rule pbmm2_align:
    input:
        bam = indir + "/{sample_id}.bam",
        fasta = fasta
    output:
        bam = outdir + "/{sample_id}/{sample_id}.bam"
    log:
        logdir + "/{sample_id}/pbmm2_align.log"
    threads: 16
    conda:
        "pbmm2.yaml"
    params:
        pbmm2 = config.get("Procedure", {}).get("pbmm2") or "pbmm2"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        {params.pbmm2} align \\
            --sort \\
            --num-threads {threads} \\
            {input.fasta} \\
            {input.bam} \\
            {output.bam} \\
            > {log} 2>&1
        """

rule pbmm2_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.bam"
