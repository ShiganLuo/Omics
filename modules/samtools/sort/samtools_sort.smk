from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
fasta = config.get("genome", {}).get("fasta")

rule bam_sort:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = outdir + "/sort/{sample_id}/{sample_id}.sorted.bam"
    log:
        logdir + "/{sample_id}/samtools_sort.log"
    threads: 8
    conda:
        "../samtools.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        {params.samtools} sort \\
            -@ {threads} \\
            -o {output.bam} \\
            {input.bam} \\
            > {log} 2>&1
        """

rule bam_index:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bai = outdir + "/{sample_id}/{sample_id}.bam.bai"
    log:
        logdir + "/{sample_id}/samtools_index.log"
    threads: 4
    conda:
        "../samtools.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    shell:
        """
        {params.samtools} index \\
            -@ {threads} \\
            {input.bam} \\
            {output.bai} \\
            > {log} 2>&1
        """

rule samtools_sort_index_result:
    input:
        bam = outdir + "/sort/{sample_id}/{sample_id}.sorted.bam",
        bai = outdir + "/{sample_id}/{sample_id}.bam.bai"
