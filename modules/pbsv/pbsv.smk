from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")

rule pbsv_discover:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        bai = indir + "/{sample_id}/{sample_id}.bam.bai"
    output:
        svsig = outdir + "/discover/{sample_id}/{sample_id}.svsig.gz"
    log:
        logdir + "/{sample_id}/pbsv_discover.log"
    threads: 4
    conda:
        "pbsv.yaml"
    params:
        pbsv = config.get("Procedure", {}).get("pbsv") or "pbsv"
    shell:
        """
        mkdir -p $(dirname {output.svsig})
        {params.pbsv} discover \\
            --num-threads {threads} \\
            {input.bam} \\
            {output.svsig} \\
            > {log} 2>&1
        """

rule pbsv_call:
    input:
        svsig = outdir + "/discover/{sample_id}/{sample_id}.svsig.gz",
        fasta = fasta
    output:
        vcf = outdir + "/call/{sample_id}/{sample_id}.vcf"
    log:
        logdir + "/{sample_id}/pbsv_call.log"
    threads: 8
    conda:
        "pbsv.yaml"
    params:
        pbsv = config.get("Procedure", {}).get("pbsv") or "pbsv"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        {params.pbsv} call \\
            --num-threads {threads} \\
            {input.fasta} \\
            {input.svsig} \\
            {output.vcf} \\
            > {log} 2>&1
        """

rule pbsv_result:
    input:
        vcf = outdir + "/call/{sample_id}/{sample_id}.vcf"
