from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai")

rule deepvariant_run:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        bai = indir + "/{sample_id}/{sample_id}.bam.bai",
        fasta = fasta,
        fai = fai
    output:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi",
        gvcf = outdir + "/{sample_id}/{sample_id}.g.vcf.gz"
    log:
        logdir + "/{sample_id}/deepvariant.log"
    threads: 8
    conda:
        "deepvariant.yaml"
    params:
        deepvariant = config.get("Procedure", {}).get("deepvariant") or "run_deepvariant",
        model_type = config.get("Params", {}).get("deepvariant", {}).get("model_type") or "PACBIO",
        outdir_sample = outdir + "/{sample_id}"
    shell:
        """
        mkdir -p {params.outdir_sample}
        {params.deepvariant} \\
            --num_shards {threads} \\
            --model_type {params.model_type} \\
            --ref {input.fasta} \\
            --reads {input.bam} \\
            --output_vcf {output.vcf} \\
            --output_gvcf {output.gvcf} \\
            > {log} 2>&1
        """

rule deepvariant_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
