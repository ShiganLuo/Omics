from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule tabix_bgzip:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf"
    output:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.vcf.gz"
    log:
        logdir + "/{sample_id}/bgzip.log"
    conda:
        "tabix.yaml"
    params:
        bgzip = config.get("Procedure", {}).get("bgzip") or "bgzip"
    shell:
        """
        mkdir -p $(dirname {output.vcf_gz})
        cp {input.vcf} {output.vcf_gz}.tmp.vcf
        {params.bgzip} -f {output.vcf_gz}.tmp.vcf
        mv {output.vcf_gz}.tmp.vcf.gz {output.vcf_gz}
        """

rule tabix_index:
    input:
        vcf_gz = indir + "/{sample_id}/{sample_id}.vcf.gz"
    output:
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
    log:
        logdir + "/{sample_id}/tabix.log"
    conda:
        "tabix.yaml"
    params:
        tabix = config.get("Procedure", {}).get("tabix") or "tabix"
    shell:
        """
        {params.tabix} -p vcf {input.vcf_gz} > {log} 2>&1
        """

rule tabix_result:
    input:
        vcf_gz = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
