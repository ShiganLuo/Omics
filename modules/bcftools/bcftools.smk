from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule bcftools_sort:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf"
    output:
        vcf = outdir + "/sort/{sample_id}/{sample_id}.sorted.vcf.gz"
    log:
        logdir + "/{sample_id}/bcftools_sort.log"
    conda:
        "bcftools.yaml"
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        {params.bcftools} sort \\
            -O z \\
            -o {output.vcf} \\
            {input.vcf} \\
            > {log} 2>&1
        """

rule bcftools_index_csi:
    input:
        vcf = indir + "/{sample_id}/{sample_id}.vcf.gz"
    output:
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
    log:
        logdir + "/{sample_id}/bcftools_index.log"
    conda:
        "bcftools.yaml"
    params:
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
    shell:
        """
        {params.bcftools} index {input.vcf} > {log} 2>&1
        """

rule bcftools_result:
    input:
        vcf = outdir + "/sort/{sample_id}/{sample_id}.sorted.vcf.gz",
        csi = outdir + "/{sample_id}/{sample_id}.vcf.gz.csi"
