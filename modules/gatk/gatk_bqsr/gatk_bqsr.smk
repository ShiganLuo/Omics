from snakemake.logging import logger
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
interval = config.get("reference", {}).get("interval")
known_sites = config.get("reference", {}).get("known_sites", [])


rule BaseRecalibrator:
    input:
        bam = indir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        bai = indir + "/{sample_id}/{sample_id}.sorted_markdup.bai",
        ref = fasta,
        interval = interval
    output:
        table = outdir + "/bqsr/{sample_id}.recal_data.table"
    log:
        logdir + "/{sample_id}/gatk_bqsr.log"
    threads: 8
    conda:
        "../gatk.yaml"
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        known_sites = " ".join([f"--known-sites {ks}" for ks in known_sites])
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" BaseRecalibrator \
            -R {input.ref} --input {input.bam} \
            {params.known_sites} \
            {('-L ' + input.interval) if input.interval else ''} \
            -O {output.table}
        """


rule ApplyBQSR:
    input:
        bam = indir + "/{sample_id}/{sample_id}.sorted_markdup.bam",
        table = outdir + "/{sample_id}/{sample_id}.recal_data.table"
        ref = fasta
    output:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam"
    log:
        logdir + "/{sample_id}/ApplyBQSR.log"
    conda:
        "../gatk.yaml"
    threads:
        8
    params:
        javaOptions = "-Xmx30G",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        interval = interval
    shell:
        """
        {params.gatk}  --java-options "{params.javaOptions}" ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            -bqsr {input.table} \
            {('-L ' + params.interval) if params.interval else ''} \
            -O {output.bam}
        """
rule gatk_bqsr_result:
    input:
        bam = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam",
        bai = outdir + "/{sample_id}/{sample_id}.sorted_markdup.bqsr.bam.bai"
