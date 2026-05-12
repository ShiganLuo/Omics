interval = config.get("reference", {}).get("interval")
known_sites = config.get("reference", {}).get("known_sites", [])

################################################################
# SNV (BaseRecalibrator, ApplyBQSR) no longer used
################################################################

rule BaseRecalibrator:
    input:
        bam = outdir + "/{genome}/gatk/bam-sorted-Markdup/{sample_id}.bam",
        bai = outdir + "/{genome}/gatk/bam-sorted-Markdup/{sample_id}.bai",
        ref = fasta,
        interval = interval
    output:
        table = outdir + "/{genome}/gatk/bqsr/{sample_id}.recal_data.table"
    log:
        outdir + "/log/gatk/{genome}/{sample_id}/gatk_bqsr.log"
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
        bam = outdir + "/{genome}/gatk/bam-sorted-Markdup/{sample_id}.bam",
        table = outdir + "/{genome}/gatk/bqsr/{sample_id}.recal_data.table",
        ref = fasta
    output:
        bam = outdir + "/{genome}/gatk/bqsr/{sample_id}.sorted.markdup.BQSR.bam"
    log:
        outdir + "/log/gatk/{genome}/{sample_id}/ApplyBQSR.log"
    conda:
        "gatk.yaml"
    threads:
        8
    params:
        javaOptions = f"--java-options -Xmx30G -Djava.io.tmpdir={TMP_DIR}",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        interval = interval
    shell:
        """
        {params.gatk} {params.javaOptions} ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            -bqsr {input.table} \
            {('-L ' + params.interval) if params.interval else ''} \
            -O {output.bam}
        """
rule gatk_bqsr_result:
    input:
        bam = outdir + "/{genome}/gatk/bqsr/{sample_id}.sorted.markdup.BQSR.bam",
        bai = outdir + "/{genome}/gatk/bqsr/{sample_id}.sorted.markdup.BQSR.bam.bai"
