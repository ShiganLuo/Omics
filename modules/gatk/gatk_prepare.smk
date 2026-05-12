outdir = config.get("outdir") or "output"
fasta = config.get("fasta")
################################################################
# SNV (indexing, addReadsGroup, MarkDuplicates)
################################################################
rule gatk_index:
    input:
        fasta = fasta
    output:
        fai_index = outdir + "/gatk/index/{genome}/{genome}.fa.fai",
        dict_index = outdir + "/gatk/index/{genome}/{genome}.dict",
        fasta_link = outdir + "/gatk/index/{genome}/{genome}.fa"
    log:
        outdir + "/log/gatk/{genome}/{genome}_indexing.log"
    threads: 4
    conda:
        "gatk.yaml"
    params:
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
    shell:
        """
        ln -s {input.fasta} {output.fasta_link}
        {params.gatk} CreateSequenceDictionary -R {input.fasta} -O {output.dict_index} > {log} 2>&1
        samtools faidx {input.fasta} -o {output.fai_index} >> {log} 2>&1
        """

rule addReadsGroup:
    input:
        bam = outdir + "/samtools/bam-sorted/{genome}/{sample_id}.bam",
    output:
        bam = temp(outdir + "/gatk/RG/{genome}/{sample_id}.bam"),
        bai = temp(outdir + "/gatk/RG/{genome}/{sample_id}.bam.bai")
    log:
        outdir + "/log/gatk/{genome}/{sample_id}/addReadsGroup.log"
    threads: 8
    conda:
        "gatk.yaml"
    params:
        id = "{sample_id}",
        javaOptions = "--java-options -Xmx15G",
        RGLB = config.get("addReadsGroup", {}).get("RGLB") or "lib1",
        RGPL = config.get("addReadsGroup", {}).get("RGPL") or "illumina",
        RGPU = config.get("addReadsGroup", {}).get("RGPU") or "unit1",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    shell:
        """
        echo "sample_id: {wildcards.sample_id}" > {log}
        {params.gatk} AddOrReplaceReadGroups {params.javaOptions} \
            --INPUT {input.bam} --OUTPUT {output.bam} \
            -SO coordinate --RGLB {params.RGLB} --RGPL {params.RGPL} --RGPU {params.RGPU} --RGSM {params.id} >> {log} 2>&1
        {params.samtools} index -@ {threads} {output.bam} >> {log} 2>&1
        """


rule MarkDuplicates:
    input:
        bam = outdir + "/gatk/RG/{genome}/{sample_id}.bam",
        bai = outdir + "/gatk/RG/{genome}/{sample_id}.bam.bai"
    output:
        bam = outdir + "/gatk/bam-sorted-Markdup/{genome}/{sample_id}.bam",
        bai = outdir + "/gatk/bam-sorted-Markdup/{genome}/{sample_id}.bai",
        metrics = outdir + "/gatk/bam-sorted-Markdup/{genome}/{sample_id}_Markdup-metrics.txt"
    log:
        outdir + "/log/gatk/{genome}/{sample_id}/MarkDuplicates.log"
    threads: 8
    conda:
        "gatk.yaml"
    params:
        javaOptions = "-Xms20g -Xmx30g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk = config.get("Procedure", {}).get("gatk") or "gatk"
    shell:
        """
        {params.gatk} --java-options "{params.javaOptions}" MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam}   \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT \
            --METRICS_FILE {output.metrics} > {log} 2>&1
        """

rule gatk_prepare_result:
    input:
        bam = outdir + "/gatk/bam-sorted-Markdup/{genome}/{sample_id}.bam",
        bai = outdir + "/gatk/bam-sorted-Markdup/{genome}/{sample_id}.bai",
        metrics = outdir + "/gatk/bam-sorted-Markdup/{genome}/{sample_id}_Markdup-metrics.txt"




