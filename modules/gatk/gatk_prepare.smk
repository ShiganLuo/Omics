from snakemake.logging import logger
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
fasta = config.get("genome",{}).get("fasta")

rule gatk_index:
    input:
        fasta = fasta
    output:
        fai_index = outdir + "/index/genome.fa.fai",
        dict_index = outdir + "/index/genome.dict",
        fasta_link = outdir + "/index/genome.fa"
    log:
        logdir + "/index/gatk_index.log"
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
        bam = indir + "/{sample_id}/{sample_id}.bam",
    output:
        bam = temp(outdir + "/RG/{sample_id}/{sample_id}.bam"),
        bai = temp(outdir + "/RG/{sample_id}/{sample_id}.bam.bai")
    log:
        logdir + "/{sample_id}/addReadsGroup.log"
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
        bam = outdir + "/RG/{sample_id}/{sample_id}.bam",
        bai = outdir + "/RG/{sample_id}/{sample_id}.bam.bai"
    output:
        bam = outdir + "/bam-sorted-Markdup/{sample_id}/{sample_id}.bam",
        bai = outdir + "/bam-sorted-Markdup/{sample_id}/{sample_id}.bai",
        metrics = outdir + "/bam-sorted-Markdup/{sample_id}/{sample_id}_Markdup-metrics.txt"
    log:
        logdir + "/{sample_id}/MarkDuplicates.log"
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
        bam = outdir + "/bam-sorted-Markdup/{sample_id}/{sample_id}.bam",
        bai = outdir + "/bam-sorted-Markdup/{sample_id}/{sample_id}.bai",
        metrics = outdir + "/bam-sorted-Markdup/{sample_id}/{sample_id}_Markdup-metrics.txt"




