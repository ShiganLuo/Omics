from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
# first col: target(human) genome,second col: contaminating genome. human sample may contaminated by mouse genome

def get_inputFile_for_XenofilterR(wildcards):
    logger.info(f"[get_inputFile_for_XenofilterR] called with wildcards: {wildcards}")
    row = [
        f"{indir}/{wildcards.pollution_source_genome}/{wildcards.sample_id}.bam",
        f"{indir}/{wildcards.host_genome}/{wildcards.sample_id}.bam"    
        ]
    return row

rule XenofilterR:
    input:
        bams = get_inputFile_for_XenofilterR
    output:
        csvIn = outdir + "/xenofilterR/{sample_id}/{sample_id}.csv",
        outBam = temp(outdir + "/xenofilterR/{sample_id}/{sample_id}_Filtered.bam"),
        outBai = temp(outdir + "/xenofilterR/{sample_id}/{sample_id}_Filtered.bam.bai")
    log:
        outdir + "/log/XenofilterR/{sample_id}/XenofilterR.log"
    threads: 8
    params:
        csv_content = lambda wildcards, input: ",".join(input.bams),
        outdir = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.sample_id}",
        outSampleName = lambda wildcards: wildcards.sample_id,
        tempBam = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.host_genome}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam",
        tempBai = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.host_genome}/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam.bai",
        MM = 8,
        script = SNAKEFILE_DIR + "/utils/XenofilteR.r",xenofilterR/{wildcards.sample_id}/Filtered_bams
        Rscript = config.get('Procedure',{}).get('Rscript') or 'Rscript'
    conda:
        "XenofilterR.yaml"
    shell:
        """
        echo "{params.csv_content}" > {output.csvIn}
        # rename ignorme .bam
        {params.Rscript} {params.script} \
            --inputFile {output.csvIn} \
            --outputDir {params.outdir} \
            --renameSamples {params.outSampleName} \
            --MM {params.MM} \
            --workers 1 > {log} 2>&1
        mv {params.tempBam} {output.outBam} # XenofilteR would run failed if it find Filtered_bams dir exist
        mv {params.tempBai} {output.outBai}
        """

