SNAKEFILE_FULL_PATH_XenofilterR = workflow.snakefile
SNAKEFILE_DIR_XenofilterR = os.path.dirname(SNAKEFILE_FULL_PATH_XenofilterR)
import csv
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_XenofilterR ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
XenofilterRYaml = get_yaml_path("XenofilterR")
configfile: XenofilterRYaml
logger.info(f"Include XenofilterR config: {XenofilterRYaml}")
logger.info(f"main snakefile excute path: {EXECUTION_DIR}")
logger.info(f"XenofilterR target samples: {XenofilterR_target_samples}\nXenofilterR_target_genome: {XenofilterR_target_genome}\nXenofilterR pollution source genome: {XenofilterR_pollution_source_genome}")
# first col: target(human) genome,second col: contaminating genome. human sample may contaminated by mouse genome

def get_inputFile_for_XenofilterR(wildcards):
    logger.info(f"[get_inputFile_for_XenofilterR] called with wildcards: {wildcards}")
    row = [
        f"{outdir}/2pass/{wildcards.sample_id}/{XenofilterR_target_genome}/{wildcards.sample_id}Aligned.sortedByCoord.out.bam",
        f"{outdir}/2pass/{wildcards.sample_id}/{XenofilterR_pollution_source_genome}/{wildcards.sample_id}Aligned.sortedByCoord.out.bam"
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
        tempBam = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam",
        tempBai = lambda wildcards: f"{outdir}/xenofilterR/{wildcards.sample_id}/Filtered_bams/{wildcards.sample_id}_Filtered.bam.bai",
        MM = 8,
        script = SNAKEFILE_DIR + "/utils/XenofilteR.r",
        Rscript = config["Procedure"]["Rscript"]
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

