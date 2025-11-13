SNAKEFILE_FULL_PATH_XenofilterR = workflow.snakefile
SNAKEFILE_DIR_XenofilterR = os.path.dirname(SNAKEFILE_FULL_PATH_XenofilterR)
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
# first col: target(human) genome,second col: contaminating genome. human sample may contaminated by mouse genome
rule generate_xenofilter_input:
    input:
        expand(outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam", 
               sample_id=all_samples, genome="Mus_musculus")
    output:
        csvIn = outdir + "/xenofilterR/xenofilterR_input.csv",
        # csvRe = outdir +"/xenofilterR/xenofilterR_reName.csv"
    run:
        import csv
        with open(output.csvIn, 'w', newline='') as f:
            writer = csv.writer(f)
            for sample in all_samples:
                row = [
                    f"{outdir}/2pass/{sample}/Homo_sapiens/{sample}Aligned.sortedByCoord.out.bam",
                    f"{outdir}/2pass/{sample}/Mus_musculus/{sample}Aligned.sortedByCoord.out.bam"
                ]
                writer.writerow(row)

rule XenofilterR:
    input:
        csvIn = outdir + "/xenofilterR/xenofilterR_input.csv",
    output:
        # expand(outdir + "/xenofilterR/Filtered_bams/{sample_id}_Filtered.bam",sample_id=all_samples),
        # expand(outdir + "/xenofilterR/Filtered_bams/{sample_id}_Filtered.bam.bai",sample_id=all_samples)
        outdir = directory(outdir + "/xenofilterR/bam") #XenofilteR设计不合理，没办法
    log:
        outdir + "/log/human/XenofilterR.log"
    threads: 6
    params:
        script = SNAKEFILE_DIR + "/utils/XenofilteR.r",
        threshold = 8,
        Rscript = config["Procedure"]["Rscript"]
    shell:
        """
        {params.Rscript} {params.script} \
            --inputFile {input.csvIn} \
            --outputDir {output.outdir} \
            --MM {params.threshold} \
            --workers {threads} > {log} 2>&1
        """