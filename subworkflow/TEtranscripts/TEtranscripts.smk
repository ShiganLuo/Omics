SNAKEFILE_FULL_PATH_TEtranscripts = workflow.snakefile
SNAKEFILE_DIR_TEtranscripts = os.path.dirname(SNAKEFILE_FULL_PATH_TEtranscripts)
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_TEtranscripts ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
TEtranscriptsYaml = get_yaml_path("TEtranscripts")
configfile: TEtranscriptsYaml
logging.info(f"Include TEtranscripts config: {TEtranscriptsYaml}")
rule TEtranscript_prepare:
    input:
        get_alignment_input,
        genome_index = lambda wildcards: config['genome'][wildcards.genome]['genome_index']
    output:
        outfile = temp(outdir + "/counts/bam/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
    log:
        log = outdir + "/log/TEtranscript/{genome}/{sample_id}/TEtranscript_prepare.log",
        STAR_log = outdir + "/log/TEtranscript/{genome}/{sample_id}/{sample_id}Log.out",
        STAR_progress = outdir + "/log/TEtranscript/{genome}/{sample_id}/{sample_id}Log.progress.out",
        STAR_final = outdir + "/log/TEtranscript/{genome}/{sample_id}/{sample_id}Log.final.out"
    threads: 12
    params:
        outPrefix = outdir + "/counts/bam/{genome}/{sample_id}",
        STAR = config['Procedure']['STAR'],
        # 动态判断输入参数,加上genome_index，如果三个参数，即为双端测序，两个参数即为单端测序
        input_params = lambda wildcards, input: \
            f"{input[0]} {input[1]}" if len(input) == 3 else f"{input[0]}"
    shell:
        """
        {params.STAR} --runThreadN {threads} \
            --genomeDir {input.genome_index} \
            --readFilesCommand zcat \
            --readFilesIn {params.input_params} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.outPrefix} \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 100  > {log.log} 2>&1
        mv {params.outPrefix}Log.out {log.STAR_log}
        mv {params.outPrefix}Log.progress.out {log.STAR_progress}
        mv {params.outPrefix}Log.final.out {log.STAR_final}
        """

rule TEcount:
    input:
        bam = outdir + "/counts/bam/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/counts/TEcount/{genome}/{sample_id}TEcount.cntTable"
    params:
        project = "{sample_id}TEcount",
        outdir = outdir + "/counts/TEcount/{genome}",
        TE_gtf = lambda wildcards: config['TEtranscripts'][wildcards.genome]['TE_gtf'],
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    log:
        outdir + "/log/TEtranscript/{genome}/{sample_id}/TEcount.log"
    conda:
        config['conda']['run']
    shell:
        """
        TEcount --sortByPos --format BAM --mode multi \
        -b {input.bam} --GTF {params.gtf} --TE {params.TE_gtf} \
        --project {params.project} --outdir {params.outdir} \
        > {log} 2>&1
        """

def get_cntTable_for_TEcount(wildcards):
    logging.info(f"[get_cntTable_for_TEcount] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id, genome in single_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/counts/TEcount/{genome}/{sample_id}TEcount.cntTable")
    for sample_id, genome in paired_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/counts/TEcount/{genome}/{sample_id}TEcount.cntTable")
    
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,genome: {wildcards.genome}\nsingle_sample_genome_pairs:{single_sample_genome_pairs}\npaired_sample_genome_pairs:{paired_sample_genome_pairs}")
    return cntTable

rule combine_TEcount:
    input:
        fileList = get_cntTable_for_TEcount
    output:
        outfile = outdir + "/counts/TEcount/{genome}/all_TEcount.cntTable"
    conda:
        config['conda']['run']
    params:
        combineTE = SNAKEFILE_DIR + "/utils/combineTE.py",
        indir = outdir + "/counts/TEcount/{genome}"
    log:
        outdir + "/log/TEtranscript/{genome}/combine_TEcount.log"
    shell:
        """
        python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile} > {log} 2>&1
        """

rule TElocal:
    input:
        bam = outdir + "/counts/bam/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/counts/TElocal/{genome}/{sample_id}TElocal.cntTable"
    log:
        outdir + "/log/TEtranscript/{genome}/{sample_id}/TElocal.log"
    params:
        project = "{sample_id}TElocal",
        TE = lambda wildcards: config['TElocal'][wildcards.genome]['TEind'],
        GTF = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    conda:
        config['conda']['run']
    shell:
        """
        TElocal --sortByPos -b {input.bam} \
        --GTF {params.GTF} --TE {params.TE} \
        --project {params.project} > {log} 2>&1
        mv {params.project}.cntTable {output.project}
        """

def get_cntTable_for_TElocal(wildcards):
    logging.info(f"[get_cntTable_for_TElocal] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id, genome in single_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/counts/TElocal/{genome}/{sample_id}TElocal.cntTable")
    for sample_id, genome in paired_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/counts/TElocal/{genome}/{sample_id}TElocal.cntTable")
    
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,genome: {wildcards.genome}\nsingle_sample_genome_pairs:{single_sample_genome_pairs}\npaired_sample_genome_pairs:{paired_sample_genome_pairs}")
    return cntTable

rule combine_TElocal:
    input:
        fileList = get_cntTable_for_TElocal
    output:
        outfile = outdir + "/counts/TElocal/{genome}/all_TElocal.cntTable"
    conda:
        config['conda']['run']
    params:
        combineTE = SNAKEFILE_DIR + "/utils/combineTE.py",
        indir = outdir + "/counts/TElocal/{genome}"
    log:
        outdir + "/log/TEtranscript/{genome}/combine_TElocal.log"
    shell:
        """
        python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile} > {log} 2>&1
        """