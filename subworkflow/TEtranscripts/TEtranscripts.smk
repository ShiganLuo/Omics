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

def get_alignment_input(wildcards):
    """
    function: Dynamically determines the input file type: paired-end or single-end sequencing.
    Based on the paired_samples and single_samples lists.This function is called in the star_align rule.

    param: 
        wildcards: Snakemake wildcards object containing the sample_id.
        paired_samples = ['sample1', 'sample2', ...]
        single_samples = ['sample3', 'sample4', ...]
    These lists must be defined in the Snakefile or config file.

    return: A list of input file paths for the STAR alignment step. 
    """
    logging.info(f"[get_alignment_input] called with wildcards: {wildcards}")
    # 构造可能的输入路径
    paired_r1 = f"{outdir}/cutadapt/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{outdir}/cutadapt/{wildcards.sample_id}_2.fq.gz"
    single = f"{outdir}/cutadapt/{wildcards.sample_id}Single.fq.gz"
    
    # 检查文件实际存在情况
    if wildcards.sample_id in paired_samples:
        logging.info(f"双端测序：{[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        logging.info(f"单端测序：{[single]}")
        return [single]
    else:
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample_id}\n"
            f"Checked paths:\n- {paired_r1}\n- {paired_r2}\n- {single}"
        )
rule TEtranscript_prepare:
    input:
        get_alignment_input,
        genome_index = lambda wildcards: config['STAR'][wildcards.genome]['genome_index']
    output:
        outfile = outdir + "/counts/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    log:
        log=outdir+"/log/{genome}/{sample_id}/TEtranscript_prepare.log"
    threads: 15
    params:
        outPrefix = outdir + "/counts/{sample_id}/{genome}/{sample_id}",
        STAR = config["STAR"]["procedure"],
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
        """

rule TEcount:
    input:
        bam = outdir + "/counts/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/counts/{sample_id}/{genome}/{sample_id}TEcount.cntTable"
    params:
        project = "{sample_id}TEcount",
        outdir = outdir + "/counts/{sample_id}/{genome}",
        TE_gtf = lambda wildcards: config['TEtranscripts'][wildcards.genome]['TE_gtf'],
        gtf = lambda wildcards: config['TEtranscripts'][wildcards.genome]['gtf']
    log:
        log = outdir + "/log/{genome}/{sample_id}/TEtranscripts.log"
    conda:
        config['conda']['RNA-SNP']
    shell:
        """
        TEcount --sortByPos --format BAM --mode multi \
        -b {input.bam} --GTF {params.gtf} --TE {params.TE_gtf} \
        --project {params.project} --outdir {params.outdir} \
        > {log.log} 2>&1
        """

rule combine_TEcount:
    input:
        fileList = expand(outdir + "/counts/{sample_id}/{genome}/{sample_id}TEcount.cntTable",sample_id=samples)
    output:
        outfile = outdir + "/counts/{genome}TEcount.cntTable"
    conda:
        config['conda']['RNA-SNP']
    params:
        combineTE = "scripts/combineTE.py",
        indir = outdir + "/counts"
    log:
        log = outdir + "/log/human/combine_TEcount.log"
    shell:
        """
        python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile} > {log.log} 2>&1
        """

rule TElocal:
    input:
        bam = outdir + "/counts/{sample_id}/human/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/counts/{sample_id}/human/{sample_id}TElocal.cntTable"
    params:
        project = "{sample_id}TElocal",
        TE = config['TElocal']['human']['TEind'],
        GTF = config['TElocal']['human']['gtf'],
        procedure = "/opt/TElocal/TElocal"
    log:
        log = outdir+"/log/human/{sample_id}/TElocal.log"
    conda:
        config['conda']['TElocal']
    shell:
        """
        which python
        {params.procedure} --sortByPos -b {input.bam} \
        --GTF {params.GTF} --TE {params.TE} \
        --project {params.project} > {log.log} 2>&1
        mv {params.project}.cntTable {output.project}
        """

rule combine_TElocal:
    input:
        fileList = expand(outdir + "/counts/{sample_id}/human/{sample_id}TElocal.cntTable",sample_id=samples)
    output:
        outfile = outdir + "/counts/humanTElocal.cntTable"
    conda:
        config['conda']['RNA-SNP']
    params:
        combineTE = "scripts/combineTE.py",
        indir = outdir + "/counts"
    log:
        log = outdir + "/log/human/combine_TElocal.log"
    shell:
        """
        python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile} > {log.log} 2>&1
        """