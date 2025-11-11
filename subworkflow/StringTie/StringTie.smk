SNAKEFILE_FULL_PATH_StringTie = workflow.snakefile
SNAKEFILE_DIR_StringTie = os.path.dirname(SNAKEFILE_FULL_PATH_StringTie)
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_StringTie ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
StringTieYaml = get_yaml_path("StringTie")
configfile: StringTieYaml
logging.info(f"Include Align config: {StringTieYaml}")
logging.info(f"genomes:{genomes}, samples: {samples}")
rule stringTie:
    input:
        bam = outdir + "/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam"
    output:
        gtf = outdir + "/stringTie/{sample_id}/{genome}/{sample_id}.gtf"
    log:
        log = outdir + "/log/{genome}/{sample_id}/stringTie.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        gtf = lambda wildcards: config['stringTie'][wildcards.genome]['gtf'] #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
    threads: 5
    shell:
        """
        stringtie -o {output.gtf} {input.bam} -G {params.gtf} -p {threads} > {log.log} 2>&1
        """

rule stringTieMerge:
    input:
        gtf = expand(outdir + "/stringTie/{sample_id}/{genome}/{sample_id}.gtf",sample_id=samples,genome=genomes)
    output:
        outfile = outdir + "/stringTie/{genome}.gtf"
    log:
        log = outdir + "/log/{genome}/stringTieMerge.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        gtf = lambda wildcards: config['stringTie'][wildcards.genome]['gtf'] #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
    shell:
        """
        stringtie --merge {input.gtf} -o {output.outfile} -G {params.gtf} > {log.log} 2>&1
        """

rule TEcountStringTie:
    input:
        bam = outdir + "/counts/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam",
        gtf = outdir + "/stringTie/{genome}.gtf"
    output:
        project = outdir + "/counts/{sample_id}/{genome}/{sample_id}TEcountStringTie.cntTable"
    params:
        project = "{sample_id}TEcountStringTie",
        outdir = outdir + "/counts/{sample_id}/{genome}",
        TE_gtf = lambda wildcards: config['TEtranscripts'][wildcards.genome]['TE_gtf'],
    threads:5 #防止过多并行运行爆内存
    log:
        log = outdir + "/log/{genome}/{sample_id}/TEtranscriptsStringTie.log"
    conda:
        config['conda']['RNA-SNP']
    shell:
        """
        TEcount --sortByPos --format BAM --mode multi \
        -b {input.bam} --GTF {input.gtf} --TE {params.TE_gtf} \
        --project {params.project} --outdir {params.outdir} \
        > {log.log} 2>&1
        """

rule combine_TEStringtie:
    input:
        fileList = expand(outdir + "/counts/{sample_id}/{genome}/{sample_id}TEcountStringTie.cntTable",sample_id=samples,genome=genomes)
    output:
        outfile = outdir + "/counts/{genome}TEcountStringTie.cntTable"
    log:
        log = outdir + "/log/{genome}/combine_TEcountStringTie.log"
    params:
        combineTE = "scripts/combineTE.py",
        indir = outdir + "/counts"
    conda:
        config['conda']['RNA-SNP']
    shell:
        """
        python {params.combineTE} -p TEcountStringTie -i {params.indir} -o {output.outfile} > {log.log} 2>&1
        """

rule getStringtieBed:
    input:
        gtf = outdir + "/stringTie/{genome}.gtf",
        infile = outdir + "/counts/{genome}TEcountStringTie.cntTable"
    output:
        genefile = outdir + "/stringTie/{genome}_STG.bed",
        TEfile = outdir + "/stringTie/{genome}_TE.bed"
    log:
        log = outdir + "/log/{genome}/getStringtieBed.log"
    conda:
        config['conda']['RNA-SNP']
    threads:2
    params:
        script = "scripts/SNP/getBed.py",
        TE_gtf = lambda wildcards: config['TEtranscripts'][wildcards.genome]['TE_gtf']
    shell:
        """
        python {params.script} \
            --mode StringTie \
            --input {input.infile} \
            --output {output.genefile} \
            --output {output.TEfile} \
            --Gtf {input.gtf} \
            --TEGtf {params.TE_gtf} > {log.log} 2>&1
        """
rule StgTEOverlap:
    input:
        genefile = outdir + "/stringTie/{genome}_STG.bed",
        TEfile = outdir + "/stringTie/{genome}_TE.bed"
    output:
        outfile = outdir + "/stringTie/{genome}_StgTEOverlap.bed"
    log:
        log = outdir + "/log/{genome}_StgTEOverlap.log"
    conda:
        config['conda']['RNA-SNP']
    threads:2
    shell:
        """
        bedtools intersect -a {input.genefile} -b {input.TEfile} -wa -wb > {output.outfile} 2>{log.log}
        """
