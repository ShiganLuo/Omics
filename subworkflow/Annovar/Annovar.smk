SNAKEFILE_FULL_PATH_Annovar = workflow.snakefile
SNAKEFILE_DIR_Annovar = os.path.dirname(SNAKEFILE_FULL_PATH_Annovar)
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_Annovar ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
AnnovarYaml = get_yaml_path("Annovar")
configfile: AnnovarYaml
logging.info(f"Include Align config: {AnnovarYaml}")
logging.info(f"main snakefile directory: {SNAKEFILE_DIR}")

rule TEcoutCPM:
    input:
       infile = outdir + "/counts/humanTEcount.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCPM.cntTable"
    log:
        log = outdir + "/log/human/TEcoutCPM.log"
    params:
        script = SNAKEFILE_DIR + "/scripts/SNP/run-NormCountMat.R",
        Rscript = config["Procedure"]["Rscript"]
    shell:
        """
        {params.Rscript} {params.script} -i {input.infile} -o {output.outfile} > {log.log} 2>&1
        """

rule commonExpression:
    input:
        infile = outdir + "/counts/humanTEcountCPM.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.cntTable"
    log:
        log = outdir + "/log/human/commonExpression.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = SNAKEFILE_DIR + "scripts/SNP/commonExpression.py"
    shell:
        """
            python {params.script} --input {input.infile} --output {output.outfile} --threshold 5 > {log.log} 2>&1
        """

rule getBed:
    input:
        infile = outdir + "/counts/humanTEcountCommon.cntTable"
    output:
        outfile = outdir + "/counts/humanTEcountCommon.bed"
    log:
        log = outdir + "/log/human/getBed.log"
    conda:
        config['conda']['RNA-SNP']
    params:
        script = SNAKEFILE_DIR + "scripts/SNP/getBed.py",
        gtf = config['bed']['human']['gtf'],
        TE_gtf = config['bed']['human']['TE_gtf']
    shell:
        """
            python {params.script} \
                --input {input.infile} \
                --output {output.outfile} \
                --Gtf {params.gtf} \
                --TEGtf {params.TE_gtf} > {log.log} 2>&1
        """

rule vcfIntersectBed:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}.vcf.gz",
        bed = outdir + "/counts/humanTEcountCommon.bed"
    output:
        outfile = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    log:
        log = outdir + "/log/human/{sample_id}/vcfIntersectBed.log"
    conda:
        config['conda']['RNA-SNP']
    threads:4 #防止同时执行太多，爆内存
    shell:
        """
            bedtools intersect -a {input.vcf} -b {input.bed} -wa -wb > {output.outfile} 2>{log.log}
        """

rule annovar_convert:
    input:
        vcf = outdir + "/filter/vcf/human/{sample_id}Common.vcf"
    output:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_convert.log"
    params:
        convert = "/opt/annovar/convert2annovar.pl",
        perl = config["Procedure"]["perl"]
    threads:4 #防止同时执行太多，爆内存
    shell:
        """
        {params.perl} {params.convert} \
        -format vcf4 \
        -withfreq {input.vcf} > {output.avinput} 2>{log.log}
        """

rule annovar_table:
    input:
        avinput = outdir + "/annovar/human/{sample_id}/{sample_id}.avinput"
    output:
        outfile = outdir + "/annovar/human/{sample_id}/{sample_id}.GRCh38_multianno.csv"
    log:
        log = outdir + "/log/human/{sample_id}/annovar_table.log"
    params:
        db = config['annovar']['human']['db'],
        buildver = config['annovar']['human']['buildver'],
        # annotate = "/opt/annovar/annotate_variation.pl",
        table = config["Procedure"]["table_annovar"]
        out = outdir + "/annovar/human/{sample_id}/{sample_id}"
    threads: 4 #防止同时执行太多，爆内存
    shell:
        """
        {params.perl} {params.table} \
        {input.avinput} {params.db} \
        -buildver {params.buildver} \
        -out {params.out} \
        -remove -protocol refGene \
        -operation g \
        -nastring . \
        -csvout > {log.log} 2>&1
        """