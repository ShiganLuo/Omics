SNAKEFILE_FULL_PATH_TEtranscripts = workflow.snakefile
SNAKEFILE_DIR_TEtranscripts = os.path.dirname(SNAKEFILE_FULL_PATH_TEtranscripts)
TEtranscriptsYaml = get_yaml_path("TEtranscripts",SNAKEFILE_DIR_TEtranscripts)
configfile: TEtranscriptsYaml
logger.info(f"Include TEtranscripts config: {TEtranscriptsYaml}")

rule TEcount:
    input:
        bam = outdir + "/TEtranscripts/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/TEtranscripts/TEcount/{genome}/{sample_id}TEcount.cntTable"
    params:
        project = "{sample_id}TEcount",
        outdir = outdir + "/TEtranscripts/TEcount/{genome}",
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
    logger.info(f"[get_cntTable_for_TEcount] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id, genome in single_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/TEtranscripts/TEcount/{genome}/{sample_id}TEcount.cntTable")
    for sample_id, genome in paired_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/TEtranscripts/TEcount/{genome}/{sample_id}TEcount.cntTable")
    
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,genome: {wildcards.genome}\nsingle_sample_genome_pairs:{single_sample_genome_pairs}\npaired_sample_genome_pairs:{paired_sample_genome_pairs}")
    return cntTable

rule combine_TEcount:
    input:
        fileList = get_cntTable_for_TEcount
    output:
        outfile = outdir + "/TEtranscripts/TEcount/{genome}/all_TEcount.cntTable"
    conda:
        config['conda']['run']
    params:
        combineTE = SNAKEFILE_DIR + "/utils/combineTE.py",
        indir = outdir + "/TEtranscripts/TEcount/{genome}"
    log:
        outdir + "/log/TEtranscript/{genome}/combine_TEcount.log"
    shell:
        """
        python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile} > {log} 2>&1
        """

rule TElocal:
    input:
        bam = outdir + "/TEtranscripts/{sample_id}/{genome}/{sample_id}.Aligned.sortedByCoord.out.bam"
    output:
        project = outdir + "/TEtranscripts/TElocal/{genome}/{sample_id}TElocal.cntTable"
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
    logger.info(f"[get_cntTable_for_TElocal] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id, genome in single_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/TEtranscripts/TElocal/{genome}/{sample_id}TElocal.cntTable")
    for sample_id, genome in paired_sample_genome_pairs:
        if genome == wildcards.genome:
            cntTable.append(f"{outdir}/TEtranscripts/TElocal/{genome}/{sample_id}TElocal.cntTable")
    
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,genome: {wildcards.genome}\nsingle_sample_genome_pairs:{single_sample_genome_pairs}\npaired_sample_genome_pairs:{paired_sample_genome_pairs}")
    return cntTable

rule combine_TElocal:
    input:
        fileList = get_cntTable_for_TElocal
    output:
        outfile = outdir + "/TEtranscripts/TElocal/{genome}/all_TElocal.cntTable"
    conda:
        config['conda']['run']
    params:
        combineTE = SNAKEFILE_DIR + "/utils/combineTE.py",
        indir = outdir + "/TEtranscripts/TElocal/{genome}"
    log:
        outdir + "/log/TEtranscript/{genome}/combine_TElocal.log"
    shell:
        """
        python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile} > {log} 2>&1
        """

rule TEtranscripts_result:
    input:
        TEcount = outdir + "/TEtranscripts/TEcount/{genome}/all_TEcount.cntTable",
        TElocal = outdir + "/TEtranscripts/TElocal/{genome}/all_TElocal.cntTable"


if config["Procedure"]["aligner"] == "star":
    include: "align_star.smk"
    logger.info("aligner: star, load align_star.smk")

elif config["Procedure"]["aligner"] == "hisat2":
    include: "align_hisat2.smk"
    logger.info("aligner: hisat2, load align_hisat2.smk")
else:
    # 默认使用star比对
    include: "align_star.smk"
    logger.info("default: load align_star.smk")