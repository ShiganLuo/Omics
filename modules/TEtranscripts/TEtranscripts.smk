from snakemake.logging import logger
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "output/raw_fastq")
ROOT_DIR = config.get("ROOT_DIR", "./")
single_samples = config.get("single_samples", [])
paired_samples = config.get("paired_samples", [])


rule TEcount:
    input:
        bam = indir + "/{sample_id}.bam"
    output:
        project = outdir + "/TEcount/{sample_id}.TEcount.cntTable"
    params:
        project = "{sample_id}.TEcount",
        outdir = outdir + "/TEtranscripts/TEcount",
        TE_gtf = lambda wildcards: config['genome']['TE_gtf'],
        gtf = lambda wildcards: config['genome']['gtf'],
        TEcount = config.get('Procedure',{}).get('TEcount') or 'TEcount'
    log:
        logdir + "/{sample_id}/TEcount.log"
    conda:
        "TEtranscripts.yaml"
    shell:
        """
        {params.TEcount} --sortByPos --format BAM --mode multi \
        -b {input.bam} --GTF {params.gtf} --TE {params.TE_gtf} \
        --project {params.project} --outdir {params.outdir} \
        > {log} 2>&1
        """

def get_cntTable_for_TEcount(wildcards):
    logger.info(f"[get_cntTable_for_TEcount] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in single_samples:
        cntTable.append(f"{outdir}/TEcount/{sample_id}.TEcount.cntTable")
    for sample_id in paired_samples:
        cntTable.append(f"{outdir}/TEcount/{sample_id}.TEcount.cntTable")
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,single_samples:{single_samples},paired_samples:{paired_samples}")
    return cntTable

rule combine_TEcount:
    input:
        fileList = get_cntTable_for_TEcount
    output:
        outfile = outdir + "/TEcount/all_TEcount.tsv"
    conda:
        "TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR + "/TEtranscripts/bin/combineTE.py",
        indir = outdir + "/TEtranscripts/TEcount"
    log:
        logdir + "/TEtranscripts/combine_TEcount.log"
    shell:
        """
        python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile} > {log} 2>&1
        """

rule TElocal:
    input:
        bam = indir + "/{sample_id}.bam"
    output:
        project = outdir + "/TElocal/{sample_id}.TElocal.cntTable"
    log:
        logdir + "/{sample_id}/TElocal.log"
    params:
        project = "{sample_id}.TElocal",
        TE = lambda wildcards: config['genome']['TEind'],
        GTF = lambda wildcards: config['genome']['gtf'],
        TElocal = config.get('Procedure',{}).get('TElocal') or 'TElocal'
    threads: 2
    conda:
        "TEtranscripts.yaml"
    shell:
        """
        {params.TElocal} --sortByPos -b {input.bam} \
        --GTF {params.GTF} --TE {params.TE} \
        --project {params.project} > {log} 2>&1
        mv {params.project}.cntTable {output.project}
        """

def get_cntTable_for_TElocal(wildcards):
    logger.info(f"[get_cntTable_for_TElocal] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in single_samples:
        cntTable.append(f"{outdir}/TElocal/{sample_id}.TElocal.cntTable")
    for sample_id in paired_samples:
        cntTable.append(f"{outdir}/TElocal/{sample_id}.TElocal.cntTable")
    
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,single_samples:{single_samples},paired_samples:{paired_samples}")
    return cntTable

rule combine_TElocal:
    input:
        fileList = get_cntTable_for_TElocal
    output:
        outfile = outdir + "/TElocal/all_TElocal.tsv"
    conda:
        "TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR + "/utils/combineTE.py",
        indir = outdir + "/TEtranscripts/TElocal"
    log:
        outdir + "/log/TEtranscript/combine_TElocal.log"
    shell:
        """
        python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile} > {log} 2>&1
        """

rule TEtranscripts_result:
    input:
        TEcount = outdir + "/TEcount/all_TEcount.tsv",
        TElocal = outdir + "/TElocal/all_TElocal.tsv"


