from typing import Tuple
from snakemake.logging import logger
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "output/disambiguate")
ROOT_DIR = config.get("ROOT_DIR", "./")
genome_pairs: Tuple[str, str] = config.get("genome_pairs", ())
genomeA, genomeB = genome_pairs
paired_samples = config.get('paired_samples', [])
single_samples = config.get('single_samples', [])

def get_input_for_TEcount(wildcards):
    logger.info(f"[get_input_for_TEcount] called with wildcards: {wildcards}")
    if wildcards.genome == genomeA:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    elif wildcards.genome == genomeB:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    else:
        raise ValueError(f"wildcards.genome {wildcards.genome} is not in genome_pairs {genome_pairs}")
    return bam

rule TEcount:
    input:
        bamA = get_input_for_TEcount
    output:
        project = outdir + "/TEcount/{genome}/{sample_id}.TEcount.cntTable"
    params:
        project = lambda wildcards: f"{wildcards.sample_id}.TEcount",
        outdir = lambda wildcards: outdir + f"/TEcount/{wildcards.genome}",
        TE_gtf = lambda wildcards: config['genome'][wildcards.genome]['TE_gtf'],
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf'],
        TEcount = config.get('Procedure',{}).get('TEcount') or 'TEcount'
    log:
        logdir + "/{sample_id}/{genome}/TEcount.log"
    threads: 2
    conda:
        "../TEtranscripts.yaml"
    shell:
        """
        {params.TEcount} --sortByPos --format BAM --mode multi \
        -b {input.bamA} --GTF {params.gtf} --TE {params.TE_gtf} \
        --project {params.project} --outdir {params.outdir} \
        > {log} 2>&1
        """

def get_cntTable_for_TEcount(wildcards):
    logger.info(f"[get_cntTable_for_TEcount] called with wildcards: {wildcards}")
    cntTable = []
    for sample_id in single_samples:
        cntTable.append(f"{outdir}/TEcount/{wildcards.genome}/{sample_id}.TEcount.cntTable")
    for sample_id in paired_samples:
        cntTable.append(f"{outdir}/TEcount/{wildcards.genome}/{sample_id}.TEcount.cntTable")
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,genome: {wildcards.genome}")
    return cntTable

rule combine_TEcount:
    input:
        get_cntTable_for_TEcount
    output:
        outfile_id = outdir + "/TEcount/{genome}/all_TEcount.tsv",
        outfile_name = outdir + "/TEcount/{genome}/all_TEcount_name.tsv"
    conda:
        "../TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR +"/modules/TEtranscripts/bin/combineTE.py",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py",
        indir = outdir + "/TEcount/{genome}",
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    threads: 1
    log:
        logdir + "/all/TEtranscripts/{genome}_combine_TEcount.log"
    shell:
        """
        python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile_id} > {log} 2>&1
        python {params.geneId2Name} -c {output.outfile_id} -g {params.gtf} -o {output.outfile_name} >> {log} 2>&1
        """

def get_input_for_TElocal(wildcards):
    logger.info(f"[get_input_for_TElocal] called with wildcards: {wildcards}")
    if wildcards.genome == genomeA:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    elif wildcards.genome == genomeB:
        bam = indir + f"/{wildcards.sample_id}/{wildcards.sample_id}.disambiguatedSpecies_{wildcards.genome}.bam"
    else:
        raise ValueError(f"wildcards.genome {wildcards.genome} is not in genome_pairs {genome_pairs}")
    return bam

rule TElocal:
    input:
        bam = get_input_for_TElocal
    output:
        project = outdir + "/TElocal/{genome}/{sample_id}.TElocal.cntTable"
    log:
        logdir + "/{sample_id}/{genome}/TElocal.log"
    params:
        project = lambda wildcards: f"{wildcards.sample_id}.TElocal",
        TE = lambda wildcards: config['genome'][wildcards.genome]['TEind'],
        GTF = lambda wildcards: config['genome'][wildcards.genome]['gtf'],
        TElocal = config.get('Procedure',{}).get('TElocal') or 'TElocal'
    threads: 2
    conda:
        "../TEtranscripts.yaml"
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
        cntTable.append(f"{outdir}/TElocal/{wildcards.genome}/{sample_id}.TElocal.cntTable")
    for sample_id in paired_samples:
        cntTable.append(f"{outdir}/TElocal/{wildcards.genome}/{sample_id}.TElocal.cntTable")
    
    if len(cntTable) == 0:
        raise ValueError(f"rule combine_TElocal didn't get any input files,genome: {wildcards.genome}")
    return cntTable

rule combine_TElocal:
    input:
        get_cntTable_for_TElocal
    output:
        outfile_id = outdir + "/TElocal/{genome}/all_TElocal.tsv",
        outfile_name = outdir + "/TElocal/{genome}/all_TElocal_name.tsv"
    conda:
        "../TEtranscripts.yaml"
    params:
        combineTE = ROOT_DIR +"/modules/TEtranscripts/bin/combineTE.py",
        geneId2Name = ROOT_DIR +"/modules/TEtranscripts/bin/geneId2Name.py",
        indir = outdir + "/TElocal/{genome}",
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    threads: 1
    log:
        logdir + "/all/TEtranscripts/{genome}_combine_TElocal.log"
    shell:
        """
        python {params.combineTE} -p TElocal -i {params.indir} -o {output.outfile_id} > {log} 2>&1
        python {params.geneId2Name} -c {output.outfile_id} -g {params.gtf} -o {output.outfile_name} >> {log} 2>&1
        """

rule TEtranscripts_result:
    input:
        TEcount = outdir + "/TEcount/{genome}/all_TEcount.tsv",
        TElocal = outdir + "/TElocal/{genome}/all_TElocal.tsv"

