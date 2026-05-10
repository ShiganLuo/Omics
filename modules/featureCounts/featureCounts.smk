from snakemake.logging import logger
indir = config.get("indir", "output")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", "./")
single_samples = config.get("single_samples", [])
paired_samples = config.get("paired_samples", [])

# need test
def get_bams_for_featureCounts_single(wildcards):
    logger.info(f"[get_bams_for_featureCounts_single] called with wildcards: {wildcards}")
    bams = []
    for sample_id in single_samples:
        bams.append(f"{indir}/{sample_id}/{sample_id}.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_single_noMultiple didn't get any input bams")
    return bams

rule featureCounts_single_noMultiple:
    input:
        bams = get_bams_for_featureCounts_single
    output:
        outfile = outdir + "/all_single_featureCounts.tsv"
    log:
        logdir + "/all/featureCounts/featureCounts_single_noMultiple.log"
    conda: 
        "featureCounts.yaml"
    threads:
        10
    params:
        featureCounts = config.get('Procedure',{}).get('featureCounts') or 'featureCounts',
        gtf = config["genome"]["gtf"]
    shell:
        """
        {params.featureCounts} -T {threads} -t exon -g gene_id -a {params.gtf} -o {output.outfile} {input.bams} > {log} 2>&1
        """

def get_bams_for_featureCounts_paired(wildcards):
    logger.info(f"[get_bams_for_featureCounts_paired] called with wildcards: {wildcards}")
    bams = []
    for sample_id in paired_samples:
        bams.append(f"{indir}/{sample_id}/{sample_id}.bam")

    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_paired_noMultiple didn't get any input bams")
    return bams

rule featureCounts_paired_noMultiple:
    input:
        bams = get_bams_for_featureCounts_paired
    output:
        outfile = outdir + "/all_paired_featureCounts.tsv",
    log:
        logdir + "/all/featureCounts/featureCounts_paired_noMultiple.log"
    conda:
        "featureCounts.yaml"
    threads:
        10
    params:
        featureCounts = config.get('Procedure',{}).get('featureCounts') or 'featureCounts',
        gtf = config["genome"]["gtf"]
    shell:
        """
        # for multiple -M -O
        {params.featureCounts} -T {threads} -B -p --countReadPairs -t exon -g gene_id -a {params.gtf} -o {output.outfile} {input.bams} > {log} 2>&1
        """


rule featureCounts_result:
    input:
        paired = outdir + "all_paired_featureCounts.tsv",
        single = outdir + "all_single_featureCounts.tsv"
