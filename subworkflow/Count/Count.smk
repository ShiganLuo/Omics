SNAKEFILE_FULL_PATH_Count = workflow.snakefile
SNAKEFILE_DIR_Count = os.path.dirname(SNAKEFILE_FULL_PATH_Count)




countYaml = get_yaml_path("Count")
configfile: countYaml
logging.info(f"Include Align config: {countYaml}")


def get_bams_for_featureCounts_single(wildcards):
    bams = []
    for sample_id, genome in single_sample_genome_pairs:
        if genome == wildcards.genome:
            bams.append(f"{outdir}/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_single_noMultiple didn't get any input bams, genome: {wildcards.genome},\nsingle_sample_genome_pairs: {single_sample_genome_pairs}")
    return bams

rule featureCounts_single_noMultiple:
    input:
        bams = get_bams_for_featureCounts_single
    output:
        outfile = outdir + "/counts/featureCounts/{genome}/{genome}_single_count.tsv"
    log:
        outdir + "/log/Align/{genome}_featureCounts_single_noMultiple.log"
    threads:
        20
    params:
        featureCounts = config['Procedure']['featureCounts'],
        gtf = lambda wildcards: config["genome"][wildcards.genome]["gtf"]
    shell:
        """
        {params.featureCounts} -T {threads} -t exon -g gene_id -a {params.gtf} -o {output.outfile} {input.bams} > {log} 2>&1
        """

def get_bams_for_featureCounts_paired(wildcards):
    bams = []
    for sample_id, genome in paired_sample_genome_pairs:
        if genome == wildcards.genome:
            bams.append(f"{outdir}/2pass/{sample_id}/{genome}/{sample_id}Aligned.sortedByCoord.out.bam")
    if len(bams) == 0:
        raise ValueError(f"rule featureCounts_paired_noMultiple didn't get any input bams, genome: {wildcards.genome},\npaired_sample_genome_pairs:{paired_sample_genome_pairs}")
    return bams

rule featureCounts_paired_noMultiple:
    input:
        bams = get_bams_for_featureCounts_paired
    output:
        outfile = outdir + "/counts/featureCounts/{genome}/{genome}_paired_count.tsv"
    log:
        outdir + "/log/Align/{genome}_featureCounts_paired_noMultiple.log"
    threads:
        20
    params:
        featureCounts = config['Procedure']['featureCounts'],
        gtf = lambda wildcards: config["genome"][wildcards.genome]["gtf"]
    shell:
        """
        # for multiple -M -O
        {params.featureCounts} -T {threads} -B -p --countReadPairs -t exon -g gene_id -a {params.gtf} -o {output.outfile} {input.bams} > {log} 2>&1
        """

rule featureCounts_result:
    input:
        paired = outdir + "/counts/featureCounts/{genome}/{genome}_paired_count.tsv",
        single = outdir + "/counts/featureCounts/{genome}/{genome}_single_count.tsv"
