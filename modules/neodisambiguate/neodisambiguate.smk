from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
# wait for test
def get_inputFile_for_neodisambiguate(wildcards):
    logger.info(f"[get_inputFile_for_neodisambiguate] called with wildcards: {wildcards}")

    return {
        "bam1": f"{indir}/{wildcards.genomeA}/{wildcards.sample_id}.bam",
        "bam2": f"{indir}/{wildcards.genomeB}/{wildcards.sample_id}.bam"
    }
rule neodisambiguate:
    input:
        **get_inputFile_for_neodisambiguate()
    output:
        clean_bam1 = outdir + "/{sample_id}/{sample_id}.{genomeA}.neodisambiguatedA.bam", # for wildcards, the order of genomeA and genomeB is not fixed, so we use .neodisambiguatedA and .neodisambiguatedB to represent the two output bam files
        clean_bam2 = outdir + "/{sample_id}/{sample_id}.{genomeB}.neodisambiguatedB.bam",
        ambiguous_bam1 = outdir + "/{sample_id}/ambiguous-alignments/{sample_id}.{genomeA}.ambiguousA.bam",
        ambiguous_bam2 = outdir + "/{sample_id}/ambiguous-alignments/{sample_id}.{genomeB}.ambiguousB.bam",
    params:
        prefix = lambda wildcards: f"{outdir}/{wildcards.sample_id}.neodisambiguated",
        neodisambiguate = config.get("Procedure", {}).get("neodisambiguate") or "neodisambiguate"
    threads: 8
    log:
        logdir + "/{sample_id}/neodisambiguate.log"
    shell:
        """
        {params.neodisambiguate} \
            -s {params.prefix} \
            -o {outdir} \
            {input.bam1} \
            {input.bam2} \
            > {output.summary} 2> {log}
        """