
import logging
from typing import Tuple
logger = logging.getLogger(__name__)
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
genome_pairs: Tuple[str, str] = config.get("genome_pairs", ())
genomeA, genomeB = genome_pairs
def get_inputFile_for_ngs_disambiguate(wildcards):
    logger.info(f"[get_inputFile_for_ngs_disambiguate] called with wildcards: {wildcards}")
    return {
        "bamA": f"{indir}/{genomeA}/{wildcards.sample_id}.bam",
        "bamB": f"{indir}/{genomeB}/{wildcards.sample_id}.bam"
    }

rule ngs_disambiguate:
    input:
        bamA = lambda wildcards: get_inputFile_for_ngs_disambiguate(wildcards)["bamA"],
        bamB = lambda wildcards: get_inputFile_for_ngs_disambiguate(wildcards)["bamB"]
    output:
        clean_bamA = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesA.bam",
        clean_bamB = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesB.bam",
        ambiguous_bamA = outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesA.bam",
        ambiguous_bamB = outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesB.bam",
        summary = outdir + "/{sample_id}/{sample_id}_summary.txt"
    params:
        outdir = lambda wildcards: f"{outdir}/{wildcards.sample_id}",
        aligner = config.get("bam", {}).get("aligner") or "hisat2",
        ngs_disambiguate = config.get("Procedure", {}).get("ngs_disambiguate") or "ngs_disambiguate"
    conda:
        "disambiguate.yaml"
    threads: 4 # for memory limit, we set threads to 4, but it can be adjusted according to the actual situation
    log:
        logdir + "/{sample_id}/ngs_disambiguate.log"
    shell:
        """
        {params.ngs_disambiguate} \
            -s {wildcards.sample_id} \
            -o {params.outdir} \
            -a {params.aligner} \
            {input.bamA} \
            {input.bamB} \
            > {log} 2>&1
        """

rule ngs_disambiguate_result:
    input:
        bamA = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesA_{genomeA}.bam",
        bamB = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesB_{genomeB}.bam"