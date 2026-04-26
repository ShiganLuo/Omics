from typing import List
from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
paired_samples =  config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
genome_pairs: List[str] = config.get("genome_pairs", [])
genomeA, genomeB = genome_pairs


def get_inputFile_for_ngs_disambiguate(wildcards):
    logger.info(f"[get_inputFile_for_ngs_disambiguate] called with wildcards: {wildcards}")
    return {
        "bamA": f"{indir}/{genomeA}/{wildcards.sample_id}.bam",
        "bamB": f"{indir}/{genomeB}/{wildcards.sample_id}.bam"
    }

rule ngs_disambiguate:
    input:
        bamA = lambda wc: get_inputFile_for_ngs_disambiguate(wc)["bamA"],
        bamB = lambda wc: get_inputFile_for_ngs_disambiguate(wc)["bamB"]
    output:
        raw_bamA = temp(outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesA.bam"),
        raw_bamB = temp(outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesB.bam"),
        raw_ambiguousA = temp(outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesA.bam"),
        raw_ambiguousB = temp(outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesB.bam"),
        summary = outdir + "/{sample_id}/{sample_id}_summary.tsv"
    params:
        bamA_sortN = temp(outdir + "/{sample_id}/{sample_id}.bamA.sortN.bam"),
        bamB_sortN = temp(outdir + "/{sample_id}/{sample_id}.bamB.sortN.bam"),
        outdir = lambda wc: f"{outdir}/{wc.sample_id}",
        aligner = config.get("Params", {}).get("ngs_disambiguate", {}).get("aligner") or "hisat2",
        ngs_disambiguate = config.get("Procedure", {}).get("ngs_disambiguate") or "ngs_disambiguate",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    threads: 4
    conda:
        "disambiguate.yaml"
    log:
        logdir + "/{sample_id}/ngs_disambiguate.log"
    shell:
        """
        {params.samtools} sort -n -@ {threads} -o {params.bamA_sortN} {input.bamA} 2>> {log}
        {params.samtools} sort -n -@ {threads} -o {params.bamB_sortN} {input.bamB} 2>> {log}

        {params.ngs_disambiguate} \
            -s {wildcards.sample_id} \
            -o {params.outdir} \
            -a {params.aligner} \
            {params.bamA_sortN} \
            {params.bamB_sortN} \
            >> {log} 2>&1
        rm -f {params.bamA_sortN} {params.bamB_sortN}
        """

rule disambiguate_sort_rename:
    input:
        summary = outdir + "/{sample_id}/{sample_id}_summary.tsv",
        raw_bamA = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesA.bam",
        raw_bamB = outdir + "/{sample_id}/{sample_id}.disambiguatedSpeciesB.bam",
        raw_ambiguousA = outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesA.bam",
        raw_ambiguousB = outdir + "/{sample_id}/{sample_id}.ambiguousSpeciesB.bam"
    output:
        clean_bamA = outdir + "/{sample_id}/{sample_id}" + f".disambiguatedSpecies_{genome_pairs[0]}.bam",
        clean_bamB = outdir + "/{sample_id}/{sample_id}" + f".disambiguatedSpecies_{genome_pairs[1]}.bam",
        ambiguous_bamA = outdir + "/{sample_id}/{sample_id}" + f".ambiguousSpecies_{genome_pairs[0]}.bam",
        ambiguous_bamB = outdir + "/{sample_id}/{sample_id}" + f".ambiguousSpecies_{genome_pairs[1]}.bam",
        clean_summary = outdir + "/{sample_id}/{sample_id}_summary_renamed.tsv"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
        speciesA = genome_pairs[0],
        speciesB = genome_pairs[1]
    threads: 4
    conda:
        "disambiguate.yaml"
    log:
        logdir + "/{sample_id}/sort_rename.log"
    shell:
        """
        sed '1s/unique species A pairs/unique species {params.speciesA} pairs/; \
            1s/unique species B pairs/unique species {params.speciesB} pairs/' {input.summary} > {output.clean_summary}
        {params.samtools} sort -@ {threads} -o {output.clean_bamA} {input.raw_bamA}
        {params.samtools} sort -@ {threads} -o {output.clean_bamB} {input.raw_bamB}

        {params.samtools} sort -@ {threads} -o {output.ambiguous_bamA} {input.raw_ambiguousA}
        {params.samtools} sort -@ {threads} -o {output.ambiguous_bamB} {input.raw_ambiguousB}

        {params.samtools} index {output.clean_bamA}
        {params.samtools} index {output.clean_bamB}
        rm -f {input.summary}
        """

rule disambiguate_report:
    input:
        reports = expand(outdir + "/{sample_id}/{sample_id}_summary_renamed.tsv", sample_id=paired_samples + single_samples)
    output:
        report = outdir + "/disambiguate_qc.tsv"
    params:
        combine_script = ROOT_DIR + "/modules/disambiguate/combineDisambiguateQC.py"
    log:
        logdir + "/disambiguate_report.log"
    conda:
        "disambiguate.yaml"
    threads: 1
    shell:
        """
        python {params.combine_script} \
            -i {input.reports} \
            -o {output.report} \
            >> {log} 2>&1
        """
rule ngs_disambiguate_result:
    input:
        bamA = lambda wc: f"{outdir}/{wc.sample_id}/{wc.sample_id}.disambiguatedSpecies_{genome_pairs[0]}.bam",
        bamB = lambda wc: f"{outdir}/{wc.sample_id}/{wc.sample_id}.disambiguatedSpecies_{genome_pairs[1]}.bam"
