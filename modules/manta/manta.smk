from snakemake.logging import logger
import time
import os

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

rule manta_config:
    """Configure Manta SV caller."""
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        fasta = config.get("genome", {}).get("fasta")
    output:
        config = outdir + "/{sample_id}/runWorkflow.py.config.pickle",
        script = outdir + "/{sample_id}/runWorkflow.py"
    log:
        logdir + "/{sample_id}/manta_config.log"
    params:
        manta_outdir = outdir + "/{sample_id}"
    shell:
        """
        configManta.py \
            --bam {input.bam} \
            --referenceFasta {input.fasta} \
            --runDir {params.manta_outdir} \
            > {log} 2>&1
        """

rule manta_run:
    """Run Manta SV calling."""
    input:
        script = outdir + "/{sample_id}/runWorkflow.py",
        config = outdir + "/{sample_id}/runWorkflow.py.config.pickle"
    output:
        vcf = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz",
        tbi = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz.tbi"
    log:
        logdir + "/{sample_id}/manta_run.log"
    threads:
        config.get("Params", {}).get("manta", {}).get("threads") or 8
    shell:
        """
        python {input.script} -j {threads} > {log} 2>&1
        """

rule manta_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        vcf = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz",
        tbi = outdir + "/{sample_id}/results/variants/candidateSV.vcf.gz.tbi"
