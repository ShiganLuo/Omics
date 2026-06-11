from snakemake.logging import logger
import time
import os

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

rule samtools_stats:
    """Get BAM statistics including fragment size distribution."""
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        stats = outdir + "/stats/{sample_id}/{sample_id}_stats.txt"
    log:
        logdir + "/{sample_id}/samtools_stats.log"
    conda:
        "fragment_size.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    shell:
        """
        {params.samtools} stats {input.bam} > {output.stats} 2>{log}
        """

rule getFragmentSize:
    """Extract and summarize fragment size distribution from all samples."""
    input:
        stats = expand(outdir + "/stats/{sample_id}/{sample_id}_stats.txt", sample_id=samples)
    output:
        hist = outdir + "/fragment/FragmentSize.txt"
    log:
        logdir + "/getFragmentSize.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/fragment_size/bin/getFragmentSize.py")
    conda:
        "fragment_size.yaml"
    shell:
        """
        python {params.script} --input {input.stats} --out {output.hist} > {log} 2>&1
        """

rule plotFragmentSize:
    """Plot fragment size distribution."""
    input:
        hist = outdir + "/fragment/FragmentSize.txt"
    output:
        png = outdir + "/fragment/FragmentSize.png"
    log:
        logdir + "/plotFragmentSize.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/fragment_size/bin/plotFragmentSize.py")
    conda:
        "fragment_size.yaml"
    shell:
        """
        python {params.script} --input {input.hist} --output {output.png} > {log} 2>&1
        """

rule fragment_size_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        hist = outdir + "/fragment/FragmentSize.txt",
        png = outdir + "/fragment/FragmentSize.png"
