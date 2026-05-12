from sanemake.logging import logger
indir = config.get("indir","input")
outdir = config.get("outdir", "output")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule bam_flagstat:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"
    log:
        logdir + "/{sample_id}/flagstat.log"
    conda:
        "samtools.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    shell:
        """
        {params.samtools} flagstat {input.bam} > {output.flagstat}
        """ 

rule samtools_result:
    input:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"

