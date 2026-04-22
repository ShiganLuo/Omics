from snakemake.logging import logger
indir = config.get('indir', "input")
outdir = config.get('outdir', "output")
logdir = config.get('logdir', "log")
wig_bin = config.get('wig',{}).get('bin', 25)
rule dedup:
    input:
        bam = indir + "/{sample_id}.bam"
    output:
        bam = outdir + "/dedup/{sample_id}.dedup.bam",
        bai = outdir + "/dedup/{sample_id}.dedup.bam.bai",
    log:
        logdir + "/{sample_id}/samtools_dedup.log"
    threads: 12
    conda:
        "igv.yaml"
    params:
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    shell:
        """
        {params.samtools} sort -n -@ {threads} {input.bam} \
        | {params.samtools} fixmate -m - - \
        | {params.samtools} sort -@ {threads} - \
        | {params.samtools} markdup -r -@ {threads} - {output.bam} 2>{log} && \
        {params.samtools} index -@ {threads} {output.bam} >> {log} 2>&1
        """

rule wig:
    input:
        bam = outdir + "/dedup/{sample_id}.dedup.bam",
        bai = outdir + "/dedup/{sample_id}.dedup.bam.bai"
    output:
        bigwig = outdir + "/{sample_id}.bigwig"
    log:
        log = logdir + "/{sample_id}/wig.log"
    conda:
        "igv.yaml"
    threads: 12 
    params:
        bin=wig_bin,
        deepTools = config.get('Procedure',{}).get('deepTools') or 'deepTools'
    shell:
        """
        bamCoverage --binSize {params.bin} --numberOfProcessors {threads} --extendReads --normalizeUsing CPM -b {input.bam} -o {output.bigwig} > {log.log} 2>&1 
        """
rule igv_result:
    input:
        bigwig = outdir + "/{sample_id}.bigwig"
