from snakemake.logging import logger
indir = config.get('indir', "input")
outdir = config.get('outdir', "output")
logdir = config.get('logdir', "log")
rule dedup_hisat2:
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
rule dedup_star:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
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
        binSize= config.get('Params',{}).get('bamCoverage',{}).get('binSize') or 50,
        bamCoverage = config.get('Procedure',{}).get('bamCoverage') or 'bamCoverage',
        normalizeUsing = config.get('Params', {}).get('bamCoverage',{}).get('normalizeUsing') or "CPM",
        offset = config.get('Params', {}).get('bamCoverage',{}).get('offset') or None,
        extendReads = config.get('Params', {}).get('bamCoverage',{}).get('extendReads') or False
    shell:
        """
        {params.bamCoverage} \
            --binSize {params.binSize} \
            --numberOfProcessors {threads} \
            --extendReads \
            --normalizeUsing {params.normalizeUsing} \
            --Offset {params.offset} \
            --extendReads {params.extendReads} \
            -b {input.bam} \
            -o {output.bigwig} > {log.log} 2>&1 
        """

rule igv_result:
    input:
        bigwig = outdir + "/{sample_id}.bigwig"
