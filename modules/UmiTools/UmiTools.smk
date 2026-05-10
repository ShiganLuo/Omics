from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")

rule umi_tools_dedup_for_hisat2:
    input:
        bam = indir + "/{sample_id}.bam"
    output:
        bam = temp(outdir + "/{sample_id}.dedup.bam"),
        bai = temp(outdir + "/{sample_id}.dedup.bam.bai"),
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        method = config.get('Params',{}).get('umi_tools',{}).get('method') or 'unique',
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    threads: 2
    conda:
        "UmiTools.yaml"
    log:
        log = logdir + "/{sample_id}/umi_tools_dedup_run.txt"
    shell:
        """
        {params.umi_tools} dedup \
            --method={params.method} \
            -I {input.bam} \
            -S {output.bam} \
            > {log} 2>&1
        {params.samtools} index {output.bam}
        """

rule umi_tools_dedup_for_star:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = temp(outdir + "/{sample_id}.dedup.bam"),
        bai = temp(outdir + "/{sample_id}.dedup.bam.bai")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        method = config.get('Params',{}).get('umi_tools',{}).get('method') or 'unique',
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    threads: 2
    conda:
        "UmiTools.yaml"
    log:
        log = logdir + "/{sample_id}/umi_tools_dedup_run.txt"
    shell:
        """
        {params.umi_tools} dedup --method={params.method} \
            -I {input.bam} -S {output.bam} > {log} 2>&1
        {params.samtools} index {output.bam} 2>{log}
        """
