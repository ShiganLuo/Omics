from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")
rule UmiTools_extract_single:
    input:
        fastq = indir + "/{sample_id}.single.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}.umi.single.fq.gz")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        extract_method = config.get('Procedure',{}).get('extract_method') or 'string',
        bc_pattern = config.get('Params',{}).get('umi_tools',{}).get('bc_pattern') or 'NNNXXXXNN'
    log:
        logdir + "/{sample_id}/umi_tools_extract_single_run.txt"
    conda:
        "../UmiTools.yaml"
    threads: 2
    shell:
        """
        {params.umi_tools} extract \
            --extract-method={params.extract_method} \
            --bc-pattern={params.bc_pattern} \
            -I {input.fastq} \
            -S {output.fastq} \
            > {log} 2>&1
        """

rule UmiTools_extract_paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        fastq1 = temp(outdir + "/{sample_id}_1.umi.fq.gz"),
        fastq2 = temp(outdir + "/{sample_id}_2.umi.fq.gz")
    params:
        umi_tools = config.get('Procedure',{}).get('umi_tools') or 'umi_tools',
        extract_method = config.get('Procedure',{}).get('extract_method') or 'string',
        bc_pattern = config.get('Params',{}).get('umi_tools',{}).get('bc_pattern') or 'NNNXXXXNN',
        bc_pattern2 = config.get('Params',{}).get('umi_tools',{}).get('bc_pattern2') or 'NNNXXXXNN'
    log:
        logdir + "/{sample_id}/umi_tools_extract_paired_run.txt"
    threads: 2
    conda:
        "../UmiTools.yaml"
    shell:
        """
        {params.umi_tools} extract \
            --extract-method={params.extract_method} \
            --bc-pattern={params.bc_pattern} \
            -S {output.fastq1} \
            --bc-pattern2={params.bc_pattern2} \
            --read2-in={input.fastq2} \
            --read2-out={output.fastq2} \
            > {log} 2>&1
        """
