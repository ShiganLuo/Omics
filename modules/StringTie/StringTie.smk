from snakemake.logging import logger
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
samples = config.get("samples") or []
rule stringTie:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        gtf = outdir + "/{sample_id}/{sample_id}.gtf"
    log:
        logdir + "/{sample_id}/stringTie.log"
    conda:
        "StringTie.yaml"
    params:
        gtf = config.get('genome', {}).get('gtf'), #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    threads: 5
    shell:
        """
        {params.stringtie} -o {output.gtf} {input.bam} -G {params.gtf} -p {threads} > {log} 2>&1
        """

rule stringTieMerge:
    input:
        gtfs = expand(outdir + "/{sample_id}/{sample_id}.gtf",sample_id=samples)
    output:
        gtf = outdir + "/stringTieMerge.gtf"
    log:
        logdir + "/all/stringtie/stringTieMerge.log"
    conda:
        "StringTie.yaml"
    params:
        gtf = config.get('genome', {}).get('gtf'), #最好使用完整的gtf文件，更有利于准确判断是否是新转录本
        stringtie = config.get("Procedure", {}).get("stringtie") or "stringtie"
    shell:
        """
        {params.stringtie} --merge {input.gtfs} -o {output.gtf} -G {params.gtf} > {log} 2>&1
        """
