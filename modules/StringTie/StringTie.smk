from snakemake.logging import logger
import time
indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
samples = config.get("samples") or []
ROOT_DIR = config.get("ROOT_DIR", ".")
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
rule TEChimericTranscripts:
    input:
        gtf = outdir + "/{sample_id}/{sample_id}.gtf"
    output:
        txt = outdir + "/{sample_id}/{sample_id}_TE_chimeric_transcripts.txt"
    log:
        logdir + "/{sample_id}/TEChimericTranscripts.log"
    conda:
        "StringTie.yaml"
    params:
        te_gtf = config.get('genome', {}).get('TE_gtf'),
        TEChimericTranscripts = ROOT_DIR + "/modules/StringTie/bin/TEChimericTranscripts.py"
    threads: 5
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/{wildcards.sample_id}/TEChimericTranscripts.{current_time}.sh"
        cmd = f"python {params.TEChimericTranscripts} -s {input.gtf} -t {params.te_gtf} -o {output.txt} > {log} 2>&1"
        with open(script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(cmd + "\n")
        shell(f"bash {script}")


rule stringTieMerge:
    input:
        gtfs = expand(outdir + "/{sample_id}/{sample_id}.gtf",sample_id=samples)
    output:
        gtf = outdir + "/stringtie_merged.gtf"
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
