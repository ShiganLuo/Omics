from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")
mode = config.get("mode") or None
def get_input_for_trimming_Paired(wildcards):
    logger.info(f"Getting input for trimming_Paired with mode: {mode}")
    if mode == "UMI":
        return [
            f"{indir}/{wildcards.sample_id}_1.umi.fq.gz",
            f"{indir}/{wildcards.sample_id}_2.umi.fq.gz"
        ]
    else:
        return [
            f"{indir}/{wildcards.sample_id}_1.fq.gz",
            f"{indir}/{wildcards.sample_id}_2.fq.gz"
        ]

rule trimming_Paired:
    input:
        get_input_for_trimming_Paired
    output:
        fastq1 = temp(outdir + "/{sample_id}/{sample_id}_1.fq.gz"),
        fastq2 = temp(outdir + "/{sample_id}/{sample_id}_2.fq.gz"),
        report1 = outdir + "/{sample_id}/trimming_statistics_1.txt",
        report2 = outdir + "/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir = outdir + "/{sample_id}",
        quality = config.get('Params',{}).get("trim_galore", {}).get('quality') or 25,
        trim_galore = config.get('Procedure',{}).get('trim_galore') or 'trim_galore'
    threads: 6
    conda:
        "cutadapt.yaml"
    log:
        log = logdir + "/{sample_id}/trimming.txt"
    shell:
        """
        # trim_galore can automatically judge the fq quality scoring system,it's no need to add such as --phred33 --phred64
        {params.trim_galore} --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input[0]} {input[1]} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1} 2>{log.log}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2} 2>{log.log}
        suffix1=$(basename {input[0]})
        suffix2=$(basename {input[1]})
        mv {params.outdir}/${{suffix1}}_trimming_report.txt {output.report1} 2>{log.log}
        mv {params.outdir}/${{suffix2}}_trimming_report.txt {output.report2} 2>{log.log}
        """

def get_input_for_trimming_Single(wildcards):
    logger.info(f"Getting input for trimming_Single with mode: {mode}")
    if mode == "UMI":
        return f"{indir}/{wildcards.sample_id}.umi.single.fq.gz",
    else:
        return f"{indir}/{wildcards.sample_id}.single.fq.gz",

rule trimming_Single:
    input:
        fastq = get_input_for_trimming_Single
    output:
        fastq = temp(outdir + "/{sample_id}/{sample_id}.single.fq.gz"),
        report = outdir + "/{sample_id}/trimming_statistics.txt"
    params:
        outdir = outdir + "/{sample_id}",
        quality = config.get('Params',{}).get("cutadapt", {}).get('quality') or 25,
        trim_galore = config.get('Procedure',{}).get('trim_galore') or 'trim_galore'
    threads: 6
    conda:
        "cutadapt.yaml"
    log:
        log = logdir + "/{sample_id}/trimming.txt"
    shell:
        """
        {params.trim_galore} --phred33  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_trimmed.fq.gz {output.fastq} 2>{log.log}
        suffix=$(basename {input[0]})
        mv {params.outdir}/${{suffix}}_trimming_report.txt {output.report} 2>{log.log}
        """



