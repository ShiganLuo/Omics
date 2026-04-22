from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")
rule trimming_Paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        fastq1 = temp(outdir + "/{sample_id}_1.fq.gz"),
        fastq2 = temp(outdir + "/{sample_id}_2.fq.gz"),
        report1 = logdir + "/{sample_id}/trimming_statistics_1.txt",
        report2 = logdir + "/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir = outdir,
        quality = 35,
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
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fq.gz_trimming_report.txt {output.report2}
        """

rule trimming_Single:
    input:
        fastq = indir + "/{sample_id}.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}.single.fq.gz"),
        report = logdir + "/{sample_id}/trimming_statistics.txt"
    params:
        outdir = outdir,
        quality = 25,
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
        mv {params.outdir}/{wildcards.sample_id}_trimmed.fq.gz {output.fastq}
        mv {params.outdir}/{wildcards.sample_id}.fq.gz_trimming_report.txt {output.report}
        """


rule triming_paired_reslut:
    input:
        fastq1 = outdir + "/{sample_id}_1.fq.gz",
        fastq2 = outdir + "/{sample_id}_2.fq.gz",
        report1 = logdir + "/{sample_id}/trimming_statistics_1.txt",
        report2 = logdir + "/{sample_id}/trimming_statistics_2.txt",

rule triming_single_result:
    input:
        fastq = outdir + "/{sample_id}.fq.gz",
        report = logdir + "/{sample_id}/trimming_statistics.txt"
