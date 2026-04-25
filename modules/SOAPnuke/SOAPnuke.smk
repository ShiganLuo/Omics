from snakemake.logging import logger
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir= config.get("indir", "output/raw_fastq")

rule soapnuke_filter_paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        clean_fastq1 = outdir + "/{sample_id}_1.fq.gz",
        clean_fastq2 = outdir + "/{sample_id}_2.fq.gz"
    log:
        logdir + "/{sample_id}/soapnuke_filter_paired.log"
    params:
        min_len = 15, # 最短保留长度：15
        low_qual = 0.2, # 低质量碱基比例：0.2
        max_n = 0.05, # N 碱基比例：0.05
        qual_type = 2, # Phred 33 (-Q 2) default
        SOAPnuke = config.get('Procedure',{}).get('SOAPnuke') or 'SOAPnuke',
        workdir = lambda wildcards: f"{outdir}/{wildcards.sample_id}",
        clean_fastq1 = lambda wildcards: f"{wildcards.sample_id}_1.fq.gz",
        clean_fastq2 = lambda wildcards: f"{wildcards.sample_id}_2.fq.gz",
    conda:
        "SOAPnuke.yaml"
    threads: 8
    shell:
        """
        mkdir -p {params.workdir}
        {params.SOAPnuke} filter \
            -l {params.min_len} \
            -q {params.low_qual} \
            -n {params.max_n} \
            -1 {input.fastq1} \
            -2 {input.fastq2} \
            -C {params.clean_fastq1} \
            -D {params.clean_fastq2} \
            -Q {params.qual_type} \
            -o {params.workdir} \
            -T {threads} \
            2>{log}
        mv {params.workdir}/{params.clean_fastq1} {output.clean_fastq1}
        mv {params.workdir}/{params.clean_fastq2} {output.clean_fastq2}
        """

rule soapnuke_filter_single:
    input:
        fastq = indir + "/{sample_id}.single.fq.gz"
    output:
        clean_fastq = outdir + "/{sample_id}.single.fq.gz"
    log:
        logdir + "/{sample_id}/soapnuke_filter_single.log"
    params:
        min_len = 15,
        low_qual = 0.2,
        max_n = 0.05,
        qual_type = 2,
        SOAPnuke = config.get('Procedure',{}).get('SOAPnuke') or 'SOAPnuke',
        workdir = lambda wildcards: f"{outdir}/{wildcards.sample_id}",
        clean_fastq = lambda wildcards: f"{wildcards.sample_id}.single.fq.gz"
    conda:
        "SOAPnuke.yaml"
    threads: 8
    shell:
        """
        mkdir -p {params.workdir}
        {params.SOAPnuke} filter \
            -l {params.min_len} \
            -q {params.low_qual} \
            -n {params.max_n} \
            -Q {params.qual_type} \
            -1 {input.fastq} \
            -C {params.clean_fastq} \
            -Q {params.qual_type} \
            -o {params.workdir} \
            -T {threads} \
            2>{log}
        mv {params.workdir}/{params.clean_fastq} {output.clean_fastq}
        """
rule soapnuke_filter_paired_result:
    input:
        clean_fastq1 = outdir + "/{sample_id}_1.fq.gz",
        clean_fastq2 = outdir + "/{sample_id}_2.fq.gz",
rule soapnuke_filter_single_result:
    input:
        clean_fastq = outdir + "/{sample_id}.single.fq.gz",
