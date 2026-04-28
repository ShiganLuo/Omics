from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")
# neet for test
rule trimmomatic_Paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fq.gz",
        fastq2 = indir + "/{sample_id}_2.fq.gz"
    output:
        fastq1 = temp(outdir + "/{sample_id}_1.fq.gz"),
        fastq2 = temp(outdir + "/{sample_id}_2.fq.gz"),
        report = outdir + "/{sample_id}/trimmomatic_report.txt"
    params:
        outdir = outdir,
        trimmomatic = config.get('Procedure',{}).get('trimmomatic') or 'trimmomatic',
        adapter = config.get('Params',{}).get("trimmomatic", {}).get('adapter_pe'),
        cmd = lambda wildcards: (
            f"java -jar {config.get('Procedure',{}).get('trimmomatic')}"
            if str(config.get('Procedure',{}).get('trimmomatic')).endswith(".jar")
            else config.get('Procedure',{}).get('trimmomatic') or "trimmomatic"
        )
    threads: 6
    conda:
        "trimmomatic.yaml"
    log:
        log = logdir + "/{sample_id}/trimmomatic.txt"
    shell:
        """
        {params.cmd} PE -threads {threads} \
            {input.fastq1} {input.fastq2} \
            -summary {output.report} \
            {output.fastq1} {params.outdir}/{wildcards.sample_id}_1.unpaired.fq.gz \
            {output.fastq2} {params.outdir}/{wildcards.sample_id}_2.unpaired.fq.gz \
            ILLUMINACLIP:{params.adapter}:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:80 \
            2> {log.log}
        """

rule trimmomatic_Single:
    input:
        fastq = indir + "/{sample_id}.single.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}.single.fq.gz"),
        report = outdir + "/{sample_id}/trimmomatic_report.txt"
    params:
        outdir = outdir,
        trimmomatic = config.get('Procedure',{}).get('trimmomatic') or 'trimmomatic',
        adapter = config.get('Params',{}).get("trimmomatic", {}).get('adapter_se'),
        cmd = lambda wildcards: (
            f"java -jar {config.get('Procedure',{}).get('trimmomatic')}"
            if str(config.get('Procedure',{}).get('trimmomatic')).endswith(".jar")
            else config.get('Procedure',{}).get('trimmomatic') or "trimmomatic"
        )
    threads: 6
    conda:
        "trimmomatic.yaml"
    log:
        log = logdir + "/{sample_id}/trimmomatic.txt"
    shell:
        """
        {params.cmd} SE -threads {threads} \
            {input.fastq} {output.fastq} \
            -summary {output.report} \
            ILLUMINACLIP:{params.adapter}:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:80 \
            2> {log.log}
        """