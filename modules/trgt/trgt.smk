from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai")
karyotype = config.get("Params", {}).get("trgt", {}).get("karyotype") or "XX"
repeat_bed = config.get("reference", {}).get("repeat_bed")

rule trgt_genotype:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        bai = indir + "/{sample_id}/{sample_id}.bam.bai",
        fasta = fasta,
        fai = fai,
        bed = repeat_bed
    output:
        vcf = outdir + "/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        bam = outdir + "/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam"
    log:
        logdir + "/{sample_id}/trgt_genotype.log"
    threads: 4
    conda:
        "trgt.yaml"
    params:
        trgt = config.get("Procedure", {}).get("trgt") or "trgt",
        karyotype = karyotype,
        prefix = outdir + "/genotype/{sample_id}/{sample_id}.trgt"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        {params.trgt} genotype \\
            --num-threads {threads} \\
            --genome {input.fasta} \\
            --repeats {input.bed} \\
            --karyotype {params.karyotype} \\
            --reads {input.bam} \\
            --output-prefix {params.prefix} \\
            > {log} 2>&1
        """

rule trgt_plot:
    input:
        vcf = outdir + "/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        bam = outdir + "/genotype/{sample_id}/{sample_id}.trgt.spanning.sorted.bam",
        fasta = fasta,
        fai = fai,
        bed = repeat_bed
    output:
        png = outdir + "/plot/{sample_id}/{sample_id}.trgt.repeat.png"
    log:
        logdir + "/{sample_id}/trgt_plot.log"
    conda:
        "trgt.yaml"
    params:
        trgt = config.get("Procedure", {}).get("trgt") or "trgt",
        repeat_id = config.get("Params", {}).get("trgt", {}).get("repeat_id") or "HTT"
    shell:
        """
        mkdir -p $(dirname {output.png})
        {params.trgt} plot \\
            --genome {input.fasta} \\
            --repeats {input.bed} \\
            --vcf {input.vcf} \\
            --spanning-reads {input.bam} \\
            --repeat-id {params.repeat_id} \\
            --image {output.png} \\
            > {log} 2>&1
        """

rule trgt_result:
    input:
        vcf = outdir + "/genotype/{sample_id}/{sample_id}.trgt.vcf.gz",
        png = outdir + "/plot/{sample_id}/{sample_id}.trgt.repeat.png"
