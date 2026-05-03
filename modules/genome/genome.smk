from snakemake.logging import logger
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
rule chromosome_sizes:
    input:
        fasta = config.get('genome',{}).get('fasta')
    output:
        chrom_sizes = outdir + "/genome/chrom.sizes"
    log:
        logdir + "/genome/chromosome_sizes.log"
    threads: 1
    conda:
        "genome.yaml"
    params:
        samtools = config.get('Procedure',{}).get('samtools') or 'samtools'
    shell:
        """
        {params.samtools} faidx {input.fasta} 2> {log}
        cut -f1,2 {input.fasta}.fai > {output.chrom_sizes}
        """
