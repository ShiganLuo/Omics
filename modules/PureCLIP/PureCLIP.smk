from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output")
logdir = config.get("logdir", "log")
rule pureclip:
    input:
        bam = indir + "/{sample_id}.dedup.bam",
        fasta = config.get('genome',{}).get('fasta')
    output:
        sites = outdir + "/{sample_id}.pureclip.sites.bed",
        log = logdir + "/{sample_id}/pureclip.txt"
    params:
        pureclip = config.get('Procedure',{}).get('pureclip') or 'pureclip',
        ld = '--ld' if config.get('Params',{}).get('pureclip',{}).get('ld', True) else '',
        outdir = outdir
    threads: 2
    conda:
        "PureCLIP.yaml"
    log:
        log = logdir + "/{sample_id}/pureclip_run.txt"
    shell:
        """
        {params.pureclip} {params.ld} -nt {threads} \
            -i {input.bam} -bai {input.bam}.bai -g {input.fasta} \
            -o {output.sites} > {log} 2>&1
        """
