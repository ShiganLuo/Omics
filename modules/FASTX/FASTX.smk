from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")

# not test
rule fastx_quality_filter_single:
    input:
        fastq = indir + "/{sample_id}.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}.fq.gz"),
    params:
        fastx_toolkit = config.get('Procedure',{}).get('fastx_toolkit') or 'fastq_quality_filter',
        q = config.get('Params',{}).get("fastx_toolkit", {}).get('q') or 10,
        p = config.get('Params',{}).get("fastx_toolkit", {}).get('p') or 100,
        Q = config.get('Params',{}).get("fastx_toolkit", {}).get('Q') or 33,
        outdir = outdir
    threads: 2
    conda:
        "FASTX.yaml"
    log:
        logdir + "/{sample_id}/fastx.txt"
    shell:
        """
        {params.fastx_toolkit} -Q {params.Q} -q {params.q} -p {params.p} \
            -i {input.fastq} -o {output.fastq} > {log} 2>&1
        """
