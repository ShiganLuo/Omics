from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")

rule fastx_quality_filter:
    input:
        fastq = indir + "/{sample_id}.fq.gz"
    output:
        fastq = temp(outdir + "/{sample_id}.filtered.fq.gz"),
        report = logdir + "/{sample_id}/fastx_statistics.txt"
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
        log = logdir + "/{sample_id}/fastx.txt"
    shell:
        """
        mkdir -p {params.outdir}
        {params.fastx_toolkit} -Q {params.Q} -q {params.q} -p {params.p} \
            -i {input.fastq} -o {output.fastq} > {log.log} 2>&1
        echo "fastx_toolkit finished for {wildcards.sample_id}" > {output.report}
        """
