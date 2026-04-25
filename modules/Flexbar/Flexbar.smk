from snakemake.logging import logger
outdir = config.get("outdir", "output")
indir = config.get("indir", "output/raw_fastq")
logdir = config.get("logdir", "log")

# need test
rule flexbar_demultiplex:
    input:
        fastq = indir + "/{sample_id}.fq.gz",
        barcodes = config.get('Params',{}).get('flexbar',{}).get('barcode_file','barcodes.fasta')
    output:
        demux = temp(outdir + "/{sample_id}.demux.fq.gz"),
        log = logdir + "/{sample_id}/flexbar_demux.txt"
    params:
        flexbar = config.get('Procedure',{}).get('flexbar') or 'flexbar',
        min_overlap = config.get('Params',{}).get('flexbar',{}).get('barcode_min_overlap', 1),
        allow_mismatch = config.get('Params',{}).get('flexbar',{}).get('barcode_allow_mismatch', 0),
        outdir = outdir
    threads: 4
    conda:
        "Flexbar.yaml"
    log:
        log = logdir + "/{sample_id}/flexbar_demux_run.txt"
    shell:
        """
        {params.flexbar} -r {input.fastq} -b {input.barcodes} \
            -n {params.min_overlap} -m {params.allow_mismatch} \
            -t {params.outdir}/{wildcards.sample_id}.demux \
            > {log} 2>&1
        pigz -c {params.outdir}/{wildcards.sample_id}.demux_1.fastq > {output.demux}
        """
