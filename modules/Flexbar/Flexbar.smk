include: "../common/common.smk"

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
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("flexbar_demultiplex", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start flexbar_demultiplex for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(params.outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"flexbar_demultiplex_{current_time}.sh")
            with open(script, "w") as f:
                f.write(f"{params.flexbar} -r {input.fastq} -b {input.barcodes} \\\n")
                f.write(f"    -n {params.min_overlap} -m {params.allow_mismatch} \\\n")
                f.write(f"    -t {params.outdir}/{wildcards.sample_id}.demux \\\n")
                f.write(f"    > {log_path} 2>&1\n")
                f.write(f"pigz -c {params.outdir}/{wildcards.sample_id}.demux_1.fastq > {output.demux}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
            rule_logger.info(f"flexbar_demultiplex for sample {wildcards.sample_id} completed")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"flexbar_demultiplex failed for sample {wildcards.sample_id}: {e}\n")
            raise e
