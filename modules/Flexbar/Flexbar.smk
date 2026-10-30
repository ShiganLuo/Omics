include: "../common/common.smk"

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
    container:
        sif("Flexbar.yaml")
    log:
        log = logdir + "/{sample_id}/flexbar_demux_run.txt"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("flexbar_demultiplex", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start flexbar_demultiplex for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.join(params.outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"flexbar_demultiplex_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"flexbar_demultiplex.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"{params.flexbar} -r {input.fastq} -b {input.barcodes} \\\n")
                f.write(f"    -n {params.min_overlap} -m {params.allow_mismatch} \\\n")
                f.write(f"    -t {params.outdir}/{wildcards.sample_id}.demux \\\n")
                f.write(f"pigz -c {params.outdir}/{wildcards.sample_id}.demux_1.fastq > {output.demux}\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
            rule_logger.info(f"flexbar_demultiplex for sample {wildcards.sample_id} completed")
        except Exception as e:
            rule_logger.error(f"Error occurred: {e}")
            raise RuntimeError(f"Error occurred: {e}, you can check the log file {log_path} for more details.")
