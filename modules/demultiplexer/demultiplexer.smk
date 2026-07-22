include: "../common/common.smk"

import time
import gzip

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
ranmer_len = config.get("Params", {}).get("demultiplexer", {}).get("ranmer_len", 10)
max_ham = config.get("Params", {}).get("demultiplexer", {}).get("max_ham", 1)

rule demultiplex_trim_dedup:
    input:
        r1 = indir + "/{sample_id}/{sample_id}_1.fq.gz",
        r2 = indir + "/{sample_id}/{sample_id}_2.fq.gz",
    output:
        r1 = outdir + "/{sample_id}/{sample_id}_1.fq.gz",
        r2 = outdir + "/{sample_id}/{sample_id}_2.fq.gz",
    log:
        logdir + "/{sample_id}/demultiplexer.log"
    threads: 1
    conda:
        "demultiplexer.yaml"
    params:
        ranmer_len = ranmer_len,
        max_ham = max_ham,
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("demultiplex_trim_dedup", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start demultiplex_trim_dedup for sample {wildcards.sample_id} at {current_time}")

            outdir_sample = os.path.dirname(str(output.r1))
            os.makedirs(outdir_sample, exist_ok=True)

            # Decompress gzipped FASTQ to temporary files (jla-trim requires uncompressed)
            tmp_r1 = os.path.join(outdir_sample, f"{wildcards.sample_id}_1.tmp.fq")
            tmp_r2 = os.path.join(outdir_sample, f"{wildcards.sample_id}_2.tmp.fq")

            with gzip.open(input.r1, 'rt') as f_in, open(tmp_r1, 'w') as f_out:
                f_out.write(f_in.read())
            with gzip.open(input.r2, 'rt') as f_in, open(tmp_r2, 'w') as f_out:
                f_out.write(f_in.read())

            rule_logger.info(f"Decompressed FASTQ files")

            # Run jla-trim: trims randommer+AG from read1 3' end, removes PCR duplicates
            cmd = f"jla-trim -r1 {tmp_r1} -r2 {tmp_r2} -r {params.ranmer_len} -maxHam {params.max_ham} >> {log_path} 2>&1"
            rule_logger.info(f"Running: {cmd}")
            shell(cmd)

            # jla-trim outputs {input}.trimmed.fastq - rename and compress
            trimmed_r1 = tmp_r1 + ".trimmed.fastq"
            trimmed_r2 = tmp_r2 + ".trimmed.fastq"

            shell(f"gzip -c {trimmed_r1} > {output.r1}")
            shell(f"gzip -c {trimmed_r2} > {output.r2}")

            # Cleanup temporary files
            for f in [tmp_r1, tmp_r2, trimmed_r1, trimmed_r2]:
                if os.path.exists(f):
                    os.remove(f)

            rule_logger.info(f"demultiplex_trim_dedup completed for sample {wildcards.sample_id}")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during demultiplex_trim_dedup for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during demultiplex_trim_dedup for sample {wildcards.sample_id}: {e}")
            raise e
