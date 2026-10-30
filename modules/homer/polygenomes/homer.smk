include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])

rule homer_annotatepeaks:
    """
    HOMER annotatePeaks.pl: annotate peaks with genomic features
    (promoter, intron, intergenic, etc.) and nearest gene information.
    """
    input:
        peak = indir + "/{sample_id}/{sample_id}_peaks.narrowPeak",
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta'],
        gtf = lambda wildcards: config['genome'][wildcards.genome]['gtf']
    output:
        annotation = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.annotatePeaks.txt"
    log:
        logdir + "/{genome}/{sample_id}/homer_annotatepeaks.log"
    threads: 1
    conda:
        "../homer.yaml"
    container:
        sif("../homer.yaml")
    params:
        annotatePeaks = config.get("Procedure", {}).get("annotatePeaks") or "annotatePeaks.pl"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
rule_logger = setup_logger("homer_annotatepeaks", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start HOMER annotatePeaks for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.annotation))
            script = os.path.join(RECEIPT_DIR, f"homer_annotatepeaks_{wildcards.genome}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"homer_annotatepeaks.sh")
            cmd = [
                params.annotatePeaks,
                str(input.peak),
                str(input.fasta),
                "-gtf", str(input.gtf),
                ">", str(output.annotation)
            ]
            success_echo = f'echo "HOMER annotatePeaks for sample {wildcards.sample_id} successfully completed !"'
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(success_echo + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} > {log_path} 2>&1")
        except Exception as e:
            logger.error(f"Error occurred during HOMER annotatePeaks for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred: 2, you can check the log file for more details.")

rule homer_annotatepeaks_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        annotation = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.annotatePeaks.txt"
