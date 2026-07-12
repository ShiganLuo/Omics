include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])

# Genome references
fasta = config.get("genome", {}).get("fasta", "")
gtf = config.get("genome", {}).get("gtf", "")

if not fasta:
    raise ValueError("homer module requires 'genome.fasta' in config, got None or empty.")
if not gtf:
    raise ValueError(
        "homer module requires 'genome.gtf' in config for peak annotation. "
        "Please provide a valid GTF/GFF path (e.g. gencode.vM35.annotation.gtf)."
    )

rule homer_annotatepeaks:
    """
    HOMER annotatePeaks.pl: annotate peaks with genomic features
    (promoter, intron, intergenic, etc.) and nearest gene information.
    """
    input:
        peak = indir + "/{sample_id}/{sample_id}_peaks.narrowPeak",
        fasta = fasta,
        gtf = gtf
    output:
        annotation = outdir + "/{sample_id}/{sample_id}_peaks.annotatePeaks.txt"
    log:
        logdir + "/{sample_id}/homer_annotatepeaks.log"
    threads: 1
    conda:
        "homer.yaml"
    params:
        annotatePeaks = config.get("Procedure", {}).get("annotatePeaks") or "annotatePeaks.pl"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("homer_annotatepeaks", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start HOMER annotatePeaks for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.annotation))
            script = os.path.join(sample_outdir, f"homer_annotatepeaks_{current_time}.sh")
            cmd = [
                params.annotatePeaks,
                str(input.peak),
                str(input.fasta),
                "-gtf", str(input.gtf),
                ">", str(output.annotation)
            ]
            success_echo = f'echo "HOMER annotatePeaks for sample {wildcards.sample_id} successfully completed !"'
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(success_echo + "\n")
            shell(f"bash {script} > {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during HOMER annotatePeaks for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during HOMER annotatePeaks for sample {wildcards.sample_id}: {e}")
            raise e

rule homer_annotatepeaks_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        annotation = outdir + "/{sample_id}/{sample_id}_peaks.annotatePeaks.txt"
