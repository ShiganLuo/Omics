include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})

def get_macs3_input(wildcards):
    """
    Get treatment (IP) and optional control (Input) BAM files for MACS3.
    sample_ip_input_map: dict mapping IP sample_id -> input sample_id (or None)
    """
    sample = wildcards.sample_id
    bam_treatment = os.path.join(indir,f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
    
    # Check if there's a matched input control
    input_sample = sample_ip_input_map.get(sample)
    if input_sample:
        bam_control = os.path.join(indir,f"{input_sample}/{input_sample}.bam")
        return {
            "bam_treatment": bam_treatment,
            "bam_control": bam_control
        }
    
    return {"bam_treatment": bam_treatment}

rule macs3_callpeak:
    """
    MACS3 peak calling for ChIP-seq/DIP-seq data.
    Supports both with-control and without-control modes.
    """
    input:
        unpack(get_macs3_input)
    output:
        peak = outdir + "/{sample_id}/{sample_id}_peaks.narrowPeak",
        xls = outdir + "/{sample_id}/{sample_id}_peaks.xls"
    log:
        logdir + "/{sample_id}/macs3.log"
    threads: 4
    conda:
        "macs3.yaml"
    params:
        macs3 = config.get("Procedure", {}).get("macs3") or "macs3",
        name = lambda wildcards: wildcards.sample_id,
        bw = config.get("Params", {}).get("macs3", {}).get("bw") or 200,
        pvalue = config.get("Params", {}).get("macs3", {}).get("pvalue") or "1e-5",
        genome_size = config.get("Params", {}).get("macs3", {}).get("genome_size") or "mm",
        seed = 2346
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            logger = setup_logger("macs3_callpeak",log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start macs3 call peak for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.peak))
            script = os.path.join(sample_outdir,f"/macs3_callpeak_{current_time}.sh")
            cmd = [
                params.macs3, callpeak,
                "--bw", params.bw,
                "-p",params.pvalue,
                "-g", params.genome_size,
                "--outdir", sample_outdir,
                "--name", params.name,
                "--seed", params.seed,
                "-t", input.bam_treatment
            ]
            if "bam_control" in input:
                cmd += ["-c", input.bam_control]
            success_echo = f'echo "macs3 call peak for sample {wildcards.sample_id}"'
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(success_echo + "\n")
            shell(f"bash {script} > {log_path} 2>&1")

        except Exception:
            with open(log_path,"a") as f:
                f.write(f"Error occurred during macs3 call peak for sample {wildcards.sample_id}: {e}\n")

rule macs3_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        peak = outdir + "/{sample_id}/{sample_id}_peaks.narrowPeak",
        xls = outdir + "/{sample_id}/{sample_id}_peaks.xls"
