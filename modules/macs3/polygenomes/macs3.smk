include: "../../common/common.smk"
import shlex
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
ip_samples = config.get("ip_samples", [])
input_samples = config.get("input_samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})
bam_substring = config.get("bam_substring", "sorted_markdup.bam")
MODULE_DIR = os.path.join(config.get("ROOT_DIR", "."), "modules", "macs3")
CUTOFF_SCRIPT = os.path.join(MODULE_DIR, "bin", "plot_cutoff_analysis.py")

def get_macs3_input(wildcards):
    """
    Get treatment (IP) and optional control (Input) BAM files for MACS3.
    sample_ip_input_map: dict mapping IP sample_id -> input sample_id (or None)
    """
    bam_treatment = os.path.join(indir,f"{wildcards.sample_id}/{wildcards.sample_id}.{bam_substring}")
    
    # Check if there's a matched input control
    input_sample = sample_ip_input_map.get(wildcards.sample_id)
    if input_sample:
        bam_control = os.path.join(indir,f"{input_sample}/{input_sample}.{bam_substring}")
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
        narrow_peak = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.narrowPeak",
        narrow_xls = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.xls",
        narrow_cutoff = outdir + "/{genome}/{sample_id}/{sample_id}_cutoff_analysis.txt",
        broad_peak = outdir + "/{genome}/{sample_id}/{sample_id}_broad_peaks.broadPeak",
        broad_xls = outdir + "/{genome}/{sample_id}/{sample_id}_broad_peaks.xls",
        broad_cutoff = outdir + "/{genome}/{sample_id}/{sample_id}_broad_cutoff_analysis.txt"
    log:
        logdir + "/{genome}/{sample_id}/macs3.log"
    threads: 4
    conda:
        "../macs3.yaml"
    container:
        sif("../macs3.yaml")
    params:
        macs3 = config.get("Procedure", {}).get("macs3") or "macs3",
        name = lambda wildcards: wildcards.sample_id,
        bw = config.get("Params", {}).get("macs3", {}).get("bw") or 200,
        pvalue = config.get("Params", {}).get("macs3", {}).get("pvalue") or "1e-5",
        genome_size = lambda wildcards: config.get("Params", {}).get("macs3", {}).get("genome_size", {}).get(wildcards.genome, "mm") if isinstance(config.get("Params", {}).get("macs3", {}).get("genome_size"), dict) else (config.get("Params", {}).get("macs3", {}).get("genome_size") or "mm"),
        seed = 2346,
        broad_cutoff = config.get("Params",{}).get("macs3", {}).get("broad_cutoff") or "0.1"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("macs3_callpeak",log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start macs3 call peak for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.narrow_peak))
            script = os.path.join(sample_outdir, f"macs3_callpeak_{current_time}.sh")
            cmd = [
                params.macs3, "callpeak",
                "--bw", str(params.bw),
                "-p", str(params.pvalue),
                "-g", str(params.genome_size),
                "--outdir", sample_outdir,
                "--name", params.name,
                "--seed", str(params.seed),
                "-t", input.bam_treatment
            ]
            if hasattr(input, "bam_control") and input.bam_control:
                rule_logger.info(f"Using control BAM: {input.bam_control}")
                cmd += ["-c", input.bam_control]
            cmd.append("--cutoff-analysis")
            # ── Broad peak calling ──
            broad_name = f"{params.name}_broad"
            broad_cmd = [
                params.macs3, "callpeak",
                "--broad", "--broad-cutoff", str(params.broad_cutoff),
                "--bw", str(params.bw),
                "-p", str(params.pvalue),
                "-g", str(params.genome_size),
                "--outdir", sample_outdir,
                "--name", broad_name,
                "--seed", str(params.seed),
                "-t", input.bam_treatment
            ]
            if hasattr(input, "bam_control") and input.bam_control:
                broad_cmd += ["-c", input.bam_control]
            broad_cmd.append("--cutoff-analysis")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                f.write(" ".join(cmd) + "\n")
                f.write(" ".join(broad_cmd) + "\n")
                f.write(f'echo "macs3 call peak (narrow + broad) for sample {wildcards.sample_id} completed !"'+ "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path,"a") as f:
                f.write(f"Error occurred during macs3 call peak for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during macs3 call peak for sample {wildcards.sample_id}: {e}")
            raise e

rule macs3_cutoff_plot:
    """
    Plot combined MACS3 cutoff analysis curve for all IP samples.
    """
    input:
        cutoffs = expand(outdir + "/{genome}/{sample_id}/{sample_id}_cutoff_analysis.txt", sample_id=ip_samples, genome="{genome}"),
    output:
        plot = outdir + "/{genome}/cutoff_analysis.png",
    log:
        logdir + "/{genome}/../group/macs3/macs3_cutoff_plot.log"
    threads: 1
    conda:
        "../macs3.yaml"
    container:
        sif("../macs3.yaml")
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("macs3_cutoff_plot", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cutoff plot for {ip_samples} at {current_time}")

            cmd = ["python", CUTOFF_SCRIPT]
            for path, name in zip(input.cutoffs, ip_samples):
                cmd += ["--input-files", path, "--sample-names", name]
            cmd += ["--output", output.plot]

            plot_dir = os.path.dirname(str(output.plot))
            os.makedirs(plot_dir, exist_ok=True)
            script = os.path.join(plot_dir, f"plot_cutoff_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during cutoff plot: {e}\n")
            logger.error(f"Error during cutoff plot: {e}")
            raise e

rule macs3_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        narrow_peak = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.narrowPeak",
        narrow_xls = outdir + "/{genome}/{sample_id}/{sample_id}_peaks.xls",
        narrow_cutoff = outdir + "/{genome}/{sample_id}/{sample_id}_cutoff_analysis.txt",
        broad_peak = outdir + "/{genome}/{sample_id}/{sample_id}_broad_peaks.broadPeak",
        broad_xls = outdir + "/{genome}/{sample_id}/{sample_id}_broad_peaks.xls",
        broad_cutoff = outdir + "/{genome}/{sample_id}/{sample_id}_broad_cutoff_analysis.txt",
        plot = outdir + "/{genome}/cutoff_analysis.png"
