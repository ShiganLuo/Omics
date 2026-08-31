include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
peaks_indir = config.get("peaks_indir", "peaks")
samples = config.get("samples", [])
te_gtf = config.get("genome", {}).get("te_gtf")
MODULE_DIR = os.path.join(config.get("ROOT_DIR", "."), "modules", "peak_te_overlap")
CHART_SCRIPT = os.path.join(MODULE_DIR, "bin", "plot_te_overlap.py")

if not te_gtf:
    raise ValueError(
        "peak_te_overlap module requires 'genome.te_gtf' in config. "
        "Please provide a valid TE annotation GTF path."
    )


rule peak_te_overlap:
    """
    Find overlaps between narrowPeak and TE annotation using bedtools.
    Produces per-sample overlap BED + overlap count TSV per TE class.
    """
    input:
        peak = peaks_indir + "/{sample_id}/{sample_id}_peaks.narrowPeak",
    output:
        overlap_bed = outdir + "/{sample_id}/{sample_id}_peak_te_overlap.bed",
        overlap_tsv = outdir + "/{sample_id}/{sample_id}_te_overlap_counts.tsv",
    log:
        logdir + "/{sample_id}/peak_te_overlap.log"
    threads: 1
    conda:
        "peak_te_overlap.yaml"
    container:
        sif("peak_te_overlap.yaml")
    params:
        bedtools = config.get("Procedure", {}).get("bedtools") or "bedtools",
        te_gtf = te_gtf,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("peak_te_overlap", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start peak-TE overlap for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.overlap_bed))
            os.makedirs(sample_outdir, exist_ok=True)

            script = os.path.join(sample_outdir, f"peak_te_overlap_{current_time}.sh")
            cmd = [
                "python", os.path.join(MODULE_DIR, "bin", "intersect_te.py"),
                "--peak", input.peak,
                "--te-gtf", params.te_gtf,
                "--bedtools", params.bedtools,
                "--overlap-bed", output.overlap_bed,
                "--overlap-tsv", output.overlap_tsv,
                "--sample-id", wildcards.sample_id,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "peak-TE overlap for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during peak-TE overlap for {wildcards.sample_id}: {e}\n")
            logger.error(f"Error during peak-TE overlap for {wildcards.sample_id}: {e}")
            raise e


rule peak_te_overlap_chart:
    """
    Aggregate per-sample TE overlap counts and generate a grouped bar chart
    comparing IP vs Input overlap ratios across TE families.
    """
    input:
        tsvs = expand(outdir + "/{sid}/{sid}_te_overlap_counts.tsv", sid=samples),
    output:
        chart = outdir + "/te_family_overlap.png",
        combined_tsv = outdir + "/te_family_overlap.tsv",
    log:
        logdir + "/peak_te_overlap_chart.log"
    threads: 1
    conda:
        "peak_te_overlap.yaml"
    container:
        sif("peak_te_overlap.yaml")
    params:
        sample_ip_input_map = config.get("sample_ip_input_map", {}),
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("peak_te_overlap_chart", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start peak-TE overlap chart generation at {current_time}")
            os.makedirs(os.path.dirname(str(output.chart)), exist_ok=True)

            cmd = [
                "python", CHART_SCRIPT,
                "--tsvs", ",".join(input.tsvs),
                "--output", output.chart,
                "--combined-tsv", output.combined_tsv,
            ]
            if params.sample_ip_input_map:
                pairs = ";".join(f"{ip}:{inp}" for ip, inp in params.sample_ip_input_map.items())
                cmd += ["--ip-input-pairs", pairs]

            script_path = os.path.join(outdir, f"peak_te_overlap_chart_{current_time}.sh")
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "peak-TE overlap chart completed at {current_time}"\n')
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during peak-TE overlap chart: {e}\n")
            logger.error(f"Error during peak-TE overlap chart: {e}")
            raise e


rule peak_te_overlap_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        chart = outdir + "/te_family_overlap.png",
        combined_tsv = outdir + "/te_family_overlap.tsv",
