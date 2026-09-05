include: "../../common/common.smk"
import shlex

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
peaks_indir = config.get("peaks_indir", "peaks")
bam_indir = config.get("bam_indir", "")
samples = config.get("samples", [])
sample_ip_input_map = config.get("sample_ip_input_map", {})
MODULE_DIR = os.path.join(config.get("ROOT_DIR", "."), "modules", "peak_te_overlap")
CHART_SCRIPT = os.path.join(MODULE_DIR, "bin", "plot_te_overlap.py")


def get_input_bam(wildcards):
    """Get the matched Input BAM for an IP sample."""
    input_sample = sample_ip_input_map.get(wildcards.sample_id)
    if input_sample and bam_indir:
        return bam_indir + f"/{input_sample}/{input_sample}.sorted_markdup.bam"
    return ""


rule peak_te_overlap:
    """
    Find overlaps between narrowPeak and TE annotation using bedtools.
    Produces per-sample overlap BED + overlap count TSV per TE class.
    Also counts reads from IP and Input BAMs in overlap regions.
    """
    input:
        peak = peaks_indir + "/{sample_id}/{sample_id}_peaks.narrowPeak",
        ip_bam = bam_indir + "/{sample_id}/{sample_id}.sorted_markdup.bam" if bam_indir else [],
        input_bam = get_input_bam,
    output:
        overlap_bed = outdir + "/{genome}/{sample_id}/{sample_id}_peak_te_overlap.bed",
        subfamily_tsv = outdir + "/{genome}/{sample_id}/{sample_id}_te_subfamily_overlap.tsv",
    log:
        logdir + "/{genome}/{sample_id}/peak_te_overlap.log"
    threads: 1
    conda:
        "../peak_te_overlap.yaml"
    container:
        sif("../peak_te_overlap.yaml")
    params:
        bedtools = config.get("Procedure", {}).get("bedtools") or "bedtools",
        te_gtf = lambda wildcards: config['genomes'][wildcards.genome]['te_gtf'],
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("peak_te_overlap", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start peak-TE overlap for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.overlap_bed))
            os.makedirs(sample_outdir, exist_ok=True)

            cmd = [
                "python", os.path.join(MODULE_DIR, "bin", "intersect_te.py"),
                "--peak", input.peak,
                "--te-gtf", params.te_gtf,
                "--bedtools", params.bedtools,
                "--overlap-bed", output.overlap_bed,
                "--subfamily-tsv", output.subfamily_tsv,
                "--sample-id", wildcards.sample_id,
            ]
            if input.ip_bam:
                cmd += ["--ip-bam", input.ip_bam]
            if input.input_bam:
                cmd += ["--input-bam", input.input_bam]

            script = os.path.join(sample_outdir, f"peak_te_overlap_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write(f'echo "peak-TE overlap for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during peak-TE overlap for {wildcards.sample_id}: {e}\n")
            logger.error(f"Error during peak-TE overlap for {wildcards.sample_id}: {e}")
            raise e


rule peak_centric_overlap:
    """
    Peak-centric TE overlap analysis.
    For each peak: TE count, TE types, TE coverage fraction.
    Useful for experimental validation — identify which peaks contain TEs.
    """
    input:
        peak = peaks_indir + "/{sample_id}/{sample_id}_peaks.narrowPeak",
    output:
        tsv = outdir + "/{genome}/{sample_id}/{sample_id}_peak_centric_te.tsv",
    log:
        logdir + "/{genome}/{sample_id}/peak_centric_overlap.log"
    threads: 1
    conda:
        "../peak_te_overlap.yaml"
    container:
        sif("../peak_te_overlap.yaml")
    params:
        bedtools = config.get("Procedure", {}).get("bedtools") or "bedtools",
        te_gtf = lambda wildcards: config['genomes'][wildcards.genome]['te_gtf'],
        script = os.path.join(MODULE_DIR, "bin", "peak_centric_overlap.py"),
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("peak_centric_overlap", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start peak-centric TE overlap for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)

            cmd = [
                "python3", params.script,
                "--peak", input.peak,
                "--te-gtf", params.te_gtf,
                "--output", output.tsv,
                "--bedtools", params.bedtools,
                "--sample-id", wildcards.sample_id,
            ]

            script = os.path.join(sample_outdir, f"peak_centric_overlap_{current_time}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write(f'echo "peak-centric TE overlap for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error during peak-centric TE overlap for {wildcards.sample_id}: {e}\n")
            logger.error(f"Error during peak-centric TE overlap for {wildcards.sample_id}: {e}")
            raise e


rule peak_te_overlap_fig:
    """
    Aggregate per-sample TE subfamily overlap data and generate enrichment plots.
    Default: separate plot per IP:Input pair.
    With Params.peak_te_overlap.combine=true: single combined plot.
    """
    input:
        tsvs = expand(outdir + "/{genome}/{sid}/{sid}_te_subfamily_overlap.tsv", sid=samples, genome="{genome}"),
    output:
        chart = outdir + "/{genome}/te_subfamily_overlap_combined.png",
        combined_tsv = outdir + "/{genome}/te_subfamily_overlap_combined.tsv",
        done = touch(outdir + "/{genome}/.peak_te_overlap_fig.done"),
    log:
        logdir + "/{genome}/../group/peak_te_overlap/peak_te_overlap.log"
    threads: 1
    conda:
        "../peak_te_overlap.yaml"
    container:
        sif("../peak_te_overlap.yaml")
    params:
        sample_ip_input_map = config.get("sample_ip_input_map", {}),
        method = config.get("Params", {}).get("peak_te_overlap", {}).get("method", "count"),
        sort_by = config.get("Params", {}).get("peak_te_overlap", {}).get("sort_by", "te_length"),
        top_n = config.get("Params", {}).get("peak_te_overlap", {}).get("top_n", 30),
        combine = config.get("Params", {}).get("peak_te_overlap", {}).get("combine", False),
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("peak_te_overlap_chart", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start peak-TE overlap chart generation at {current_time}")

            chart_dir = os.path.dirname(str(output.chart))
            os.makedirs(chart_dir, exist_ok=True)

            if params.combine:
                output_target = output.chart
            else:
                output_target = chart_dir

            cmd = [
                "python", CHART_SCRIPT,
                "--tsvs", ",".join(input.tsvs),
                "--output", output_target,
                "--combined-tsv", output.combined_tsv,
                "--method", params.method,
                "--sort-by", params.sort_by,
                "--top-n", str(params.top_n),
            ]
            if params.combine:
                cmd += ["--combine"]
            if params.sample_ip_input_map:
                pairs = ";".join(f"{ip}:{inp}" for ip, inp in params.sample_ip_input_map.items())
                cmd += ["--ip-input-pairs", pairs]

            script_path = os.path.join(chart_dir, f"peak_te_overlap_chart_{current_time}.sh")
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
            shell(f"bash {script_path} >> {log_path} 2>&1")

            # For separate mode: create the expected single chart as symlink to first pair
            if not params.combine and params.sample_ip_input_map:
                first_ip = list(params.sample_ip_input_map.keys())[0]
                first_inp = params.sample_ip_input_map[first_ip]
                first_chart = os.path.join(chart_dir, f"{first_ip}_vs_{first_inp}_enrichment.png")
                if os.path.exists(first_chart) and not os.path.exists(output.chart):
                    os.symlink(first_chart, output.chart)
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
        chart = outdir + "/{genome}/te_subfamily_overlap_combined.png",
        peak_centric = expand(outdir + "/{genome}/{sid}/{sid}_peak_centric_te.tsv", sid=samples, genome="{genome}"),
