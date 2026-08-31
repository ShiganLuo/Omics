include: "../common/common.smk"
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
bigwig_dir = config.get("bigwig_dir", "")

# --- Regions resolution ---
# Three modes:
#   1. "peaks" (default): per-sample MACS3 narrowPeak → reference-point, center
#   2. path/to/regions.bed: user-provided BED → auto mode (reference-point or scale-regions)
#   3. "tss": auto-generate TSS BED from GTF → reference-point, TSS
regions_cfg = config.get("Params", {}).get("computeMatrix", {}).get("regions", "peaks")
gtf = config.get("genome", {}).get("gtf")
MODULE_DIR = os.path.join(config.get("ROOT_DIR", "."), "modules", "deeptools_heatmap")
GENERATE_TSS_SCRIPT = os.path.join(MODULE_DIR, "bin", "generate_tss_bed.py")


def get_bigwig_for_sample(wildcards):
    """Return BigWig file path for a given sample."""
    return os.path.join(bigwig_dir, wildcards.sample_id, wildcards.sample_id + ".bigwig")


def get_regions(wildcards):
    """Resolve BED regions file based on config."""
    if regions_cfg == "peaks":
        return os.path.join(indir, wildcards.sample_id, wildcards.sample_id + "_peaks.narrowPeak")
    elif regions_cfg == "tss":
        if not gtf:
            raise ValueError("computeMatrix regions='tss' requires genome.gtf in config")
        return os.path.join(outdir, "_tss_regions.bed")
    else:
        if not os.path.isfile(regions_cfg):
            raise ValueError(f"computeMatrix regions BED not found: {regions_cfg}")
        return regions_cfg


# Pre-generate TSS BED from GTF (runs once, shared across all samples)
# Only triggered when regions='tss' via get_regions input function dependency.
rule generate_tss_bed:
    input:
        gtf = gtf or "/dev/null",
    output:
        bed = outdir + "/_tss_regions.bed",
    log:
        logdir + "/generate_tss_bed.log"
    threads: 1
    conda:
        "deeptools_heatmap.yaml"
    container:
        sif("deeptools_heatmap.yaml")
    params:
        script = GENERATE_TSS_SCRIPT,
        flank = config.get("Params", {}).get("computeMatrix", {}).get("tss_flank") or 1000,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("generate_tss_bed", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start generating TSS BED from GTF at {current_time}")
            os.makedirs(os.path.dirname(str(output.bed)), exist_ok=True)
            script_path = os.path.join(outdir, f"generate_tss_bed_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--gtf", input.gtf,
                "--output", output.bed,
                "--flank", str(params.flank),
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "TSS BED generation completed at {current_time}"\n')
            shell(f"bash {script_path} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error generating TSS BED: {e}\n")
            logger.error(f"Error generating TSS BED: {e}")
            raise e


rule computeMatrix:
    """
    Compute signal matrix around specified regions for enrichment heatmap.

    Regions can be:
      - MACS3 peaks (default): reference-point mode, center
      - User BED file: auto-detected mode (reference-point or scale-regions)
      - TSS regions: auto-generated from GTF, reference-point mode, TSS
    """
    input:
        bigwig = get_bigwig_for_sample,
        regions = get_regions,
    output:
        matrix = outdir + "/{sample_id}/{sample_id}_matrix.gz",
    log:
        logdir + "/{sample_id}/computeMatrix.log"
    threads: 4
    conda:
        "deeptools_heatmap.yaml"
    container:
        sif("deeptools_heatmap.yaml")
    params:
        computeMatrix = config.get("Procedure", {}).get("computeMatrix") or "computeMatrix",
        mode = config.get("Params", {}).get("computeMatrix", {}).get("mode") or "reference-point",
        referencePoint = config.get("Params", {}).get("computeMatrix", {}).get("referencePoint") or "center",
        before = config.get("Params", {}).get("computeMatrix", {}).get("before") or 3000,
        after = config.get("Params", {}).get("computeMatrix", {}).get("after") or 3000,
        upstream = config.get("Params", {}).get("computeMatrix", {}).get("upstream") or 3000,
        downstream = config.get("Params", {}).get("computeMatrix", {}).get("downstream") or 3000,
        bodyLength = config.get("Params", {}).get("computeMatrix", {}).get("bodyLength") or 5000,
        binSize = config.get("Params", {}).get("computeMatrix", {}).get("binSize") or 10,
        sortUsing = config.get("Params", {}).get("computeMatrix", {}).get("sortUsing") or "mean",
        missingDataAsZero = config.get("Params", {}).get("computeMatrix", {}).get("missingDataAsZero") or True,
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("computeMatrix", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start computeMatrix for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.matrix))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"computeMatrix_{current_time}.sh")

            # Auto-select mode based on regions source
            mode = params.mode
            if regions_cfg == "peaks":
                mode = "reference-point"
                ref_point = "center"
            elif regions_cfg == "tss":
                mode = "reference-point"
                ref_point = "TSS"
            else:
                ref_point = params.referencePoint

            cmd = [
                params.computeMatrix, mode,
                "--binSize", str(params.binSize),
                "--sortUsing", params.sortUsing,
                "--numberOfProcessors", str(threads),
                "-R", input.regions,
                "-S", input.bigwig,
                "-o", output.matrix,
            ]

            if mode == "reference-point":
                cmd += [
                    "--referencePoint", ref_point,
                    "--beforeRegionStartLength", str(params.before),
                    "--afterRegionStartLength", str(params.after),
                ]
            elif mode == "scale-regions":
                cmd += [
                    "--regionBodyLength", str(params.bodyLength),
                    "--upstream", str(params.upstream),
                    "--downstream", str(params.downstream),
                ]

            if params.missingDataAsZero:
                cmd.append("--missingDataAsZero")

            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "computeMatrix for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during computeMatrix for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during computeMatrix for sample {wildcards.sample_id}: {e}")
            raise e


rule plotHeatmap:
    """
    Plot enrichment heatmap from computeMatrix output.
    Produces heatmap PNG and optional color bar / profile plot.
    """
    input:
        matrix = outdir + "/{sample_id}/{sample_id}_matrix.gz",
    output:
        heatmap = outdir + "/{sample_id}/{sample_id}_heatmap.png",
    log:
        logdir + "/{sample_id}/plotHeatmap.log"
    threads: 1
    conda:
        "deeptools_heatmap.yaml"
    container:
        sif("deeptools_heatmap.yaml")
    params:
        plotHeatmap = config.get("Procedure", {}).get("plotHeatmap") or "plotHeatmap",
        colorMap = config.get("Params", {}).get("plotHeatmap", {}).get("colorMap") or "YlOrRd",
        heatmapHeight = config.get("Params", {}).get("plotHeatmap", {}).get("heatmapHeight") or 15,
        heatmapWidth = config.get("Params", {}).get("plotHeatmap", {}).get("heatmapWidth") or 8,
        whatToShow = config.get("Params", {}).get("plotHeatmap", {}).get("whatToShow") or "heatmap, colorbar, metagene",
        plotTitle = "",
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("plotHeatmap", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start plotHeatmap for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.heatmap))
            os.makedirs(sample_outdir, exist_ok=True)
            plot_title = params.plotTitle or f"{wildcards.sample_id} Peak Enrichment Heatmap"
            script = os.path.join(sample_outdir, f"plotHeatmap_{current_time}.sh")
            cmd = [
                params.plotHeatmap,
                "-m", input.matrix,
                "-o", output.heatmap,
                "--colorMap", params.colorMap,
                "--heatmapHeight", str(params.heatmapHeight),
                "--heatmapWidth", str(params.heatmapWidth),
                "--whatToShow", params.whatToShow,
                "--plotTitle", plot_title,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "plotHeatmap for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during plotHeatmap for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during plotHeatmap for sample {wildcards.sample_id}: {e}")
            raise e


rule deeptools_heatmap_result:
    """
    Result aggregation rule for subworkflow use rule import.
    """
    input:
        matrix = outdir + "/{sample_id}/{sample_id}_matrix.gz",
        heatmap = outdir + "/{sample_id}/{sample_id}_heatmap.png",
