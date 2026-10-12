"""Decoil ecDNA reconstruction module for PacBio HiFi long reads.

Reconstructs extrachromosomal circular DNA (eccDNA/ecDNA) from long-read
sequencing BAM files using Decoil's ``decoil reconstruct`` sub-command
(Giurgiu et al. 2024, Genome Research, DOI: 10.1101/gr.279123.124).

This module is **reconstruction-only**. SV calling (Sniffles, pbsv, ...) and
coverage track generation (deeptools bamCoverage) must be supplied as
upstream inputs — they are NOT performed by this module. The recommended
upstream chain is::

    pbsv  ->  {vcf_dir}/{sample_id}/{sample_id}.sv.vcf.gz
    igv.wig -> {bw_dir}/{sample_id}/{sample_id}.bigwig

Decoil needs only three files to reconstruct ecDNA:
- aligned, coordinate-sorted BAM (+ BAI)
- an SV VCF from sniffles1 / sniffles2 / lumpy / multivcf
- a coverage bigwig (genome-wide read depth)

Inputs:
    - PacBio HiFi BAM (markdup or coordinate-sorted) and BAI index
    - Reference FASTA
    - Gene annotation GTF
    - External SV VCF (``vcf_dir`` + ``vcf_suffix``)
    - External coverage bigwig (``bw_dir`` + ``.bigwig``)

Outputs:
    - reconstruct.bed: all genomic fragments in order composing all reconstructions
    - reconstruct.ecDNA.bed: fragments in order composing ecDNA-labeled cycles
    - reconstruct.ecDNA.filtered.bed: ecDNA reconstructions passing --filter-score
    - summary.txt: one row per reconstructed cycle with topology and proportion

Citation:
    Giurgiu et al., "Reconstructing extrachromosomal DNA structural
    heterogeneity from long-read sequencing data using Decoil",
    Genome Research, 2024. doi:10.1101/gr.279123.124
"""

include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""
vcf_dir = config.get("vcf_dir")
bw_dir = config.get("bw_dir") or outdir
vcf_suffix = config.get("vcf_suffix") or ".sv.vcf.gz"
bw_suffix = config.get("bw_suffix") or ".bigwig"

# Decoil rule params — per-rule subfields under Params.decoil
# (reconstruct / plot / compare / igv), shared params at the parent level
# (e.g. formats). Rule-level keys fall back to the parent field.
decoil_params = config.get("Params", {}).get("decoil", {})


def reconstruct_param(key: str, default=None):
    """Resolve Params.decoil.reconstruct.<key> → Params.decoil.<key> → default.

    Parameters
    ----------
    key : str
        Reconstruction parameter name.
    default : object, optional
        Value when neither the rule subfield nor the parent field defines it.

    Returns
    -------
    object
        The resolved config value.
    """
    return decoil_params.get("reconstruct", {}).get(key, decoil_params.get(key, default))


PLOT_FORMATS = (decoil_params.get("plot", {}).get("formats")
                or decoil_params.get("formats") or ["png"])
COMPARE_FORMATS = (decoil_params.get("compare", {}).get("formats")
                   or decoil_params.get("formats") or ["png"])
IGV_LOCUS = decoil_params.get("igv", {}).get("locus", "")

# pbsv→decoil VCF converter — bin/ script located via ROOT_DIR (project convention)
PBSV2DECOIL_SCRIPT = os.path.join(ROOT_DIR, "modules", "decoil", "bin", "pbsv2decoil.py")


def get_input_for_decoil_reconstruct(wildcards):
    """Resolve BAM + SV VCF + coverage bigwig + reference inputs.

    Pure path resolution — no ``os.path.exists`` checks. Snakemake raises
    ``MissingInputException`` automatically for missing sample-derived inputs
    (BAM, VCF, BW); reference files (ref, gtf) are pre-validated at module
    load time so users see a clear error if ``genome.fasta`` / ``genome.gtf``
    are absent from config.
    """
    sid = wildcards.sample_id
    in_dict = {}

    # BAM + BAI from indir
    if bam_substring:
        in_dict["bam"] = os.path.join(indir, sid, f"{sid}.{bam_substring}.bam")
        in_dict["bai"] = os.path.join(indir, sid, f"{sid}.{bam_substring}.bai")
    else:
        in_dict["bam"] = os.path.join(indir, sid, f"{sid}.bam")
        in_dict["bai"] = os.path.join(indir, sid, f"{sid}.bai")

    # External SV VCF (pbsv / sniffles upstream)
    in_dict["vcf"] = os.path.join(vcf_dir, sid, f"{sid}{vcf_suffix}")

    # External coverage bigwig (decoil/coverage submodule upstream)
    in_dict["bw"] = os.path.join(bw_dir, sid, f"{sid}{bw_suffix}")
    fasta = config.get("genome", {}).get("fasta")
    gtf = config.get("genome", {}).get("gtf")
    # Reference files (config-validated at module load)
    in_dict["ref"] = fasta
    in_dict["gtf"] = gtf
    if not fasta or not os.path.exists(fasta):
        raise ValueError("decoil module requires 'genome.fasta' in config. Please provide a valid reference genome FASTA path.")
    if not gtf or not os.path.exists(gtf):
        raise ValueError("decoil module requires 'genome.gtf' in config. Please provide a valid gene annotation GTF path.")
    return in_dict

rule decoil_reconstruct:
    """Run ``decoil reconstruct`` on one sample.
    Consumes BAM + external SV VCF + external coverage bigwig and produces
    reconstruct.bed / reconstruct.ecDNA.bed / reconstruct.ecDNA.filtered.bed
    / summary.txt. Does NOT perform SV calling internally.
    """
    input:
        unpack(get_input_for_decoil_reconstruct),
    output:
        reconstruct_all      = outdir + "/{sample_id}/decoil/reconstruct.bed",
        reconstruct_ecdna    = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.bed",
        reconstruct_filtered = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed",
        summary              = outdir + "/{sample_id}/decoil/summary.txt",
    log:
        logdir + "/{sample_id}/decoil_reconstruct.log"
    conda:
        "decoil.yaml"
    container:
        sif("decoil.yaml")
    params:
        # Tool executable
        decoil           = config.get("Procedure", {}).get("decoil") or "decoil",
        # Decoil `reconstruct` flags (each parameter flat in params, not dict-packed)
        sv_caller        = reconstruct_param("sv_caller") or "sniffles1",
        min_sv_len       = reconstruct_param("min_sv_len"),
        fragment_min_cov = reconstruct_param("fragment_min_cov"),
        fragment_max_cov = reconstruct_param("fragment_max_cov"),
        fragment_min_size= reconstruct_param("fragment_min_size"),
        min_vaf          = reconstruct_param("min_vaf"),
        min_cov_alt      = reconstruct_param("min_cov_alt"),
        min_cov          = reconstruct_param("min_cov"),
        max_explog       = reconstruct_param("max_explog_threshold"),
        filter_score     = reconstruct_param("filter_score"),
        extend_chr       = reconstruct_param("extend_allowed_chromosomes") or "",
        skip_fasta       = bool(reconstruct_param("skip")),
    threads: 16
    run:
        log_path = str(log)
        rule_logger = setup_logger(logger_name="decoil_reconstruct", log_file=log_path)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start decoil reconstruct for sample {wildcards.sample_id} "
                f"at {current_time})"
            )
            sample_outdir = os.path.dirname(str(output.reconstruct_all))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(
                sample_outdir,
                f"decoil_reconstruct_{wildcards.sample_id}_{current_time}.sh",
            )

            # `decoil reconstruct` has no --threads flag (Decoil 2.x argparse
            # rejects it), reads VCFs as plain text (no bgzip), and requires
            # sniffles1-style fields (STRANDS/CHR2/END + FORMAT DR/DV) that
            # pbsv does not emit. Convert the input VCF first — verified
            # against decoil.sif 2.0.2 with a full reconstruction run.
            converted_vcf = os.path.join(sample_outdir, "decoil_input.vcf")

            # `decoil reconstruct` is reconstruction-only — consumes pre-computed
            # SV VCF and coverage bigwig. No internal sniffles call.
            cmd = [
                params.decoil, "reconstruct",
                "-b", input.bam,
                "-i", converted_vcf,
                "-c", input.bw,
                "-r", input.ref,
                "-g", input.gtf,
                "-o", sample_outdir,
                "--name", wildcards.sample_id,
                "--sv-caller", params.sv_caller,
                "--min-sv-len",        str(params.min_sv_len),
                "--fragment-min-cov",  str(params.fragment_min_cov),
                "--fragment-max-cov",  str(params.fragment_max_cov),
                "--fragment-min-size", str(params.fragment_min_size),
                "--min-vaf",           str(params.min_vaf),
                "--min-cov-alt",       str(params.min_cov_alt),
                "--min-cov",           str(params.min_cov),
                "--max-explog-threshold", str(params.max_explog),
                "--filter-score",      str(params.filter_score),
            ]
            if params.extend_chr:
                cmd += ["--extend-allowed-chr", params.extend_chr]
            if params.skip_fasta:
                cmd.append("--skip")

            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(
                    f"python3 {PBSV2DECOIL_SCRIPT} "
                    f"{input.vcf} {converted_vcf}\n"
                )
                f.write(" ".join(cmd) + "\n")
                f.write(
                    f'echo "decoil_reconstruct for {wildcards.sample_id} at '
                    f'{current_time} completed successfully"\n'
                )
            rule_logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"decoil_reconstruct failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"decoil_reconstruct failed for sample {wildcards.sample_id}: {e}")


# ============================================================
# Result aggregation rule
# ============================================================
rule decoil_result:
    """Result aggregation — touch sentinel for Snakemake DAG tracking."""
    input:
        reconstruct_all      = outdir + "/{sample_id}/decoil/reconstruct.bed",
        reconstruct_ecdna    = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.bed",
        reconstruct_filtered = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed",
        summary              = outdir + "/{sample_id}/decoil/summary.txt",


# ============================================================
# Visualization: per-sample statistics plots + cross-sample compare
# Scripts located via ROOT_DIR (project convention, see pitfall #55).
# ============================================================
DECOIL_PLOT_SCRIPT = os.path.join(ROOT_DIR, "modules", "decoil", "bin", "decoil_plot.py")
DECOIL_IGV_CONFIG_SCRIPT = os.path.join(ROOT_DIR, "modules", "decoil", "bin", "decoil_igv_config.py")
TRACK_SCRIPT = os.path.join(ROOT_DIR, "modules", "track", "bin", "track.py")

# IGV base config (reference fasta/index/cytoband/GTF + publicPathMap) —
# same structure as modules/track/track.json 'igv' section.
igv_cfg = config.get("Params", {}).get("igv", {})


rule decoil_plot:
    """Per-sample ecDNA statistics: size distribution, topology, region table.

    Consumes the sample's own summary.txt + filtered ecDNA BED and produces
    plots + TSV tables under {outdir}/{sample_id}/visualization/.
    """
    input:
        summary              = outdir + "/{sample_id}/decoil/summary.txt",
        reconstruct_filtered = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed",
    output:
        png                  = outdir + "/{sample_id}/visualization/decoil_stats.png",
        regions_tsv          = outdir + "/{sample_id}/visualization/decoil_stats_ecdna_regions.tsv",
    log:
        logdir + "/{sample_id}/decoil_plot.log"
    conda:
        "decoil.yaml"
    container:
        sif("decoil.yaml")
    params:
        plot_script = DECOIL_PLOT_SCRIPT,
        formats     = PLOT_FORMATS,
    run:
        log_path = str(log)
        rule_logger = setup_logger(logger_name="decoil_plot", log_file=log_path)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start decoil_plot for sample {wildcards.sample_id} at {current_time}")
            sample_vis_dir = os.path.dirname(str(output.png))
            os.makedirs(sample_vis_dir, exist_ok=True)
            script = os.path.join(sample_vis_dir, f"decoil_plot_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.plot_script,
                "--mode", "single",
                "--sample", f"{wildcards.sample_id}:{input.summary}",
                "--ecdna-bed", f"{wildcards.sample_id}:{input.reconstruct_filtered}",
                "-o", os.path.join(sample_vis_dir, "decoil_stats"),
            ]
            for fmt in params.formats:
                cmd += ["-f", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            rule_logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"decoil_plot failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"decoil_plot failed for sample {wildcards.sample_id}: {e}")


rule decoil_compare:
    """Cross-sample ecDNA comparison panel + shared/unique region table.

    Aggregates every sample's summary.txt + filtered ecDNA BED into one
    comparison figure and a shared/unique TSV under {outdir}/visualization/.
    """
    input:
        summaries = expand(outdir + "/{sample_id}/decoil/summary.txt", sample_id=samples),
        beds      = expand(outdir + "/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed", sample_id=samples),
    output:
        png            = outdir + "/visualization/decoil_compare.png",
        shared_tsv     = outdir + "/visualization/decoil_compare_shared_unique.tsv",
        regions_tsv    = outdir + "/visualization/decoil_compare_ecdna_regions.tsv",
    log:
        logdir + "/decoil_compare.log"
    conda:
        "decoil.yaml"
    container:
        sif("decoil.yaml")
    params:
        plot_script = DECOIL_PLOT_SCRIPT,
        formats     = COMPARE_FORMATS,
    run:
        log_path = str(log)
        rule_logger = setup_logger(logger_name="decoil_compare", log_file=log_path)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start decoil_compare for samples {samples} at {current_time}")
            compare_dir = os.path.dirname(str(output.png))
            os.makedirs(compare_dir, exist_ok=True)
            script = os.path.join(compare_dir, f"decoil_compare_{current_time}.sh")
            cmd = [
                "python", params.plot_script,
                "--mode", "compare",
            ]
            for sid, summary_path in zip(samples, list(input.summaries)):
                cmd += ["--sample", f"{sid}:{summary_path}"]
            for sid, bed_path in zip(samples, list(input.beds)):
                cmd += ["--ecdna-bed", f"{sid}:{bed_path}"]
            cmd += ["-o", os.path.join(compare_dir, "decoil_compare")]
            for fmt in params.formats:
                cmd += ["-f", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            rule_logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"decoil_compare failed: {e}")
            raise RuntimeError(f"decoil_compare failed: {e}")


rule decoil_igv:
    """Interactive IGV.js HTML per sample: coverage bigwig + ecDNA BED + GTF.

    Reuses modules/track/bin/track.py (``--mode igv``) with a Decoil-specific
    browser config generated by bin/decoil_igv_config.py. The HTML opens on
    the sample's largest ecDNA region by default.
    """
    input:
        reconstruct_filtered = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed",
        bw                   = os.path.join(bw_dir, "{sample_id}", "{sample_id}" + bw_suffix),
    output:
        html                 = outdir + "/{sample_id}/visualization/igv_ecdna.html",
    log:
        logdir + "/{sample_id}/decoil_igv.log"
    conda:
        "decoil.yaml"
    container:
        sif("decoil.yaml")
    params:
        igv_config_script = DECOIL_IGV_CONFIG_SCRIPT,
        track_script      = TRACK_SCRIPT,
        base_igv_config   = igv_cfg,
        config_json       = outdir + "/{sample_id}/visualization/igv_ecdna_config.json",
        locus             = IGV_LOCUS,
    run:
        log_path = str(log)
        rule_logger = setup_logger(logger_name="decoil_igv", log_file=log_path)
        open(log_path, "w").close()
        try:
            import json
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start decoil_igv for sample {wildcards.sample_id} at {current_time}")
            if not params.base_igv_config:
                raise ValueError(
                    "decoil_igv requires 'igv' config (reference fastaURL/"
                    "indexURL/GTF tracks + publicPathMap). Add an 'igv' "
                    "section to the workflow config, see modules/track/track.json."
                )
            sample_vis_dir = os.path.dirname(str(output.html))
            os.makedirs(sample_vis_dir, exist_ok=True)

            # Render the config template with this sample's concrete paths,
            # then let decoil_igv_config.py add bigwig/bed tracks + locus.
            base_cfg = dict(params.base_igv_config)
            config_json_path = str(params.config_json)
            tmp_base_json = os.path.join(sample_vis_dir, f"igv_base_{wildcards.sample_id}.json")
            with open(tmp_base_json, "w") as f:
                json.dump(base_cfg, f, indent=2)

            script = os.path.join(sample_vis_dir, f"decoil_igv_{wildcards.sample_id}_{current_time}.sh")
            cmd1 = [
                "python", params.igv_config_script,
                "--sample", f"{wildcards.sample_id}:{input.reconstruct_filtered}",
                "--bw", f"{wildcards.sample_id}:{input.bw}",
                "--igv-config", tmp_base_json,
                "-o", config_json_path,
            ]
            if params.locus:
                cmd1 += ["--locus", str(params.locus)]
            cmd2 = [
                "python", params.track_script,
                "--mode", "igv",
                "--config", config_json_path,
                "-i", input.bw,
                "-o", str(output.html),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd1) + "\n")
                f.write(" ".join(cmd2) + "\n")
            rule_logger.info("Executing: " + " ".join(cmd1) + " && " + " ".join(cmd2))
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"decoil_igv failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"decoil_igv failed for sample {wildcards.sample_id}: {e}")


rule decoil_ucsc:
    """UCSC Genome Browser track line file per sample.

    Same input set as decoil_igv (coverage bigwig + ecDNA BED) but emits a
    UCSC-style track hub text (``track type=... bigDataUrl=...`` lines) via
    track.py ``--mode ucsc``, for loading the sample tracks into a UCSC
    Genome Browser session or any hub-compatible viewer.
    """
    input:
        reconstruct_filtered = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed",
        bw                   = os.path.join(bw_dir, "{sample_id}", "{sample_id}" + bw_suffix),
    output:
        track                = outdir + "/{sample_id}/visualization/ucsc_track.txt",
    log:
        logdir + "/{sample_id}/decoil_ucsc.log"
    conda:
        "decoil.yaml"
    container:
        sif("decoil.yaml")
    params:
        track_script = TRACK_SCRIPT,
    run:
        log_path = str(log)
        rule_logger = setup_logger(logger_name="decoil_ucsc", log_file=log_path)
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start decoil_ucsc for sample {wildcards.sample_id} at {current_time}")
            sample_vis_dir = os.path.dirname(str(output.track))
            os.makedirs(sample_vis_dir, exist_ok=True)
            script = os.path.join(sample_vis_dir, f"decoil_ucsc_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                "python", params.track_script,
                "--mode", "ucsc",
                "-i", input.bw, input.reconstruct_filtered,
                "-o", str(output.track),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(" ".join(cmd) + "\n")
            rule_logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"decoil_ucsc failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"decoil_ucsc failed for sample {wildcards.sample_id}: {e}")