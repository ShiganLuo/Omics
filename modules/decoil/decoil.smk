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

# Reference files — fail fast on missing config paths (snakemake's own
# MissingInputException would otherwise report a "None" path which is
# confusing for users). These are external resources, not derived inputs.
fasta = config.get("genome", {}).get("fasta")
gtf = config.get("genome", {}).get("gtf")
if not fasta:
    raise ValueError(
        "decoil module requires 'genome.fasta' in config. "
        "Please provide a valid reference genome FASTA path."
    )
if not gtf:
    raise ValueError(
        "decoil module requires 'genome.gtf' (gene annotation) in config. "
        "Decoil reconstruct uses annotated gene boundaries."
    )

# Upstream file locations (sample-derived paths — no validation here;
# snakemake's MissingInputException is the correct error if the upstream
# module did not produce them)
vcf_dir = config.get("vcf_dir")
bw_dir = config.get("bw_dir") or outdir
vcf_suffix = config.get("vcf_suffix") or ".sv.vcf.gz"
bw_suffix = config.get("bw_suffix") or ".bigwig"


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

    # Reference files (config-validated at module load)
    in_dict["ref"] = fasta
    in_dict["gtf"] = gtf

    return in_dict


# ============================================================
# Rule: Run Decoil reconstruct (ecDNA reconstruction only)
# ============================================================
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
        sv_caller        = config.get("Params", {}).get("decoil", {}).get("sv_caller") or "sniffles1",
        min_sv_len       = config.get("Params", {}).get("decoil", {}).get("min_sv_len", 500),
        fragment_min_cov = config.get("Params", {}).get("decoil", {}).get("fragment_min_cov", 5),
        fragment_max_cov = config.get("Params", {}).get("decoil", {}).get("fragment_max_cov", 100000),
        fragment_min_size= config.get("Params", {}).get("decoil", {}).get("fragment_min_size", 500),
        min_vaf          = config.get("Params", {}).get("decoil", {}).get("min_vaf", 0.01),
        min_cov_alt      = config.get("Params", {}).get("decoil", {}).get("min_cov_alt", 6),
        min_cov          = config.get("Params", {}).get("decoil", {}).get("min_cov", 8),
        max_explog       = config.get("Params", {}).get("decoil", {}).get("max_explog_threshold", 0.1),
        filter_score     = config.get("Params", {}).get("decoil", {}).get("filter_score", 5),
        extend_chr       = config.get("Params", {}).get("decoil", {}).get("extend_allowed_chromosomes") or "",
        skip_fasta       = bool(config.get("Params", {}).get("decoil", {}).get("skip", False)),
        threads          = config.get("Params", {}).get("decoil", {}).get("threads", 16),
    threads: 16
    resources:
        mem_mb = 30000
    run:
        log_path = str(log)
        # `logger` is bound at module scope by common.smk
        # (`from snakemake.logging import logger`). DO NOT rebind it here:
        # snakemake's own logger handlers must remain intact. Use a distinct
        # local name for our per-sample file-handler logger so the run body
        # and the except branch share the same reference without clashing
        # with snakemake's logger.
        rule_logger = setup_logger(
            logger_name="decoil_reconstruct", log_file=log_path
        )
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start decoil reconstruct for sample {wildcards.sample_id} "
                f"at {current_time} (threads={params.threads})"
            )
            sample_outdir = os.path.dirname(str(output.reconstruct_all))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(
                sample_outdir,
                f"decoil_reconstruct_{wildcards.sample_id}_{current_time}.sh",
            )

            # `decoil reconstruct` is reconstruction-only — consumes pre-computed
            # SV VCF and coverage bigwig. No internal sniffles call.
            cmd = [
                params.decoil, "reconstruct",
                "--threads", str(params.threads),
                "-b", input.bam,
                "-i", input.vcf,
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
                f.write(" ".join(cmd) + "\n")
                f.write(
                    f'echo "decoil_reconstruct for {wildcards.sample_id} at '
                    f'{current_time} completed successfully"\n'
                )
            rule_logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(
                    f"decoil_reconstruct failed for sample {wildcards.sample_id}: {e}\n"
                )
            rule_logger.error(
                f"decoil_reconstruct failed for sample {wildcards.sample_id}: {e}"
            )
            raise RuntimeError(
                f"decoil_reconstruct failed for sample {wildcards.sample_id}: {e}"
            )


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