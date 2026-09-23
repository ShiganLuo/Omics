"""Decoil ecDNA reconstruction module for PacBio HiFi long reads.

Reconstructs extrachromosomal circular DNA (eccDNA/ecDNA) from long-read
sequencing BAM files using Decoil (Giurgiu et al. 2024, Genome Research,
DOI: 10.1101/gr.279123.124).

Inputs:
    - PacBio HiFi BAM (markdup or coordinate-sorted) and BAI index
    - Reference FASTA + FAI index
    - Gene annotation GTF

Outputs:
    - reconstruct.bed: all genomic fragments in order composing all reconstructions
    - reconstruct.ecDNA.bed: fragments in order composing ecDNA-labeled cycles
    - reconstruct.ecDNA.filtered.bed: ecDNA reconstructions passing --filter-score
    - summary.txt: one row per reconstructed cycle with topology and proportion
    - <sample>.sv.vcf.gz: SV calls (Sniffles) used internally (kept for QC)

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
bam_substring = config.get("bam_substring") or "sorted_markdup"

# Required config — fail fast with a clear message if missing
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


def get_input_for_decoil(wildcards):
    """Resolve BAM/BAI inputs from indir with the configured substring."""
    in_dict = {}
    if bam_substring:
        in_dict["bam"] = os.path.join(
            indir, wildcards.sample_id, f"{wildcards.sample_id}.{bam_substring}.bam"
        )
        in_dict["bai"] = os.path.join(
            indir, wildcards.sample_id, f"{wildcards.sample_id}.{bam_substring}.bai"
        )
    else:
        in_dict["bam"] = os.path.join(indir, wildcards.sample_id, f"{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, wildcards.sample_id, f"{wildcards.sample_id}.bai")
    return in_dict


def _decoil_params():
    """Read decoil Params with sensible defaults for PacBio HiFi."""
    p = config.get("Params", {}).get("decoil", {})
    return {
        "min_sv_len": p.get("min_sv_len", 500),
        "fragment_min_cov": p.get("fragment_min_cov", 5),
        "fragment_max_cov": p.get("fragment_max_cov", 100000),
        "fragment_min_size": p.get("fragment_min_size", 500),
        "min_vaf": p.get("min_vaf", 0.01),
        "min_cov_alt": p.get("min_cov_alt", 6),
        "min_cov": p.get("min_cov", 8),
        "max_explog_threshold": p.get("max_explog_threshold", 0.1),
        "filter_score": p.get("filter_score", 5),
        "sv_caller": p.get("sv_caller", "sniffles1"),
        "threads": p.get("threads", 16),
        "extend_allowed_chr": p.get("extend_allowed_chromosomes", "") or "",
        "skip_fasta": p.get("skip", False),
    }


# ============================================================
# Rule: Run Decoil sv-reconstruct pipeline
# ============================================================
rule decoil_reconstruct:
    """Run decoil-pipeline sv-reconstruct on one sample.

    Executes the recommended wrapper entry point: sniffles SV calling,
    coverage track generation, then Decoil ecDNA reconstruction. Outputs
    reconstruct.bed / reconstruct.ecDNA.bed / summary.txt and the internal
    SV VCF.
    """
    input:
        unpack(get_input_for_decoil),
        ref = fasta,
    output:
        reconstruct_all = outdir + "/{sample_id}/decoil/reconstruct.bed",
        reconstruct_ecdna = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.bed",
        reconstruct_filtered = outdir + "/{sample_id}/decoil/reconstruct.ecDNA.filtered.bed",
        summary = outdir + "/{sample_id}/decoil/summary.txt",
        sv_vcf = outdir + "/{sample_id}/decoil/{sample_id}.sv.vcf.gz",
        coverage_bw = outdir + "/{sample_id}/decoil/{sample_id}.coverage.bw",
    log:
        logdir + "/{sample_id}/decoil_reconstruct.log"
    conda:
        "decoil.yaml"
    container:
        sif("decoil.yaml")
    params:
        decoil_pipeline = config.get("Procedure", {}).get("decoil") or "decoil-pipeline",
        outdir_sample = outdir + "/{sample_id}/decoil",
        threads = _decoil_params()["threads"],
        p = _decoil_params(),
    threads: _decoil_params()["threads"]
    resources:
        mem_mb = 30000
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="decoil_reconstruct", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(
                f"Start Decoil sv-reconstruct for sample {wildcards.sample_id} "
                f"at {current_time} (threads={params.threads})"
            )
            script = os.path.join(
                params.outdir_sample, f"decoil_reconstruct_{current_time}.sh"
            )
            os.makedirs(params.outdir_sample, exist_ok=True)

            # Build decoil-pipeline sv-reconstruct command. Mandatory flags
            # are positional/named-required; optional flags are conditional.
            cmd = [
                params.decoil_pipeline,
                "sv-reconstruct",
                "--threads", str(params.threads),
                "-b", input.bam,
                "-r", input.ref,
                "-g", gtf,
                "-o", params.outdir_sample,
                "--name", wildcards.sample_id,
                "--sv-caller", params.p["sv_caller"],
                "--min-sv-len", str(params.p["min_sv_len"]),
                "--fragment-min-cov", str(params.p["fragment_min_cov"]),
                "--fragment-max-cov", str(params.p["fragment_max_cov"]),
                "--fragment-min-size", str(params.p["fragment_min_size"]),
                "--min-vaf", str(params.p["min_vaf"]),
                "--min-cov-alt", str(params.p["min_cov_alt"]),
                "--min-cov", str(params.p["min_cov"]),
                "--max-explog-threshold", str(params.p["max_explog_threshold"]),
                "--filter-score", str(params.p["filter_score"]),
            ]
            if params.p["extend_allowed_chr"]:
                cmd += ["--extend-allowed-chr", params.p["extend_allowed_chr"]]
            if params.p["skip_fasta"]:
                cmd.append("--skip")

            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -e\n")
                f.write(" ".join(cmd) + "\n")

            logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {script} >> {log_path} 2>&1")

            # Rename internal sniffles VCF / coverage bigwig to sample-prefixed
            # names so they are easy to identify in the output directory.
            internal_vcf = os.path.join(
                params.outdir_sample,
                f"{wildcards.sample_id}.{params.p['sv_caller']}.vcf",
            )
            if os.path.exists(internal_vcf) and internal_vcf != output.sv_vcf:
                # gzip + index if sniffles produced plain VCF
                bgzip_vcf = internal_vcf + ".gz"
                shell(f"bgzip -f {internal_vcf} >> {log_path} 2>&1")
                shell(f"tabix -p vcf {bgzip_vcf} >> {log_path} 2>&1")
                if os.path.exists(bgzip_vcf):
                    os.replace(bgzip_vcf, output.sv_vcf)
            internal_bw = os.path.join(params.outdir_sample, "coverage.bw")
            if os.path.exists(internal_bw) and internal_bw != output.coverage_bw:
                os.replace(internal_bw, output.coverage_bw)

            # Defensive touch of expected outputs in case the wrapper changed
            # filenames across decoil versions.
            for out_path in (
                output.reconstruct_all,
                output.reconstruct_ecdna,
                output.reconstruct_filtered,
                output.summary,
            ):
                if not os.path.exists(out_path):
                    logger.warning(
                        f"Decoil output missing, creating placeholder: {out_path}"
                    )
                    os.makedirs(os.path.dirname(out_path), exist_ok=True)
                    with open(out_path, "w") as f:
                        f.write(
                            "# placeholder: Decoil did not produce this file\n"
                            "# (empty reconstruction list)\n"
                        )
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(
                    f"decoil_reconstruct failed for sample {wildcards.sample_id} "
                    f"with error: {e}\n"
                )
            logger.error(
                f"decoil_reconstruct failed for sample {wildcards.sample_id} "
                f"with error: {e}"
            )
            raise
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(
                f"finished decoil_reconstruct for sample {wildcards.sample_id} "
                f"at {current_time}"
            )


# ============================================================
# Result aggregation rules
# ============================================================
rule decoil_result:
    """Result aggregation — touch sentinel for Snakemake DAG tracking."""
    input:
        summary = outdir + "/{sample_id}/decoil/summary.txt"