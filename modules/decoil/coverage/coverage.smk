"""Coverage bigwig generation for Decoil ecDNA reconstruction.

Generates a per-sample genome-wide coverage bigwig (``bamCoverage`` from
deeptools) directly from the markdup BAM. The bigwig is consumed by the
``decoil`` module as the ``-c / --coverage`` argument to ``decoil reconstruct``.

Decoil needs a single bigwig per sample; we do not perform samtools dedup
before bamCoverage because PacBio HiFi markdup BAMs are already duplicate-free
in practice and a second dedup pass would drop signal.
"""

include: "../../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

bam_substring = config.get("bam_substring") or ""
bw_suffix = config.get("bw_suffix") or ".coverage.bw"


def get_bam_input(wildcards):
    """Resolve markdup BAM + BAI inputs."""
    sid = wildcards.sample_id
    if bam_substring:
        return {
            "bam": os.path.join(indir, sid, f"{sid}.{bam_substring}.bam"),
            "bai": os.path.join(indir, sid, f"{sid}.{bam_substring}.bai"),
        }
    return {
        "bam": os.path.join(indir, sid, f"{sid}.bam"),
        "bai": os.path.join(indir, sid, f"{sid}.bai"),
    }


rule decoil_coverage:
    """Generate per-sample genome-wide coverage bigwig for Decoil."""
    input:
        unpack(get_bam_input),
    output:
        bw = outdir + "/{sample_id}/{sample_id}" + bw_suffix,
    log:
        logdir + "/{sample_id}/decoil_coverage.log"
    threads: 8
    conda:
        "coverage.yaml"
    container:
        sif("coverage.yaml")
    params:
        bamCoverage = config.get("Procedure", {}).get("bamCoverage") or "bamCoverage",
        bin_size    = config.get("Params", {}).get("coverage", {}).get("bin_size", 1000),
        normalize   = config.get("Params", {}).get("coverage", {}).get("normalize", "RPGC"),
        outdir_sample = outdir + "/{sample_id}",
    run:
        log_path = str(log)
        rule_logger = setup_logger(
            logger_name="decoil_coverage", log_file=log_path
        )
        open(log_path, "w").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start decoil_coverage for sample {wildcards.sample_id} "
                f"at {current_time}"
            )
            os.makedirs(str(params.outdir_sample), exist_ok=True)
            script = os.path.join(
                str(params.outdir_sample),
                f"decoil_coverage_{wildcards.sample_id}_{current_time}.sh",
            )
            cmd = [
                params.bamCoverage,
                "--numberOfProcessors", str(threads),
                "--binSize", str(params.bin_size),
                "--normalizeUsing", params.normalize,
                "-b", input.bam,
                "-o", output.bw,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(
                    f'echo "decoil_coverage for {wildcards.sample_id} at '
                    f'{current_time} completed successfully"\n'
                )
            rule_logger.info("Executing: " + " ".join(cmd))
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(
                    f"decoil_coverage failed for sample {wildcards.sample_id}: {e}\n"
                )
            rule_logger.error(
                f"decoil_coverage failed for sample {wildcards.sample_id}: {e}"
            )
            raise RuntimeError(
                f"decoil_coverage failed for sample {wildcards.sample_id}: {e}"
            )


rule decoil_coverage_result:
    """Result aggregation rule."""
    input:
        bw = outdir + "/{sample_id}/{sample_id}" + bw_suffix,