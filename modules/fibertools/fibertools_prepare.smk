"""Fiber-seq prepare module: m6A prediction + nucleosome + FIRE calling.

Auto-detects whether input BAM already has m6A calls (e.g. SPRQ chemistry):
  - No m6A → ft predict-m6a (adds m6A + nuc) → ft fire
  - Has m6A → skip predict-m6a → ft add-nucleosomes → ft fire

Detection logic runs in shell script (inside container) for reliable path/samtools access.

Requires fibertools-rs (ft) CLI tool.
See: https://fiberseq.github.io/
"""
include: "../common/common.smk"
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])


def _detect_m6a_sh(bam_path: str) -> str:
    """Return shell snippet that sets HAS_M6A=1 or HAS_M6A=0.

    Checks @RG header for SPRQ binding kit (103-496-900),
    falls back to checking first 100 reads for MM:Z: tags.
    """
    return f"""\
HAS_M6A=0
if samtools view -H {bam_path} 2>/dev/null | grep -q '103-496-900'; then
    HAS_M6A=1
elif samtools view {bam_path} 2>/dev/null | head -100 | grep -q 'MM:Z:'; then
    HAS_M6A=1
fi
"""


# Pre-detect m6A at DAG time (runs on host) for fire input routing
_fibertools_sif = sif("fibertools.yaml")

def _has_m6a(sample_id):
    """Check header for SPRQ binding kit via container's samtools."""
    bam = os.path.join(indir, sample_id, f"{sample_id}.bam")
    if not os.path.exists(bam):
        return False
    try:
        real_bam = os.path.realpath(bam)
        bam_dir = os.path.dirname(real_bam)
        cmd = ["apptainer", "exec", "--bind", f"{bam_dir}:{bam_dir}", _fibertools_sif, "samtools", "view", "-H", real_bam]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        return "103-496-900" in r.stdout
    except Exception:
        return False


rule ft_predict_m6a:
    """Predict m6A positions from PacBio HiFi CCS BAM with kinetics.

    Auto-detects: if BAM already has m6A calls (SPRQ chemistry), symlinks input as output.
    Detection runs in shell script inside container for reliable access.
    """
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        bam = outdir + "/{sample_id}/{sample_id}.fiberseq.bam"
    log:
        logdir + "/{sample_id}/ft_predict_m6a.log"
    threads: 16
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft"
    run:
        log_path = str(log)
        sample_outdir = os.path.join(outdir, wildcards.sample_id)
        os.makedirs(sample_outdir, exist_ok=True)
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        command_script = os.path.join(sample_outdir, f"ft_predict_m6a_{current_time}.sh")

        detect_sh = _detect_m6a_sh(input.bam)
        cmd_predict = " ".join([params.ft, "predict-m6a", "-t", str(threads), input.bam, output.bam])

        with open(command_script, "w") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(f'echo "Detecting m6A in {wildcards.sample_id}..."\n')
            f.write(detect_sh)
            f.write('if [ "$HAS_M6A" -eq 1 ]; then\n')
            f.write(f'    echo "BAM already has m6A calls (SPRQ?), skipping predict-m6a"\n')
            f.write(f'    ln -sf "$(realpath {input.bam})" {output.bam}\n')
            f.write(f'    echo "Symlinked {input.bam} -> {output.bam}"\n')
            f.write('else\n')
            f.write(f'    echo "Predicting m6A for sample {wildcards.sample_id}"\n')
            f.write(f'    {cmd_predict}\n')
            f.write(f'    echo "m6A prediction completed for sample {wildcards.sample_id}"\n')
            f.write('fi\n')

        shell(f"bash {command_script} >> {log_path} 2>&1")


def get_input_for_ft_fire(wildcards):
    """Determine fire input based on m6A detection and manual config.

    Priority: auto-detect > manual config
    - SPRQ (has m6A): fire reads nuc.bam
    - Non-SPRQ (no m6A): fire reads fiberseq.bam (predict-m6a already has nuc)
    - add_nucleosomes_manual=true: force fire reads nuc.bam
    """
    add_nucleosomes_manual = config.get("Params", {}).get("fibertools", {}).get("add_nucleosomes_manual", None)
    if _has_m6a(wildcards.sample_id) or add_nucleosomes_manual:
        return outdir + "/{sample_id}/{sample_id}.fiberseq.nuc.bam"
    else:
        return outdir + "/{sample_id}/{sample_id}.fiberseq.bam"


rule ft_add_nucleosomes:
    """Add nucleosome calls to Fiber-seq BAM with m6A predictions.

    Input: Fiber-seq BAM with m6A calls (from ft_predict_m6a or SPRQ).
    Output: Fiber-seq BAM with nucleosome and MSP calls added.
    """
    input:
        bam = outdir + "/{sample_id}/{sample_id}.fiberseq.bam"
    output:
        bam = outdir + "/{sample_id}/{sample_id}.fiberseq.nuc.bam"
    log:
        logdir + "/{sample_id}/ft_add_nucleosomes.log"
    threads: 16
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft"
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("ft_add_nucleosomes", log_file=log_path)
            rule_logger.info(f"Adding nucleosomes for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(sample_outdir, f"ft_add_nucleosomes_{current_time}.sh")
            cmd = [
                params.ft, "add-nucleosomes",
                "-t", str(threads),
                input.bam,
                output.bam,
            ]
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Nucleosome calling completed for sample {wildcards.sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"ft_add_nucleosomes failed: {e}\n")
            raise RuntimeError(f"ft_add_nucleosomes failed: {e}\n")


rule ft_fire:
    """Call Fiber-seq Inferred Regulatory Elements (FIREs).

    Input: Fiber-seq BAM with m6A and nucleosome calls.
    Output: Fiber-seq BAM with FIRE calls in aq tags.
    """
    input:
        bam = get_input_for_ft_fire
    output:
        bam = outdir + "/{sample_id}/{sample_id}.fiberseq.fire.bam"
    log:
        logdir + "/{sample_id}/ft_fire.log"
    threads: 8
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
        is_ont = config.get("Params", {}).get("fibertools", {}).get("ont", False)
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("ft_fire", log_file=log_path)
            rule_logger.info(f"Calling FIREs for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(sample_outdir, f"ft_fire_{current_time}.sh")
            cmd = [
                params.ft, "fire",
                "-t", str(threads),
                input.bam,
                output.bam,
            ]
            if params.is_ont:
                cmd.extend(["--ont"])
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "FIRE calling completed for sample {wildcards.sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"ft_fire failed: {e}\n")
            raise RuntimeError(f"ft_fire failed: {e}\n")