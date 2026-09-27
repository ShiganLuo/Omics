"""Fiber-seq prepare module: m6A prediction + nucleosome + FIRE calling.

Auto-detects whether input BAM already has m6A calls (e.g. SPRQ chemistry):
  - No m6A → ft predict-m6a (adds m6A + nuc) → ft fire
  - Has m6A → skip predict-m6a → ft add-nucleosomes → ft fire

Requires fibertools-rs (ft) CLI tool.
See: https://fiberseq.github.io/
"""
include: "../common/common.smk"
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])


def _detect_m6a(bam_path: str) -> bool:
    """Check if BAM already has m6A MM tags in first 100 reads."""
    import subprocess
    try:
        r = subprocess.run(
            ["samtools", "view", bam_path],
            capture_output=True, text=True, timeout=30
        )
        lines = r.stdout.split("\n")[:100]
        m6a_count = sum(1 for line in lines if "MM:Z:" in line)
        return m6a_count > 50
    except Exception:
        return False


rule ft_predict_m6a:
    """Predict m6A positions from PacBio HiFi CCS BAM with kinetics.

    Auto-detects: if BAM already has m6A calls (SPRQ chemistry), copies input as output.
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
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("ft_predict_m6a", log_file=log_path)
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)

            # Auto-detect m6A
            has_m6a = _detect_m6a(input.bam)
            if has_m6a:
                rule_logger.info(f"BAM already has m6A calls (SPRQ?), skipping predict-m6a")
                # Copy BAM to output
                import shutil
                shutil.copy2(input.bam, output.bam)
                rule_logger.info(f"Copied {input.bam} → {output.bam}")
            else:
                rule_logger.info(f"Predicting m6A for sample {wildcards.sample_id}")
                current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
                command_script = os.path.join(sample_outdir, f"ft_predict_m6a_{current_time}.sh")
                cmd = [
                    params.ft, "predict-m6a",
                    "-t", str(threads),
                    input.bam,
                    output.bam,
                ]
                with open(command_script, "w") as f:
                    f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                    f.write(" ".join(cmd) + "\n")
                    f.write(f'echo "m6A prediction completed for sample {wildcards.sample_id}"\n')
                shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"ft_predict_m6a failed: {e}\n")
            raise RuntimeError(f"ft_predict_m6a failed: {e}\n")


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
        bam = outdir + "/{sample_id}/{sample_id}.fiberseq.nuc.bam"
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