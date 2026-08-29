"""Fiber-seq module for single-molecule chromatin accessibility analysis.

Provides rules for:
  - ft_predict_m6a: Predict m6A modifications from PacBio HiFi kinetics
  - ft_add_nucleosomes: Add nucleosome calls to Fiber-seq BAM
  - ft_fire: Call Fiber-seq Inferred Regulatory Elements (FIREs)
  - ft_extract: Extract Fiber-seq data to BED/bigBed format

Requires fibertools-rs (ft) CLI tool.
See: https://fiberseq.github.io/
"""
include: "../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta") or ""


rule ft_predict_m6a:
    """Predict m6A positions from PacBio HiFi CCS BAM with kinetics.

    Input: PacBio CCS BAM with polymerase kinetics (IPD/PL) tags.
    Output: Fiber-seq BAM with m6A calls in MM/ML tags.
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
            rule_logger.info(f"Predicting m6A for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
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

    Input: Fiber-seq BAM with m6A calls (from ft_predict_m6a or dorado).
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


rule ft_extract:
    """Extract Fiber-seq data to BED format.

    Input: Fiber-seq BAM with m6A, nucleosome, and optionally FIRE calls.
    Output: Compressed BED12 files for m6A, nucleosomes, MSPs, and FIREs.
    """
    input:
        bam = outdir + "/{sample_id}/{sample_id}.fiberseq.fire.bam"
    output:
        m6a = outdir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        nuc = outdir + "/{sample_id}/{sample_id}.nuc.bed.gz",
        msp = outdir + "/{sample_id}/{sample_id}.msp.bed.gz",
        fire = outdir + "/{sample_id}/{sample_id}.fire.bed.gz",
    log:
        logdir + "/{sample_id}/ft_extract.log"
    threads: 8
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
            rule_logger = setup_logger("ft_extract", log_file=log_path)
            rule_logger.info(f"Extracting Fiber-seq data for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(sample_outdir, f"ft_extract_{current_time}.sh")
            cmd = [
                params.ft, "extract",
                "-t", str(threads),
                "--m6a", output.m6a,
                "--nuc", output.nuc,
                "--msp", output.msp,
                "--fire", output.fire,
                input.bam,
            ]
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Fiber-seq extraction completed for sample {wildcards.sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"ft_extract failed: {e}\n")
            raise RuntimeError(f"ft_extract failed: {e}\n")
