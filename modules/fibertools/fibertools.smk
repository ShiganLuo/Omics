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
upstream_outdir = config.get("upstream_outdir", indir)
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
        bam = upstream_outdir + "/{sample_id}/{sample_id}.bam"
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
        bam = upstream_outdir + "/{sample_id}/{sample_id}.fiberseq.nuc.bam"
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

    Input: Fiber-seq BAM with m6A, nucleosome, and FIRE calls.
    Output: Compressed BED file with all Fiber-seq annotations.
    Note: ft extract v0.13.1 --m6a/--nuc/--msp individual outputs are buggy;
          use --all to get all data in one file.
    """
    input:
        bam = upstream_outdir + "/{sample_id}/{sample_id}.fiberseq.fire.bam"
    output:
        all_bed = outdir + "/{sample_id}/{sample_id}.fiberseq.all.bed.gz",
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
                "--all", output.all_bed,
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


rule ft_call_peaks:
    """Call FIRE peaks using FDR-based peak calling.

    Input: Fiber-seq BAM with FIRE calls.
    Output: BED file with called FIRE peaks.
    """
    input:
        bam = upstream_outdir + "/{sample_id}/{sample_id}.fiberseq.fire.bam"
    output:
        peaks = outdir + "/{sample_id}/{sample_id}.fire_peaks.bed",
    log:
        logdir + "/{sample_id}/ft_call_peaks.log"
    threads: 1
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
        max_fdr = config.get("Params", {}).get("fibertools", {}).get("max_fdr", 0.05),
        min_fire_frac = config.get("Params", {}).get("fibertools", {}).get("min_fire_frac", None),
        sd_cov = config.get("Params", {}).get("fibertools", {}).get("sd_cov", 5.0),
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("ft_call_peaks", log_file=log_path)
            rule_logger.info(f"Calling FIRE peaks for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(sample_outdir, f"ft_call_peaks_{current_time}.sh")
            cmd = [
                params.ft, "call-peaks",
                "-o", output.peaks,
                "--sd-cov", str(params.sd_cov),
                "--max-fdr", str(params.max_fdr),
            ]
            if params.min_fire_frac is not None:
                cmd.extend(["--min-fire-frac", str(params.min_fire_frac)])
            cmd.append(input.bam)
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"samtools index {input.bam}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "FIRE peak calling completed for sample {wildcards.sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"ft_call_peaks failed: {e}\n")
            raise RuntimeError(f"ft_call_peaks failed: {e}\n")


rule ft_qc:
    """Collect QC metrics from a Fiber-seq BAM.

    Input: Fiber-seq BAM with m6A, nucleosome, MSP, and FIRE calls.
    Output: TSV file with QC metrics.
    """
    input:
        bam = upstream_outdir + "/{sample_id}/{sample_id}.fiberseq.fire.bam"
    output:
        qc = outdir + "/{sample_id}/{sample_id}.qc.tsv",
    log:
        logdir + "/{sample_id}/ft_qc.log"
    threads: 1
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
        use_acf = config.get("Params", {}).get("fibertools", {}).get("acf", False),
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger("ft_qc", log_file=log_path)
            rule_logger.info(f"Collecting QC metrics for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(sample_outdir, f"ft_qc_{current_time}.sh")
            cmd = [
                params.ft, "qc",
                input.bam,
                output.qc,
            ]
            if params.use_acf:
                cmd.append("--acf")
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "QC metrics collected for sample {wildcards.sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"ft_qc failed: {e}\n")
            raise RuntimeError(f"ft_qc failed: {e}\n")
