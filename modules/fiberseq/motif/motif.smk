include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/motif")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
_formats_profile = config.get("Params", {}).get("fiberseq_motif", {}).get("fiberseq_motif_profile", {}).get("formats") or ["pdf"]
if isinstance(_formats_profile, str):
    _formats_profile = [_formats_profile]
_formats_result = config.get("Params", {}).get("fiberseq_motif", {}).get("fiberseq_motif_result", {}).get("formats") or ["pdf"]
if isinstance(_formats_result, str):
    _formats_result = [_formats_result]

rule fiberseq_motif_scan:
    """Scan a MEME motif file in FIRE peak sequences (any species)."""
    input:
        peaks = indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        fasta = lambda wildcards: _require(genome_cfg.get("fasta"), "genome.fasta", "fiberseq_motif"),
        motif_file = lambda wildcards: _require(genome_cfg.get("motif_file"), "genome.motif_file", "fiberseq_motif"),
    output:
        sites = outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
    log:
        logdir + "/{sample_id}/fiberseq_motif_scan.log"
    threads: rule_res("fiberseq_motif").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_motif_scan.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_motif", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 4000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_motif_scan")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/motif_scan.py",
        pwm_threshold = config.get("Params", {}).get("fiberseq_motif", {}).get("fiberseq_motif_scan", {}).get("pwm_threshold", 0.8),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_motif_scan", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start motif scan for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.sites))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"fiberseq_motif_scan_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"fiberseq_motif_scan_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.motif_file),
                "-r", str(input.peaks),
                "-f", str(input.fasta),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-t", str(params.pwm_threshold),
                "--threads", str(threads),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_motif_scan for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_motif_scan for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_motif_scan for {wildcards.sample_id}: {e}\n")

rule fiberseq_motif_profile:
    """Aggregate m6A/accessibility/nucleosome profiles around motif sites."""
    input:
        sites = outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        profile = outdir + "/{sample_id}/{sample_id}.motif_profile.tsv",
        scores = outdir + "/{sample_id}/{sample_id}.motif_footprint_scores.tsv",
        plots = [outdir + "/{sample_id}/{sample_id}.motif_profile." + fmt
                 for fmt in _formats_profile],
    log:
        logdir + "/{sample_id}/fiberseq_motif_profile.log"
    threads: rule_res("fiberseq_motif").get("threads", 1)
    benchmark: logdir + "/{sample_id}/fiberseq_motif_profile.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_motif", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_motif_profile")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/motif_profile.py",
        flank = config.get("Params", {}).get("fiberseq_motif", {}).get("fiberseq_motif_profile", {}).get("flank", 250),
        min_fibers = config.get("Params", {}).get("fiberseq_motif", {}).get("fiberseq_motif_profile", {}).get("min_fibers", 10),
        formats = _formats_profile,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_motif_profile", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start motif profile for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.profile))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"fiberseq_motif_profile_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"fiberseq_motif_profile_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.sites),
                "-a", str(input.m6a),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-f", str(params.flank),
                "-i", str(params.min_fibers),
                "-l", log_path,
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_motif_profile for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_motif_profile for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_motif_profile for {wildcards.sample_id}: {e}\n")

rule fiberseq_motif_result:
    input:
        sites = expand(outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
                       sample_id=samples),
        profile = expand(outdir + "/{sample_id}/{sample_id}.motif_profile.tsv",
                         sample_id=samples),
        plots = expand(outdir + "/{sample_id}/{sample_id}.motif_profile.{fmt}",
                       sample_id=samples, fmt=_formats_result),
