"""Telomere analysis module.

This module provides multiple approaches for telomere length measurement:

1. telogator2_run: Per-chromosome-arm telomere length (seed-extend matching)
2. assembly_telomere_scan: Scan assembly contig ends for telomeric repeats (Approach A)
3. read_density_telomere: Genome-wide average from read-level k-mer density (Approach B)
4. tidk_scan: Assembly-based telomere scan using tidk community tool (Approach C)

Note: For mouse telomeres (30-150kb), Approach A is recommended as HiFi reads
(~15-25kb) cannot span the full telomere for telogator2.
"""

include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", logdir)
indir = config.get("indir", "input")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""
assembly_dir = config.get("assembly_dir", "")
inputs = config.get("inputs", {}) or {}
formats = config.get("Params", {}).get("telomere", {}).get("formats") or ["png"]


def get_telogator2_ref():
    """Resolve the telogator2 subtelomere reference FASTA input.

    Prefers the explicit config path, then the tool's bundled non-human
    reference, then the repo-relative builtin path for shipped species.
    Returns "" when telogator2's built-in human reference applies.

    Returns:
        Reference FASTA path, or "" for the built-in human reference.

    Raises:
        ValueError: no reference resolves for the configured species.
    """
    explicit = config.get("genome", {}).get("subtel_ref") or ""
    if explicit:
        if not os.path.exists(explicit):
            raise ValueError(
                f"telogator2 -t reference '{explicit}' does not exist; set "
                f"'genome.references.<genome>.subtel_ref' to an existing file.")
        return explicit
    species = (config.get("Params", {}).get("telomere", {})
               .get("telogator2", {}).get("species", "human"))
    if species == "human":
        return ""
    try:
        result = subprocess.run(["which", "telogator2"], capture_output=True, text=True)
        telogator2_bin = result.stdout.strip()
        if telogator2_bin:
            pyver = f"python{'.'.join(map(str, __import__('sys').version_info[:2]))}"
            resource_dir = os.path.join(os.path.dirname(telogator2_bin), "..",
                "lib", pyver, "site-packages", "source", "resources")
            ref_file = os.path.join(resource_dir, "non-human", f"telogator-ref-{species}.fa.gz")
            if os.path.exists(ref_file):
                return ref_file
    except Exception as e:
        logger.warning(f"Error finding telogator2 ref: {e}")
    builtin_species = {"mouse", "maize"}
    if species in builtin_species:
        return f"resources/non-human/telogator-ref-{species}.fa.gz"
    raise ValueError(
        f"telogator2 has no built-in subtelomere reference for species "
        f"'{species}' (shipped non-human refs: {sorted(builtin_species)}); "
        f"set 'genome.references.<genome>.subtel_ref' in the config JSON.")


def _inp(method: str, key: str, sample_id: str):
    d = inputs.get(method, {})
    if isinstance(d, dict):
        k = d.get(key, {})
        if isinstance(k, dict):
            return k.get(sample_id)
    return None


def _shapeA_bam(sample_id: str, want_bai: bool = False) -> str:
    base = (os.path.join(indir, f"{sample_id}/{sample_id}." + bam_substring)
            if bam_substring != "" else os.path.join(indir, f"{sample_id}/{sample_id}"))
    return base + (".bai" if want_bai else ".bam")


def get_input_for_telogator2(wildcards):
    bam = _inp("telogator2", "bam", wildcards.sample_id) or _shapeA_bam(wildcards.sample_id)
    bai = _inp("telogator2", "bai", wildcards.sample_id) or _shapeA_bam(wildcards.sample_id, True)
    in_dict = {"bam": bam, "bai": bai}
    ref = get_telogator2_ref()
    if ref:
        in_dict["ref"] = ref
    return in_dict


def get_input_for_read_density(wildcards):
    bam = _inp("read_density", "bam", wildcards.sample_id) or _shapeA_bam(wildcards.sample_id)
    return {"bam": bam}


def get_assembly_fasta(wildcards):
    fa = (_inp("assembly_scan", "fasta", wildcards.sample_id)
          or _inp("tidk", "fasta", wildcards.sample_id))
    if fa:
        return fa
    if assembly_dir:
        return os.path.join(assembly_dir, f"{wildcards.sample_id}/assembly/asm.bp.p_ctg.fa")
    return _require(None, "inputs.assembly_scan.fasta (per-sample dict) / "
                    "inputs.tidk.fasta / assembly_dir", "telomere")


def get_depth_qc_bed(wildcards):
    return _require(config.get("Params", {}).get("telomere", {})
                    .get("depth_qc", {}).get("telomere_bed"),
                    "Params.telomere.depth_qc.telomere_bed", "telomere")


def get_tidk_scan_input(wildcards):
    in_dict = {"fasta": get_assembly_fasta(wildcards)}
    tidk_db = os.path.expanduser("~/.local/share/tidk/tidk_database.csv")
    if not os.path.exists(tidk_db):
        in_dict["db_init"] = outdir + "/.tidk_build.done"
    return in_dict


def _all_method_tsvs(method, suffix):
    return [os.path.join(outdir, method, s, s + suffix) for s in samples]


rule telogator2_run:
    """Run Telogator2 telomere length analysis on a single sample."""
    input:
        unpack(get_input_for_telogator2)
    output:
        allele_tsv = os.path.join(outdir, "telogator2/{sample_id}/tlens_by_allele.tsv"),
        allele_plot = os.path.join(outdir, "telogator2/{sample_id}/all_final_alleles.png"),
        atl_plot = os.path.join(outdir, "telogator2/{sample_id}/violin_atl.png")
    log:
        logdir + "/{sample_id}/telogator2_run.log"
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "telogator2_run")
    params:
        telogator2 = config.get("Procedure", {}).get("telogator2") or "telogator2",
    threads: rule_res("telogator2_run").get("threads", 16)
    benchmark: logdir + "/{sample_id}/telogator2_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "telogator2_run", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("telogator2_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start telogator2 run for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.allele_tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"telogator2_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"telogator2_run_{wildcards.sample_id}.sh")
            cmd = [
                params.telogator2, "-i", input.bam,
                "-o", sample_outdir,
                "-r", "hifi",
                "-p", str(threads)
            ]
            if getattr(input, "ref", ""):
                cmd.extend(["-t", input.ref])
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "telogator2_run for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during telogator2_run for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during telogator2_run for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule assembly_telomere_scan:
    """Scan hifiasm assembly contig ends for telomeric repeats.

    This is the recommended approach for mouse telomeres (30-150kb).
    Requires assembly from centromere module.
    """
    input:
        fasta = get_assembly_fasta
    output:
        # intermediate controlled by Params.telomere.keep_intermediates (modules.md temp spec)
        tsv = maybe_temp(os.path.join(outdir, "assembly_scan/{sample_id}/{sample_id}_assembly_telomere.tsv"), "telomere"),
        stats = os.path.join(outdir, "assembly_scan/{sample_id}/{sample_id}_assembly_telomere_stats.txt"),
        plots = [os.path.join(outdir, "assembly_scan/{sample_id}/{sample_id}_assembly_telomere_length." + fmt)
                 for fmt in formats],
    log:
        logdir + "/{sample_id}/assembly_telomere_scan.log"
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "assembly_telomere_scan")
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/scan_assembly_telomere.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        output_dir = os.path.join(outdir, "assembly_scan", "{sample_id}"),
        scan_length = (config.get("Params", {}).get("telomere", {})
                       .get("assembly_scan", {}).get("assembly_scan_length") or 50000),
        formats = formats
    threads: rule_res("assembly_telomere_scan").get("threads", 1)
    benchmark: logdir + "/{sample_id}/assembly_telomere_scan.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "assembly_telomere_scan", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("assembly_telomere_scan", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start assembly_telomere_scan for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"assembly_telomere_scan_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"assembly_telomere_scan_{wildcards.sample_id}.sh")
            cmd = [
                params.python, params.script,
                "--fasta", input.fasta,
                "--sample_name", wildcards.sample_id,
                "--output_dir", params.output_dir,
                "--scan_length", str(params.scan_length)
            ]
            for _f in params.formats:
                cmd += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "assembly_telomere_scan for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during assembly_telomere_scan for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during assembly_telomere_scan for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule read_density_telomere:
    """Estimate genome-wide average telomere length from read-level k-mer density.

    Counts telomeric k-mers across all HiFi reads and divides by chromosome arms.
    Gives a genome-wide average, not per-chromosome measurements.
    """
    input:
        unpack(get_input_for_read_density)
    output:
        # intermediate controlled by Params.telomere.keep_intermediates (modules.md temp spec)
        tsv = maybe_temp(os.path.join(outdir, "read_density/{sample_id}/{sample_id}_read_telomere.tsv"), "telomere"),
        stats = os.path.join(outdir, "read_density/{sample_id}/{sample_id}_read_telomere_stats.txt"),
        plots = [os.path.join(outdir, "read_density/{sample_id}/{sample_id}_read_telomere_content." + fmt)
                 for fmt in formats],
    log:
        logdir + "/{sample_id}/read_density_telomere.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/read_density_telomere.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        output_dir = os.path.join(outdir, "read_density", "{sample_id}"),
        n_chrom_arms = (config.get("Params", {}).get("telomere", {})
                        .get("read_density") or {}).get("n_chrom_arms") or 40,
        formats = formats
    threads: rule_res("read_density_telomere").get("threads", 1)
    benchmark: logdir + "/{sample_id}/read_density_telomere.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "read_density_telomere", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "read_density_telomere")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("read_density_telomere", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start read_density_telomere for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"read_density_telomere_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"read_density_telomere_{wildcards.sample_id}.sh")
            cmd = [
                params.python, params.script,
                "--bam", input.bam,
                "--sample_name", wildcards.sample_id,
                "--output_dir", params.output_dir,
                "--n_chrom_arms", str(params.n_chrom_arms)
            ]
            for _f in params.formats:
                cmd += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "read_density_telomere for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during read_density_telomere for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during read_density_telomere for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule tidk_init:
    """Build tidk reference database (one-time init)."""
    output:
        outdir + "/.tidk_build.done"
    log:
        logdir_combine + "/telomere/tidk_init.log"
    params:
        tidk = config.get("Procedure", {}).get("tidk") or "tidk"
    threads: rule_res("tidk_init").get("threads", 1)
    benchmark: logdir_combine + "/telomere/tidk_init.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "tidk_init", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "tidk_init")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("tidk_init", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start tidk_init at {current_time}")
            sample_outdir = os.path.dirname(str(output[0]))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"tidk_init_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, "tidk_init.sh")
            cmd = [params.tidk, "build"]
            cmd_touch = ["touch", output[0]]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd_touch) + "\n")
                f.write(f'echo "tidk_init was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd + cmd_touch))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during tidk_init: {e}\n")
            raise RuntimeError(f"Error occurred during tidk_init: {e}, you can check the log file {log_path} for more details.\n")


rule tidk_scan:
    """Scan assembly for telomeric repeats using tidk community tool.

    tidk is a standardized tool for telomere identification from assemblies.
    """
    input:
        unpack(get_tidk_scan_input)
    output:
        tsv = os.path.join(outdir, "tidk/{sample_id}/{sample_id}_tidk_telomeres.tsv"),
        plots = [os.path.join(outdir, "tidk/{sample_id}/{sample_id}_tidk_telomeres." + fmt)
                 for fmt in formats],
    log:
        logdir + "/{sample_id}/tidk_scan.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/tidk_scan.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        output_dir = os.path.join(outdir, "tidk", "{sample_id}"),
        tidk = config.get("Procedure", {}).get("tidk") or "tidk",
        motif = "TTAGGG",
        out_base = "{sample_id}_tidk_telomeres",
        formats = formats
    threads: rule_res("tidk_scan").get("threads", 1)
    benchmark: logdir + "/{sample_id}/tidk_scan.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "tidk_scan", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "tidk_scan")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("tidk_scan", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start tidk_scan for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"tidk_scan_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"tidk_scan_{wildcards.sample_id}.sh")
            cmd = [
                params.python, params.script,
                "--fasta", input.fasta,
                "--sample", wildcards.sample_id,
                "--outdir", params.output_dir,
                "--out_base", params.out_base,
                "--tidk", params.tidk,
                "--motif", params.motif,
            ]
            for _f in params.formats:
                cmd += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "tidk_scan for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during tidk_scan for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during tidk_scan for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule telogator2_result:
    """Result aggregation for telogator2."""
    input:
        tsv = outdir + "/telogator2/{sample_id}/tlens_by_allele.tsv"
rule assembly_telomere_result:
    """Result aggregation for assembly-based telomere scan."""
    input:
        stats = outdir + "/assembly_scan/{sample_id}/{sample_id}_assembly_telomere_stats.txt"
rule read_density_telomere_result:
    """Result aggregation for read-density telomere estimation."""
    input:
        stats = outdir + "/read_density/{sample_id}/{sample_id}_read_telomere_stats.txt"
rule tidk_result:
    """Result aggregation for tidk scan."""
    input:
        tsv = outdir + "/tidk/{sample_id}/{sample_id}_tidk_telomeres.tsv"


rule read_density_combined:
    """Combine per-sample read-density estimates into one figure."""
    input:
        tsvs = lambda w: _all_method_tsvs("read_density", "_read_telomere.tsv")
    output:
        plots = [os.path.join(outdir, "read_density", "combined." + fmt)
                 for fmt in formats],
    log:
        logdir_combine + "/telomere/read_density_combined.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/plot_combined.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        output_dir = os.path.join(outdir, "read_density"),
        formats = formats,
    threads: 1
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "read_density_combined")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("read_density_combined", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start read_density_combined at {current_time}")
            sample_outdir = os.path.dirname(str(output.plots[0]))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"read_density_combined_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, "read_density_combined.sh")
            cmd = [params.python, params.script, "-m", "read_density",
                   "-c", "tel_bp_estimate", "-o", params.output_dir,
                   "--out_base", "combined"]
            for tsv in input.tsvs:
                sample = os.path.basename(os.path.dirname(str(tsv)))
                cmd += ["-s", f"{sample}={tsv}"]
            for _f in params.formats:
                cmd += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "read_density_combined was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during read_density_combined: {e}\n")
            raise RuntimeError(f"Error occurred during read_density_combined: {e}, you can check the log file {log_path} for more details.\n")


rule telomere_method_overview:
    """Cross-method overview figure over all samples."""
    input:
        stats = _all_method_tsvs("read_density", "_read_telomere_stats.txt")
    output:
        plots = [os.path.join(outdir, "combined_methods." + fmt) for fmt in formats],
    log:
        logdir_combine + "/telomere/telomere_method_overview.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/plot_combined.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        output_dir = outdir,
        formats = formats,
    threads: 1
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "telomere_method_overview")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("telomere_method_overview", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start telomere_method_overview at {current_time}")
            sample_outdir = os.path.dirname(str(output.plots[0]))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"telomere_method_overview_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, "telomere_method_overview.sh")
            cmd = [params.python, params.script, "-m", "overview", "-o", params.output_dir,
                   "--out_base", "combined_methods"]
            for st in input.stats:
                cmd += ["-s", f"read_density={st}"]
            for _f in params.formats:
                cmd += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "telomere_method_overview was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during telomere_method_overview: {e}\n")
            raise RuntimeError(f"Error occurred during telomere_method_overview: {e}, you can check the log file {log_path} for more details.\n")


rule telomere_depth_qc:
    """Count sequencing depth over telomere regions as QC for aligned-BAM methods."""
    input:
        bam = lambda w: _inp("telogator2", "bam", w.sample_id) or _shapeA_bam(w.sample_id),
        bed = get_depth_qc_bed,
    output:
        tsv = os.path.join(outdir, "depth_qc", "{sample_id}", "{sample_id}_telomere_region_depth.tsv"),
        stats = os.path.join(outdir, "depth_qc", "{sample_id}", "{sample_id}_telomere_region_depth_stats.txt"),
    log:
        logdir + "/{sample_id}/telomere_depth_qc.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/telomere_region_depth.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        output_dir = os.path.join(outdir, "depth_qc", "{sample_id}"),
        control = (config.get("Params", {}).get("telomere", {})
                   .get("depth_qc", {}).get("control", "")),
    threads: 1
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "telomere_depth_qc")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("telomere_depth_qc", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start telomere_depth_qc for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"telomere_depth_qc_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"telomere_depth_qc_{wildcards.sample_id}.sh")
            cmd = [params.python, params.script, "-b", input.bam, "-d", input.bed,
                   "-s", wildcards.sample_id, "-o", params.output_dir]
            if params.control:
                cmd += ["-c", params.control]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "telomere_depth_qc for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during telomere_depth_qc for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during telomere_depth_qc for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")
