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
indir = config.get("indir", "input")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""
assembly_dir = config.get("assembly_dir", "")
_inputs = config.get("inputs", {}) or {}
_tel = config.get("Params", {}).get("telomere", {})
species = _tel.get("telogator2", {}).get("species", "human")
_tel_formats = _tel.get("formats") or ["png"]
if isinstance(_tel_formats, str):
    _tel_formats = [_tel_formats]

_TEL_BUILTIN_REF_SPECIES = {"mouse", "maize"}


def _validate_telogator_ref(path: str) -> str:
    if not os.path.exists(path):
        raise ValueError(
            f"telogator2 -t reference '{path}' does not exist; set "
            f"'genome.references.<genome>.subtel_ref' to an existing file.")
    return path


def get_telogator2_ref():
    explicit = config.get("genome", {}).get("subtel_ref") or ""
    if explicit:
        return _validate_telogator_ref(explicit)
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
    if species in _TEL_BUILTIN_REF_SPECIES:
        return f"resources/non-human/telogator-ref-{species}.fa.gz"
    raise ValueError(
        f"telogator2 has no built-in subtelomere reference for species "
        f"'{species}' (shipped non-human refs: {sorted(_TEL_BUILTIN_REF_SPECIES)}); "
        f"set 'genome.references.<genome>.subtel_ref' in the config JSON.")


def _inp(method: str, key: str, sample_id: str):
    d = _inputs.get(method, {})
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
    logger.info(f"telogator2_run called with {wildcards}")
    bam = _inp("telogator2", "bam", wildcards.sample_id) or _shapeA_bam(wildcards.sample_id)
    bai = _inp("telogator2", "bai", wildcards.sample_id) or _shapeA_bam(wildcards.sample_id, True)
    return {"bam": bam, "bai": bai}


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
    return _require(_tel.get("depth_qc", {}).get("telomere_bed"),
                    "Params.telomere.depth_qc.telomere_bed", "telomere")


def get_tidk_scan_input(wildcards):
    in_dict = {"fasta": get_assembly_fasta(wildcards)}
    tidk_db = os.path.expanduser("~/.local/share/tidk/tidk_database.csv")
    if not os.path.exists(tidk_db):
        in_dict["db_init"] = outdir + "/.tidk_build.done"
    return in_dict

# ============================================================
# Rule 1: Telogator2 (existing, per-chromosome-arm)
# ============================================================
rule telogator2_run:
    """Run Telogator2 telomere length analysis on a single sample."""
    input:
        unpack(get_input_for_telogator2)
    output:
        allele_tsv = os.path.join(outdir, "telogator2/{sample_id}/tlens_by_allele.tsv"),
        allele_plot = os.path.join(outdir, "telogator2/{sample_id}/all_final_alleles.png"),
        atl_plot = os.path.join(outdir, "telogator2/{sample_id}/violin_atl.png")
    log:
        logdir + "/{sample_id}/telogator2.log"
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "telogator2_run")
    params:
        telogator2 = config.get("Procedure", {}).get("telogator2") or "telogator2",
        ref = get_telogator2_ref()
    threads: rule_res("telogator2_run").get("threads", 16)
    benchmark: logdir + "/{sample_id}/telogator2.benchmark.tsv",
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
            rule_logger.info(f"Start telogator2 for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"telogator2_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            output_dir = os.path.join(outdir, "telogator2", wildcards.sample_id)
            os.makedirs(output_dir, exist_ok=True)
            plain_script = os.path.join(output_dir, f"telogator2_run_{wildcards.sample_id}.sh")
            cmd = [
                params.telogator2, "-i", input.bam,
                "-o", output_dir,
                "-r", "hifi",
                "-p", str(threads)
            ]
            if params.ref:
                cmd.extend(["-t", params.ref])
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "telogator2_run completed successfully"')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during telogator2_run for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during telogator2_run for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

        # ============================================================
        # Rule 2: Assembly-based telomere scan (Approach A)
        # ============================================================
rule assembly_telomere_scan:
    """Scan hifiasm assembly contig ends for telomeric repeats.

    This is the recommended approach for mouse telomeres (30-150kb).
    Requires assembly from centromere module.
    """
    input:
        fasta = get_assembly_fasta
    output:
        tsv = os.path.join(outdir, "assembly_scan/{sample_id}/{sample_id}_assembly_telomere.tsv"),
        stats = os.path.join(outdir, "assembly_scan/{sample_id}/{sample_id}_assembly_telomere_stats.txt"),
        plots = [os.path.join(outdir, "assembly_scan/{sample_id}/{sample_id}_assembly_telomere_length." + fmt)
                 for fmt in _tel_formats],
    log:
        logdir + "/{sample_id}/assembly_telomere_scan.log"
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "assembly_telomere_scan")
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/scan_assembly_telomere.py"),
        output_dir = os.path.join(outdir, "assembly_scan", "{sample_id}"),
        scan_length = _tel.get("assembly_scan", {}).get("assembly_scan_length") or 50000,
        formats = _tel_formats
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
            rule_logger.info(f"Start assembly telomere scan for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"assembly_telomere_scan_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "assembly_scan", wildcards.sample_id, f"assembly_telomere_scan_{wildcards.sample_id}.sh")
            os.makedirs(os.path.dirname(plain_script), exist_ok=True)
            cmd = [
                "python", params.script,
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
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "assembly_telomere_scan completed successfully"')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during assembly_telomere_scan for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during assembly_telomere_scan for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

        # ============================================================
        # Rule 3: Read-level k-mer density (Approach B)
        # ============================================================
rule read_density_telomere:
    """Estimate genome-wide average telomere length from read-level k-mer density.

    Counts telomeric k-mers across all HiFi reads and divides by chromosome arms.
    Gives a genome-wide average, not per-chromosome measurements.
    """
    input:
        unpack(get_input_for_read_density)
    output:
        tsv = os.path.join(outdir, "read_density/{sample_id}/{sample_id}_read_telomere.tsv"),
        stats = os.path.join(outdir, "read_density/{sample_id}/{sample_id}_read_telomere_stats.txt"),
        plots = [os.path.join(outdir, "read_density/{sample_id}/{sample_id}_read_telomere_content." + fmt)
                 for fmt in _tel_formats],
    log:
        logdir + "/{sample_id}/read_density_telomere.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/read_density_telomere.py"),
        output_dir = os.path.join(outdir, "read_density", "{sample_id}"),
        n_chrom_arms = (_tel.get("read_density") or {}).get("n_chrom_arms") or 40,
        formats = _tel_formats
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
            rule_logger.info(f"Start read density telomere for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"read_density_telomere_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "read_density", wildcards.sample_id, f"read_density_telomere_{wildcards.sample_id}.sh")
            os.makedirs(os.path.dirname(plain_script), exist_ok=True)
            cmd = [
                "python", params.script,
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
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "read_density_telomere completed successfully"')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during read_density_telomere for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during read_density_telomere for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

        # ============================================================
        # Rule 4: tidk scan (Approach C)
        # ============================================================
rule tidk_init:
    """Build tidk reference database (one-time init)."""
    output:
        outdir + "/.tidk_build.done"
    log:
        logdir + "/all/tidk_init.log"
    params:
        tidk = config.get("Procedure", {}).get("tidk") or "tidk"
    threads: rule_res("tidk_init").get("threads", 1)
    benchmark: logdir + "/all/tidk_init.benchmark.tsv",
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
            rule_logger.info("Building tidk reference database (one-time init)")
            script = os.path.join(RECEIPT_DIR, f"tidk_init_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "tidk_init.sh")
            cmd = [params.tidk, "build"]
            cmd_touch = ["touch", output[0]]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(" ".join(cmd_touch) + "\n")
                f.write(f'echo "tidk_init completed successfully"')
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
                 for fmt in _tel_formats],
    log:
        logdir + "/{sample_id}/tidk_scan.log"
    params:
        output_dir = os.path.join(outdir, "tidk", "{sample_id}"),
        tidk = config.get("Procedure", {}).get("tidk") or "tidk",
        motif = "TTAGGG",
        plot_script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/plot_tidk.py"),
        formats = _tel_formats
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
            rule_logger.info(f"Start tidk scan for sample {wildcards.sample_id} at {current_time}")

            script = os.path.join(RECEIPT_DIR, f"tidk_scan_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "tidk", wildcards.sample_id, f"tidk_scan_{wildcards.sample_id}.sh")
            os.makedirs(os.path.dirname(plain_script), exist_ok=True)
            cmd_search = [
                params.tidk, "search",
                "--string", params.motif,
                "--output", f"{wildcards.sample_id}_tidk_search",
                "--dir", params.output_dir,
                input.fasta
            ]
            _search_tsv = os.path.join(
                params.output_dir,
                f"{wildcards.sample_id}_tidk_search_telomeric_repeat_windows.tsv")
            cmd_plot = [
                "python", params.plot_script,
                "-i", _search_tsv,
                "-s", wildcards.sample_id,
                "-o", params.output_dir,
                "--out_base", f"{wildcards.sample_id}_tidk_telomeres",
            ]
            for _f in params.formats:
                cmd_plot += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd_search) + "\n")
                f.write(" ".join(cmd_plot) + "\n")
                f.write(f'echo "tidk_scan completed successfully"')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            search_tsv = os.path.join(params.output_dir, f"{wildcards.sample_id}_tidk_search_telomeric_repeat_windows.tsv")
            if os.path.exists(search_tsv):
                shutil.copy2(search_tsv, output.tsv)
                rule_logger.info(f"Copied tidk search output to {output.tsv}")
            else:
                rule_logger.warning(f"tidk search output not found: {search_tsv}")
                with open(output.tsv, "w") as f:
                    f.write("# No telomeric repeats found by tidk\n")

        except Exception as e:
            rule_logger.error(f"Error occurred during tidk_scan for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during tidk_scan for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

        # ============================================================
        # Result rules (aggregation points for subworkflow)
        # ============================================================
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


# ============================================================
# Phase 3: multi-sample combined figures (all samples in one plot)
#   + cross-method overview (when all methods run)
# ============================================================
def _all_method_tsvs(method, suffix):
    return [os.path.join(outdir, method, s, s + suffix) for s in samples]

rule read_density_combined:
    input:
        tsvs = lambda w: _all_method_tsvs("read_density", "_read_telomere.tsv")
    output:
        plots = [os.path.join(outdir, "read_density", "combined." + fmt)
                 for fmt in _tel_formats],
    log:
        logdir + "/telomere/read_density_combined.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/plot_combined.py"),
        outdir = os.path.join(outdir, "read_density"),
        formats = _tel_formats,
    threads: 1
    conda: "../telomere.yaml"
    run:
        log_path = str(log); open(log_path, "w").close()
        rl = setup_logger("read_density_combined", log_file=log_path)
        try:
            script = os.path.join(RECEIPT_DIR, "read_density_combined.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain = os.path.join(outdir, "read_density", "read_density_combined.sh")
            cmd = ["python", params.script, "-m", "read_density",
                   "-c", "tel_bp_estimate", "-o", params.outdir,
                   "--out_base", "combined"]
            for s in samples:
                cmd += ["-s", f"{s}={os.path.join(params.outdir, s, s + '_read_telomere.tsv')}"]
            for _f in params.formats:
                cmd += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write('echo "read_density_combined completed"\n')
            shutil.copyfile(script, plain)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain} >> {log_path} 2>&1")
        except Exception as e:
            rl.error(f"read_density_combined failed: {e}\n")
            raise RuntimeError(f"read_density_combined failed: {e}, see {log_path}\n")


rule telomere_method_overview:
    input:
        stats = _all_method_tsvs("read_density", "_read_telomere_stats.txt")
    output:
        plots = [os.path.join(outdir, "combined_methods." + fmt) for fmt in _tel_formats],
    log:
        logdir + "/telomere/telomere_method_overview.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/length/bin/plot_combined.py"),
        outdir = outdir,
        formats = _tel_formats,
    threads: 1
    conda: "../telomere.yaml"
    run:
        log_path = str(log); open(log_path, "w").close()
        rl = setup_logger("telomere_method_overview", log_file=log_path)
        try:
            script = os.path.join(RECEIPT_DIR, "telomere_method_overview.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain = os.path.join(outdir, "telomere_method_overview.sh")
            cmd = ["python", params.script, "-m", "overview", "-o", params.outdir,
                   "--out_base", "combined_methods"]
            for s in samples:
                cmd += ["-s", f"read_density={os.path.join(params.outdir, 'read_density', s, s + '_read_telomere_stats.txt')}"]
            for _f in params.formats:
                cmd += ["-F", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write('echo "overview completed"\n')
            shutil.copyfile(script, plain)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain} >> {log_path} 2>&1")
        except Exception as e:
            rl.error(f"telomere_method_overview failed: {e}\n")
            raise RuntimeError(f"overview failed: {e}, see {log_path}\n")


# ============================================================
# Phase 4: depth QC for aligned-BAM methods (⑥)
#   scan reference for telomere coords + count sequencing depth as QC
# ============================================================
rule telomere_depth_qc:
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
        outdir = os.path.join(outdir, "depth_qc", "{sample_id}"),
        control = _tel.get("depth_qc", {}).get("control", ""),
    threads: 1
    conda: "../telomere.yaml"
    run:
        log_path = str(log); open(log_path, "w").close()
        rl = setup_logger("telomere_depth_qc", log_file=log_path)
        try:
            script_p = os.path.join(RECEIPT_DIR, f"telomere_depth_qc_{wildcards.sample_id}.sh")
            os.makedirs(os.path.dirname(script_p), exist_ok=True)
            plain = os.path.join(outdir, "depth_qc", wildcards.sample_id,
                                 f"telomere_depth_qc_{wildcards.sample_id}.sh")
            os.makedirs(os.path.dirname(plain), exist_ok=True)
            cmd = ["python", params.script, "-b", input.bam, "-d", input.bed,
                   "-s", wildcards.sample_id, "-o", params.outdir]
            if params.control:
                cmd += ["-c", params.control]
            with open(script_p, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write('echo "depth_qc completed"\n')
            shutil.copyfile(script_p, plain)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain} >> {log_path} 2>&1")
        except Exception as e:
            rl.error(f"telomere_depth_qc failed: {e}\n")
            raise RuntimeError(f"depth_qc failed: {e}, see {log_path}\n")
