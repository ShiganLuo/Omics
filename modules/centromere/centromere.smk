"""Centromere analysis module.

This module performs:
1. RepeatMasker for satellite DNA annotation
2. Centromere statistics extraction

Assembly comes from the upstream hifiasm module (modules/hifiasm): the input
FASTA is expected at the shared layout {assembly_dir}/{sample_id}/assembly/
asm.bp.p_ctg.fa. Note: This module includes common.smk for shared utilities.
"""

include: "../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", logdir)
samples = config.get("samples", [])
assembly_dir = config.get("assembly_dir") or outdir
RepeatMasker_params = config.get("Params", {}).get("RepeatMasker", {})

def get_input_for_repeatmasker_init(wildcards):
    """Declare repeatmasker_init inputs.

    The custom Dfam h5 joins only when configured (optional file, pattern 1);
    without it the engine warms the cache with RepeatMasker's built-in
    libraries. Configured-but-missing fails at DAG build via the input
    declaration.
    """
    files = {}
    if RepeatMasker_params.get("dfam_h5"):
        files["dfam_h5"] = RepeatMasker_params.get("dfam_h5")
    return files

def get_input_for_repeatmasker(wildcards):
    """Declare RepeatMasker inputs.

    The lib-init marker joins only when a custom Dfam h5 is configured: its
    producer (repeatmasker_init) then runs once before all samples; without a
    custom library RepeatMasker uses its built-in defaults and no cache init
    is needed.
    """
    files = {
        "fasta": os.path.join(assembly_dir, f"{wildcards.sample_id}/assembly/asm.bp.p_ctg.fa"),
    }
    if RepeatMasker_params.get("dfam_h5"):
        files["lib_init"] = os.path.join(outdir, ".repeatmasker_lib_init.done")
    return files

rule repeatmasker_init:
    """Build RepeatMasker library cache once to avoid concurrent makeblastdb conflicts."""
    input:
        unpack(get_input_for_repeatmasker_init)
    output:
        done = outdir + "/.repeatmasker_lib_init.done"
    log:
        logdir_combine + "/centromere/repeatmasker_init.log"
    benchmark:
        logdir_combine + "/centromere/repeatmasker_init.benchmark.tsv"
    conda: "centromere.yaml"
    container:
        sif("centromere.yaml", "repeatmasker_init")
    threads: rule_res("repeatmasker_init").get("threads", 1)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "repeatmasker_init", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/centromere/bin/repeatmasker_init.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        RepeatMasker = config.get("Procedure", {}).get("RepeatMasker") or "RepeatMasker",
        species = RepeatMasker_params.get("repeatmasker_init", {}).get("species") or RepeatMasker_params.get("species"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("repeatmasker_init", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start RepeatMasker library cache init at {current_time}")
            sample_outdir = os.path.dirname(str(output.done))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"repeatmasker_init_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "repeatmasker_init.sh")
            cmd = [
                params.python, params.script,
                "--repeatmasker", params.RepeatMasker,
                "--species", params.species or "",
                "--threads", str(threads),
                "--marker", str(output.done),
            ]
            dfam_h5 = getattr(input, "dfam_h5", "")
            if dfam_h5:
                cmd += ["--dfam-h5", str(dfam_h5)]
            cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(cmd_str + "\n")
                f.write('echo "RepeatMasker library cache init was completed successfully"\n')
            rule_logger.info("Executing: " + cmd_str)
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during RepeatMasker library cache init: {e}\n")
            raise RuntimeError(f"Error occurred during RepeatMasker library cache init: {e}, you can check the log file {log_path} for more details.\n")

rule repeatmasker_run:
    """Run RepeatMasker on assembled contigs to annotate satellite DNA."""
    input:
        unpack(get_input_for_repeatmasker)
    output:
        out = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.out",
        tbl = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.tbl"
    log:
        logdir + "/{sample_id}/repeatmasker_run.log"
    benchmark:
        logdir + "/{sample_id}/repeatmasker_run.benchmark.tsv"
    conda: "centromere.yaml"
    container:
        sif("centromere.yaml", "repeatmasker_run")
    threads: rule_res("repeatmasker_run").get("threads", 24)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "repeatmasker_run", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    params:
        rm_dir = outdir + "/{sample_id}/RepeatMasker",
        RepeatMasker = config.get("Procedure", {}).get("RepeatMasker") or "RepeatMasker",
        species = RepeatMasker_params.get("repeatmasker_run", {}).get("species") or RepeatMasker_params.get("species"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("repeatmasker_run", log_file=log_path)
        try:
            if not params.species:
                raise ValueError("Missing RepeatMasker species: set 'Params.RepeatMasker.species' in the workflow config JSON (never guess a species).")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start RepeatMasker for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.out))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"repeatmasker_run_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"repeatmasker_run_{wildcards.sample_id}.sh")
            cmd = [
                params.RepeatMasker,
                "-species", params.species,
                "-pa", str(threads),
                "-dir", params.rm_dir,
                "-gff",
                str(input.fasta),
            ]
            cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(cmd_str + "\n")
                f.write(f'echo "RepeatMasker annotation for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + cmd_str)
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during RepeatMasker for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during RepeatMasker for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

rule centromere_extract:
    """Extract satellite DNA statistics from RepeatMasker output."""
    input:
        out = outdir + "/{sample_id}/RepeatMasker/asm.bp.p_ctg.fa.out"
    output:
        stats = outdir + "/{sample_id}/{sample_id}.centromere_stats.txt"
    log:
        logdir + "/{sample_id}/centromere_extract.log"
    benchmark:
        logdir + "/{sample_id}/centromere_extract.benchmark.tsv"
    conda: "centromere.yaml"
    container:
        sif("centromere.yaml", "centromere_extract")
    threads: rule_res("centromere_extract").get("threads", 1)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "centromere_extract", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/centromere/bin/extract_centromere_stats.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        species = RepeatMasker_params.get("repeatmasker_run", {}).get("species") or RepeatMasker_params.get("species"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("centromere_extract", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start centromere stats extraction for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.stats))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"centromere_extract_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"centromere_extract_{wildcards.sample_id}.sh")
            cmd = [
                params.python, params.script,
                "--rm_out", str(input.out),
                "--output", str(output.stats),
                "--species", params.species or "",
            ]
            cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(cmd_str + "\n")
                f.write(f'echo "centromere stats extraction for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + cmd_str)
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during centromere stats extraction for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during centromere stats extraction for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

rule centromere_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        stats = outdir + "/{sample_id}/{sample_id}.centromere_stats.txt"
