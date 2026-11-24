"""HiFi genome assembly module (hifiasm).

Assembles PacBio HiFi BAM into primary-contig FASTA. Downstream modules
(centromere, telomere assembly_scan/tidk, ...) consume the output file at the
layout contract {outdir}/{sample_id}/assembly/asm.bp.p_ctg.fa.
"""

include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""

def get_input_for_hifiasm(wildcards):
    """Get input files for hifiasm assembly.

    Dict-mode (inputs.hifiasm.bam / inputs.hifiasm.pbi, per-sample) wins —
    e.g. Fiberseq's staged ``{raw}/{sample}/{sample}.bam`` layout; else the
    module's standalone flat shape A ``{indir}/{sample_id}[.bam_substring].
    bam`` (+ its colocated ``.bam.pbi`` index, the same naming contract
    shape A uses).
    """
    per = (config.get("inputs", {}) or {}).get("hifiasm", {}) or {}
    bam_map = per.get("bam", {}) or {}
    pbi_map = per.get("pbi", {}) or {}
    if wildcards.sample_id in bam_map:
        bam = bam_map[wildcards.sample_id]
        return {"bam": bam,
                "bai": pbi_map.get(wildcards.sample_id) or (bam + ".pbi")}
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}." + bam_substring + ".bam.pbi")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}.bam.pbi")
    return in_dict

rule hifiasm_assemble:
    """Assemble HiFi reads with hifiasm."""
    input:
        unpack(get_input_for_hifiasm)
    output:
        fasta = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.fa"
    log:
        logdir + "/{sample_id}/hifiasm_assemble.log"
    benchmark:
        logdir + "/{sample_id}/hifiasm_assemble.benchmark.tsv"
    conda: "hifiasm.yaml"
    container:
        sif("hifiasm.yaml", "hifiasm_assemble")
    threads: rule_res("hifiasm_assemble").get("threads", 32)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "hifiasm_assemble", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 64000),
    params:
        script = os.path.join(ROOT_DIR, "modules/hifiasm/bin/hifiasm_assemble.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        hifiasm = config.get("Procedure", {}).get("hifiasm") or "hifiasm",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
        prefix = outdir + "/{sample_id}/assembly/asm",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("hifiasm_assemble", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hifiasm assembly for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.fasta))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"hifiasm_assemble_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"hifiasm_assemble_{wildcards.sample_id}.sh")
            cmd = [
                params.python, params.script,
                "--bam", str(input.bam),
                "--out-prefix", params.prefix,
                "--out-fasta", str(output.fasta),
                "--threads", str(threads),
                "--samtools", params.samtools,
                "--hifiasm", params.hifiasm,
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
                f.write(f'echo "hifiasm assembly for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + cmd_str)
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during hifiasm assembly for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during hifiasm assembly for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

rule hifiasm_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        fasta = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.fa"
