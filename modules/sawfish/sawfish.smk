include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
fasta = config.get("genome", {}).get("fasta")
fai = config.get("genome", {}).get("fai_index")
bam_substring = config.get("bam_substring") or ""
call_sv = config.get("Params", {}).get("sawfish", {}).get("call_sv") or False

def get_input_for_sawfish_discover(wildcards):
    logger.info(f"sawfish_discover called with {wildcards}")
    sa = config.get("Params", {}).get("sawfish", {}) or {}
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}." + bam_substring + ".bai")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}/{wildcards.sample_id}.bai")
    in_dict["fasta"] = _require(fasta, "genome.fasta", "sawfish")
    if fai:
        in_dict["fai"] = fai
    if sa.get("expected_cn"):
        in_dict["expected_cn"] = sa.get("expected_cn")
    if sa.get("cnv_excluded_regions"):
        in_dict["cnv_excluded_regions"] = sa.get("cnv_excluded_regions")
    return in_dict

rule sawfish_discover:
    """Per-sample sawfish discover: candidate SV assembly + binned depth
    track + draft CNV segmentation (sawfish user guide 'Analysis steps')."""
    input:
        unpack(get_input_for_sawfish_discover)
    output:
        done = outdir + "/discover/{sample_id}/sawfish_discover.done"
    log:
        logdir + "/sawfish/{sample_id}/sawfish_discover.log"
    benchmark: logdir + "/sawfish/{sample_id}/sawfish_discover.benchmark.tsv",
    threads: rule_res("sawfish_discover").get("threads", 8)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sawfish_discover", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        sawfish = config.get("Procedure", {}).get("sawfish") or "sawfish",
        sa = config.get("Params", {}).get("sawfish", {}) or {},
        out_dir = outdir + "/discover/{sample_id}"
    conda:
        "sawfish.yaml"
    container:
        sif("sawfish.yaml", "sawfish_discover")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sawfish_discover", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sawfish discover for sample {wildcards.sample_id} at {current_time}")
            script = receipt_script(f"sawfish_discover_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, f"discover/{wildcards.sample_id}/sawfish_discover.sh")
            sa = params.sa
            fast_cnv = sa.get("fast_cnv_mode")
            if fast_cnv is None:
                fast_cnv = not (sa.get("call_sv") or False)
            cmd = [
                params.sawfish, "discover",
                "--threads", str(threads),
                "--ref", str(input.fasta),
                "--bam", str(input.bam),
                "--output-dir", str(params.out_dir),
                "--clobber",
            ]
            if getattr(input, "expected_cn", ""):
                cmd += ["--expected-cn", str(input.expected_cn)]
            if getattr(input, "cnv_excluded_regions", ""):
                cmd += ["--cnv-excluded-regions", str(input.cnv_excluded_regions)]
            if fast_cnv:
                cmd += ["--fast-cnv-mode"]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write(f"echo 'sawfish discover for {wildcards.sample_id} was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            with open(str(output.done), "w") as f:
                f.write("ok\n")
        except Exception as e:
            rule_logger.error(f"sawfish discover failed for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"sawfish discover failed for sample {wildcards.sample_id}: {e}, please check the log file {log_path} for more details.")

rule sawfish_joint_call:
    """Cross-sample sawfish joint-call: merge/genotype SV haplotypes, re-run
    CNV segmentation on all samples, unify SV/CNV records into one VCF."""
    input:
        discover_done = expand(outdir + "/discover/{sample_id}/sawfish_discover.done", sample_id=samples),
        fasta = _require(fasta, "genome.fasta", "sawfish"),
        fai = [fai] if fai else []
    output:
        vcf = maybe_temp(outdir + "/joint_call/genotyped.sv.vcf.gz", "sawfish")
    log:
        logdir + "/sawfish/sawfish_joint_call.log"
    benchmark: logdir + "/sawfish/sawfish_joint_call.benchmark.tsv",
    threads: rule_res("sawfish_joint_call").get("threads", 8)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sawfish_joint_call", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    params:
        sawfish = config.get("Procedure", {}).get("sawfish") or "sawfish",
        sa = config.get("Params", {}).get("sawfish", {}) or {},
        out_dir = outdir + "/joint_call"
    conda:
        "sawfish.yaml"
    container:
        sif("sawfish.yaml", "sawfish_joint_call")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sawfish_joint_call", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sawfish joint-call at {current_time}")
            script = receipt_script(f"sawfish_joint_call_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "sawfish_joint_call.sh")
            cmd = [
                params.sawfish, "joint-call",
                "--threads", str(threads),
                "--ref", str(input.fasta),
                "--output-dir", str(params.out_dir),
                "--clobber",
            ]
            sample_dirs = [os.path.dirname(str(d)) for d in input.discover_done]
            if not sample_dirs:
                raise RuntimeError(
                    "sawfish joint-call needs at least one sample: set 'samples' "
                    "in the workflow config JSON (or the module config), then rerun.")
            for d in sample_dirs:
                cmd += ["--sample", d]
            if params.sa.get("treat_single_copy_as_haploid"):
                cmd += ["--treat-single-copy-as-haploid"]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write("echo 'sawfish joint-call was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sawfish joint-call failed: {e}")
            raise RuntimeError(f"sawfish joint-call failed: {e}, please check the log file {log_path} for more details.")

rule sawfish_cnv_filter:
    """Split the joint-call VCF by SVCLAIM into the CNV set (default
    deliverable) and, only when Params.sawfish.call_sv, the breakpoint-SV
    set. Single engine command (bin/sawfish_split_vcf.py)."""
    input:
        vcf = outdir + "/joint_call/genotyped.sv.vcf.gz"
    output:
        cnv_vcf = outdir + "/cnv/sawfish.cnv.vcf.gz",
        cnv_csi = outdir + "/cnv/sawfish.cnv.vcf.gz.csi",
        summary = outdir + "/cnv/sawfish_cnv_summary.tsv",
        sv_vcf = (outdir + "/cnv/sawfish.sv.vcf.gz") if call_sv else [],
        sv_csi = (outdir + "/cnv/sawfish.sv.vcf.gz.csi") if call_sv else []
    log:
        logdir + "/sawfish/sawfish_cnv_filter.log"
    benchmark: logdir + "/sawfish/sawfish_cnv_filter.benchmark.tsv",
    threads: rule_res("sawfish_cnv_filter").get("threads", 2)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "sawfish_cnv_filter", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        python = config.get("Procedure", {}).get("python") or "python3",
        bgzip = config.get("Procedure", {}).get("bgzip") or "bgzip",
        bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools",
        engine = os.path.join(ROOT_DIR, "modules/sawfish/bin/sawfish_split_vcf.py"),
        out_dir = outdir + "/cnv",
        keep_sv = call_sv
    conda:
        "sawfish.yaml"
    container:
        sif("sawfish.yaml", "sawfish_cnv_filter")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sawfish_cnv_filter", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sawfish CNV/SV VCF split at {current_time}")
            script = receipt_script(f"sawfish_cnv_filter_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "sawfish_cnv_filter.sh")
            cmd = [
                params.python, params.engine,
                "--vcf", str(input.vcf),
                "--out-dir", str(params.out_dir),
                "--bgzip", params.bgzip,
                "--bcftools", params.bcftools,
            ]
            if params.keep_sv:
                cmd += ["--keep-sv"]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write("echo 'sawfish cnv_filter was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"sawfish cnv_filter failed: {e}")
            raise RuntimeError(f"sawfish cnv_filter failed: {e}, please check the log file {log_path} for more details.")

rule sawfish_result:
    """Aggregate sawfish CNV (and optional SV) outputs as dependency endpoint."""
    input:
        cnv_vcf = outdir + "/cnv/sawfish.cnv.vcf.gz",
        cnv_csi = outdir + "/cnv/sawfish.cnv.vcf.gz.csi",
        summary = outdir + "/cnv/sawfish_cnv_summary.tsv",
        sv_vcf = (outdir + "/cnv/sawfish.sv.vcf.gz") if call_sv else [],
        sv_csi = (outdir + "/cnv/sawfish.sv.vcf.gz.csi") if call_sv else []
