include: "../common/common.smk"

outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"
logdir_combine = config.get("logdir_combine") or logdir
ROOT_DIR = config.get("ROOT_DIR", ".")

sample_snv_vcf_dict = config.get("sample_snv_vcf_dict") or {}
sample_sv_vcf_dict = config.get("sample_sv_vcf_dict") or {}
sample_cnv_dict = config.get("sample_cnv_dict") or {}
burden_params = config.get("Params", {}).get("burden", {}) or {}
targets_bed = config.get("genome", {}).get("targets") or None
gtf = config.get("genome", {}).get("gtf") or None
annotate = burden_params.get("annotate", True)
SNV_TSV = outdir + "/{sample_id}/{sample_id}.snv_indel_burden.tsv"
VEP_VCF = outdir + "/{sample_id}/{sample_id}.vep.vcf"
SV_TSV = outdir + "/{sample_id}/{sample_id}.sv_burden.tsv"
CNV_TSV = outdir + "/{sample_id}/{sample_id}.cnv_burden.tsv"
SUMMARY_TSV = outdir + "/burden_summary.tsv"


rule burden_snv_indel:
    """Per-sample TMB and indel burden from a filtered somatic VCF.

    VEP annotation + meaningful-variant filtering on by default (see
    Params.burden.annotate); the annotated VCF is kept beside the TSV.
    """
    input:
        vcf = lambda wildcards: _require(
            sample_snv_vcf_dict.get(wildcards.sample_id),
            "sample_snv_vcf_dict", "burden")
    output:
        tsv = SNV_TSV,
        vep_vcf = [VEP_VCF] if annotate else []
    log:
        logdir + "/{sample_id}/burden_snv_indel.log"
    threads: rule_res("burden_snv_indel").get("threads", 1)
    resources:
        mem_mb = rule_res("burden_snv_indel").get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/burden/bin/burden.py"),
        min_vaf = burden_params.get("min_vaf", 0.05),
        size_mb = burden_params.get("size_mb", 34.0),
        gtf_feature = burden_params.get("gtf_feature", "exon"),
        include_sex = burden_params.get("include_sex_chroms", False),
        annotate = annotate,
        vep_cache = burden_params.get("vep_cache") or "~/.vep",
        vep_species = burden_params.get("vep_species") or "mus_musculus",
        vep_assembly = burden_params.get("vep_assembly") or "GRCm39",
        targets_bed = targets_bed,
        gtf = gtf,
    conda:
        "burden.yaml"
    container:
        sif("burden.yaml", "burden_snv_indel")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("burden_snv_indel", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start burden snv-indel for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"burden_snv_indel_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"burden_snv_indel_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script, "snv-indel",
                "-i", str(input.vcf),
                "-s", wildcards.sample_id,
                "-o", str(output.tsv),
                "-v", str(params.min_vaf),
                "-m", str(params.size_mb),
                "-e", str(params.gtf_feature),
                "-c", str(params.vep_cache),
                "-p", str(params.vep_species),
                "-y", str(params.vep_assembly),
                "-w",
            ]
            if not params.annotate:
                cmd += ["--no-annotate"]
            if params.targets_bed:
                cmd += ["-b", str(params.targets_bed)]
            elif params.gtf:
                cmd += ["-g", str(params.gtf)]
            if params.include_sex:
                cmd += ["-x"]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "burden snv-indel for {wildcards.sample_id} completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during burden snv-indel for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during burden snv-indel for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule burden_sv:
    """Per-sample structural-variant burden from a somatic SV VCF."""
    input:
        vcf = lambda wildcards: _require(
            sample_sv_vcf_dict.get(wildcards.sample_id),
            "sample_sv_vcf_dict", "burden")
    output:
        tsv = SV_TSV
    log:
        logdir + "/{sample_id}/burden_sv.log"
    threads: rule_res("burden_sv").get("threads", 1)
    resources:
        mem_mb = rule_res("burden_sv").get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/burden/bin/burden.py"),
    conda:
        "burden.yaml"
    container:
        sif("burden.yaml", "burden_sv")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("burden_sv", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start burden sv for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"burden_sv_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"burden_sv_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script, "sv",
                "-i", str(input.vcf),
                "-s", wildcards.sample_id,
                "-o", str(output.tsv),
                "-w",
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "burden sv for {wildcards.sample_id} completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during burden sv for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during burden sv for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule burden_cnv:
    """Per-sample copy-number burden from a CNVkit .cns segment table."""
    input:
        cns = lambda wildcards: _require(
            sample_cnv_dict.get(wildcards.sample_id),
            "sample_cnv_dict", "burden")
    output:
        tsv = CNV_TSV
    log:
        logdir + "/{sample_id}/burden_cnv.log"
    threads: rule_res("burden_cnv").get("threads", 1)
    resources:
        mem_mb = rule_res("burden_cnv").get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/burden/bin/burden.py"),
        min_abs_log2 = burden_params.get("min_abs_log2", 0.3),
    conda:
        "burden.yaml"
    container:
        sif("burden.yaml", "burden_cnv")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("burden_cnv", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start burden cnv for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"burden_cnv_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"burden_cnv_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script, "cnv",
                "-i", str(input.cns),
                "-s", wildcards.sample_id,
                "-o", str(output.tsv),
                "-m", str(params.min_abs_log2),
                "-w",
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "burden cnv for {wildcards.sample_id} completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during burden cnv for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during burden cnv for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule burden_aggregate:
    """Merge the per-sample burden tables into one cross-sample summary TSV."""
    input:
        snv_tsvs = [SNV_TSV.format(sample_id=s) for s in sample_snv_vcf_dict],
        sv_tsvs = [SV_TSV.format(sample_id=s) for s in sample_sv_vcf_dict],
        cnv_tsvs = [CNV_TSV.format(sample_id=s) for s in sample_cnv_dict],
    output:
        summary = SUMMARY_TSV
    log:
        logdir_combine + "/burden/burden_aggregate.log"
    threads: rule_res("burden_aggregate").get("threads", 1)
    resources:
        mem_mb = rule_res("burden_aggregate").get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/burden/bin/burden.py"),
    conda:
        "burden.yaml"
    container:
        sif("burden.yaml", "burden_aggregate")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("burden_aggregate", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start burden aggregate at {current_time}")
            os.makedirs(os.path.dirname(str(output.summary)), exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"burden_aggregate_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(str(output.summary)), "burden_aggregate.sh")
            # One flag carries every value per table type (nargs-style option).
            cmd = ["python", params.script, "aggregate", "-o", str(output.summary), "-w"]
            if input.snv_tsvs:
                cmd += ["-n"] + [str(p) for p in input.snv_tsvs]
            if input.sv_tsvs:
                cmd += ["-v"] + [str(p) for p in input.sv_tsvs]
            if input.cnv_tsvs:
                cmd += ["-c"] + [str(p) for p in input.cnv_tsvs]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write('echo "burden aggregate completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during burden aggregate: {e}\n")
            raise RuntimeError(f"Error occurred during burden aggregate: {e}, you can check the log file {log_path} for more details.\n")


rule burden_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        snv = [SNV_TSV.format(sample_id=s) for s in sample_snv_vcf_dict],
        sv = [SV_TSV.format(sample_id=s) for s in sample_sv_vcf_dict],
        cnv = [CNV_TSV.format(sample_id=s) for s in sample_cnv_dict],
        summary = SUMMARY_TSV
    resources:
        mem_mb = rule_res("burden_result").get("mem_mb", 8000),
