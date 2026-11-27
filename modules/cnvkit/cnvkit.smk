include: "../common/common.smk"
ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
bam_suffix = config.get("bam_suffix", "")
samples = config.get("samples", [])
control_samples = config.get("control_samples", [])
_paired_design_raw = config.get("paired_design") or []
if isinstance(_paired_design_raw, list):
    paired_design = {
        rec["experimental"]: rec["control"]
        for rec in _paired_design_raw
        if rec.get("experimental") and rec.get("control")
    }
else:
    paired_design = dict(_paired_design_raw)
access = config.get("genome", {}).get("access")
targets_bed = config.get("genome", {}).get("targets")
paired_tumors = list(paired_design.keys())


_cnvkit_reference_method = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_reference", {}).get("method", "wgs")
_cnvkit_reference_use_paired = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_reference", {}).get("mode", "auto") == "paired"
_cnvkit_batch_method = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_batch", {}).get("method", "wgs")
_cnvkit_batch_use_paired = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_batch", {}).get("mode", "auto") == "paired"
_cnvkit_paired_method = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_paired", {}).get("method", "wgs")
_cnvkit_paired_use_paired = config.get("Params", {}).get("cnvkit", {}).get("cnvkit_paired", {}).get("mode", "auto") == "paired"

cnv_samples = paired_tumors if _cnvkit_paired_use_paired else list(samples)



rule cnvkit_reference:
    """Build CNVkit reference from control samples."""
    input:
        control_bams = expand(indir + "/{sid}/{sid}" + bam_suffix + ".bam", sid=control_samples),
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "cnvkit"),
        access = lambda wildcards: [access] if access else [],
        targets = lambda wildcards: [targets_bed] if (_cnvkit_reference_method == "hybrid" and targets_bed) else [],
    output:
        ref = [outdir + "/reference/reference.cnn"] if (control_samples and not _cnvkit_reference_use_paired) else []
    log:
        logdir + "/cnvkit/cnvkit_reference.log"
    threads: rule_res("cnvkit_reference").get("threads", 8)
    benchmark: logdir + "/cnvkit/cnvkit_reference.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "cnvkit_reference", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 16000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_reference")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = _cnvkit_reference_method,
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("cnvkit_reference", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit reference build at {current_time}")
            sample_outdir = os.path.dirname(str(output.ref))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"cnvkit_reference_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"cnvkit_reference.sh")
            cmd = [params.cnvkit, "batch"]
            cmd += [str(bam) for bam in input.control_bams]
            cmd += ["--method", params.method, "--fasta", str(input.fasta)]
            if input.access:
                cmd += ["--access", str(input.access[0])]
            if input.targets:
                cmd += ["--targets", str(input.targets[0])]
            cmd += ["--processes", str(threads), "--output-reference", str(output.ref)]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit reference build: {e}")
            raise RuntimeError(f"Error occurred during cnvkit reference build: {e}, you can check the log file {log_path} for more details.")

rule cnvkit_batch:
    """Run CNVkit batch analysis on all samples."""
    input:
        treat_bams = expand(indir + "/{sid}/{sid}" + bam_suffix + ".bam", sid=samples),
        ref = [outdir + "/reference/reference.cnn"] if control_samples else [],
        fasta =  _require(config.get("genome", {}).get("fasta"), "genome.fasta", "cnvkit"),
        access = [access] if access else [],
        targets = [targets_bed] if (_cnvkit_batch_method == "hybrid" and targets_bed) else [],
    output:
        cnr = expand(outdir + "/{sample_id}/{sample_id}.cnr", sample_id=samples) if not _cnvkit_batch_use_paired else [],
        cns = expand(outdir + "/{sample_id}/{sample_id}.cns", sample_id=samples) if not _cnvkit_batch_use_paired else []
    log:
        logdir + "/cnvkit/cnvkit_batch.log"
    threads: rule_res("cnvkit_batch").get("threads", 8)
    benchmark: logdir + "/cnvkit/cnvkit_batch.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "cnvkit_batch", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_batch")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = _cnvkit_batch_method,
        staging = outdir + "/batch",
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("cnvkit_batch", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit batch at {current_time}")
            script = receipt_script(f"cnvkit_batch_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(str(outdir), f"cnvkit_batch.sh")
            os.makedirs(str(outdir), exist_ok=True)
            os.makedirs(params.staging, exist_ok=True)
            cmd = [params.cnvkit, "batch"]
            cmd += [str(bam) for bam in input.treat_bams]
            if input.ref:
                cmd += ["--reference", str(input.ref[0])]
            else:
                cmd += ["-n"]
                cmd += ["--method", params.method, "--fasta", str(input.fasta)]
                if input.access:
                    cmd += ["--access", str(input.access[0])]
                if input.targets:
                    cmd += ["--targets", str(input.targets[0])]
            cmd += ["--processes", str(threads), "--output-dir", params.staging, "--scatter", "--diagram"]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                for sid in samples:
                    f.write(f'mkdir -p "{outdir}/{sid}"\n')
                    f.write(f'for f in "{params.staging}/{sid}{bam_suffix}."* "{params.staging}/{sid}{bam_suffix}-"*; do\n')
                    f.write(f'    if [ -e "$f" ]; then\n')
                    f.write(f'        base="${{f##*/}}"\n')
                    f.write(f'        mv "$f" "{outdir}/{sid}/{sid}${{base#{sid}{bam_suffix}}}"\n')
                    f.write(f'    fi\n')
                    f.write(f'done\n')
                    f.write(f'if [ ! -s "{outdir}/{sid}/{sid}.cnr" ] || [ ! -s "{outdir}/{sid}/{sid}.cns" ]; then\n')
                    f.write(f'    echo "ERROR: declared output {outdir}/{sid}/{sid}.cnr or .cns missing/empty after cnvkit batch" >&2; exit 1\n')
                    f.write(f'fi\n')
                f.write(f'rm -rf "{params.staging}"\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit batch: {e}")
            raise RuntimeError(f"Error occurred during cnvkit batch: {e}, you can check the log file {log_path} for more details.")

def get_paired_normal_bam(wildcards):
    """Resolve the matched normal BAM for a tumor sample (paired mode)."""
    normal_sid = paired_design.get(wildcards.sample_id)
    if not normal_sid:
        raise ValueError(
            f"cnvkit: no matched normal in paired_design for sample "
            f"'{wildcards.sample_id}'; fill Params paired_design or unset it "
            "to use the pooled/flat batch mode."
        )
    return f"{indir}/{normal_sid}/{normal_sid}{bam_suffix}.bam"

rule cnvkit_paired:
    """Run CNVkit batch per tumor/normal pair (paired mode).

    Each tumor is called against its own matched normal: one `cnvkit batch
    <tumor.bam> -n <normal.bam>` invocation per pair (CNVkit pools every -n
    file it receives, so pairs must NOT share one command). Output layout
    matches the pooled mode (outdir/<sample_id>/<sample_id>.{cnr,cns}), so
    cnvkit_batch declares its outputs only when paired_design is empty.
    """
    input:
        tumor_bam = lambda wildcards: f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}{bam_suffix}.bam",
        normal_bam = get_paired_normal_bam,
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "cnvkit"),
        access = lambda wildcards: [access] if access else [],
        targets = lambda wildcards: [targets_bed] if (_cnvkit_paired_method == "hybrid" and targets_bed) else [],
    output:
        cnr = [outdir + "/{sample_id}/{sample_id}.cnr"] if _cnvkit_paired_use_paired else [],
        cns = [outdir + "/{sample_id}/{sample_id}.cns"] if _cnvkit_paired_use_paired else [],
    log:
        logdir + "/{sample_id}/cnvkit_paired.log"
    threads: rule_res("cnvkit_paired").get("threads", 8)
    benchmark: logdir + "/{sample_id}/cnvkit_paired.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "cnvkit_paired", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_paired")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = _cnvkit_paired_method,
        staging = outdir + "/batch_paired/{sample_id}",
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("cnvkit_paired", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit batch (paired) for sample {wildcards.sample_id} at {current_time}")
            script = receipt_script(f"cnvkit_paired_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(str(outdir), wildcards.sample_id, f"cnvkit_paired.sh")
            os.makedirs(os.path.join(str(outdir), wildcards.sample_id), exist_ok=True)
            os.makedirs(params.staging, exist_ok=True)
            cmd = [
                params.cnvkit, "batch",
                str(input.tumor_bam),
                "--normal", str(input.normal_bam),
                "--method", params.method, "--fasta", str(input.fasta),
            ]
            if input.access:
                cmd += ["--access", str(input.access[0])]
            if input.targets:
                cmd += ["--targets", str(input.targets[0])]
            cmd += ["--processes", str(threads), "--output-dir", params.staging, "--scatter", "--diagram"]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                sid = wildcards.sample_id
                f.write(f'for f in "{params.staging}/{sid}{bam_suffix}."* "{params.staging}/{sid}{bam_suffix}-"*; do\n')
                f.write(f'    if [ -e "$f" ]; then\n')
                f.write(f'        base="${{f##*/}}"\n')
                f.write(f'        mv "$f" "{outdir}/{sid}/{sid}${{base#{sid}{bam_suffix}}}"\n')
                f.write(f'    fi\n')
                f.write(f'done\n')
                f.write(f'rm -rf "{params.staging}"\n')
                f.write(f'if [ ! -s "{outdir}/{sid}/{sid}.cnr" ] || [ ! -s "{outdir}/{sid}/{sid}.cns" ]; then\n')
                f.write(f'    echo "ERROR: declared output {outdir}/{sid}/{sid}.cnr or .cns missing/empty after cnvkit batch" >&2; exit 1\n')
                f.write(f'fi\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit batch (paired) for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during cnvkit batch (paired) for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

rule cnvkit_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        cnr = expand(outdir + "/{sample_id}/{sample_id}.cnr", sample_id=cnv_samples),
        cns = expand(outdir + "/{sample_id}/{sample_id}.cns", sample_id=cnv_samples)
