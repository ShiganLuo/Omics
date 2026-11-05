include: "../common/common.smk"
ROOT_DIR = config.get("ROOT_DIR", ".")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
bam_suffix = config.get("bam_suffix", "")
samples = config.get("samples", [])
control_samples = config.get("control_samples", [])
# CNVkit calling mode (Params.cnvkit.mode, resolved by node.py for "auto"):
#   paired - per tumor/normal pair, each tumor called against its own matched
#            normal (rule cnvkit_paired)
#   pooled - shared reference built from control_samples (rule cnvkit_reference
#            + cnvkit_batch --reference)
#   flat   - no reference at all (bare -n; uniform expected coverage)
_mode_cfg = config.get("Params", {}).get("cnvkit", {}).get("mode") or "auto"
_paired_design_raw = config.get("paired_design") or []
if isinstance(_paired_design_raw, list):
    paired_design = {
        rec["experimental"]: rec["control"]
        for rec in _paired_design_raw
        if rec.get("experimental") and rec.get("control")
    }
else:
    paired_design = dict(_paired_design_raw)
if _mode_cfg not in ("auto", "paired", "pooled", "flat"):
    raise ValueError(
        f"cnvkit: unknown Params.cnvkit.mode '{_mode_cfg}'; "
        "expected one of: auto, paired, pooled, flat."
    )
if _mode_cfg == "paired" and not paired_design:
    raise ValueError(
        "cnvkit: Params.cnvkit.mode 'paired' but paired_design is empty; "
        "fill it (or let node.py derive it from the meta design)."
    )
use_paired = (_mode_cfg == "paired")
access = config.get("genome", {}).get("access")
targets_bed = config.get("genome", {}).get("targets")
_method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs"
paired_tumors = list(paired_design.keys())
cnv_samples = paired_tumors if use_paired else list(samples)



rule cnvkit_reference:
    """Build CNVkit reference from control samples."""
    input:
        control_bams = expand(indir + "/{sid}/{sid}" + bam_suffix + ".bam", sid=control_samples),
        fasta = lambda wildcards: _require(config.get("genome", {}).get("fasta"), "genome.fasta", "cnvkit"),
        access = lambda wildcards: [access] if access else [],
        targets = lambda wildcards: [targets_bed] if (_method == "hybrid" and targets_bed) else [],
    output:
        ref = [outdir + "/reference/reference.cnn"] if (control_samples and not use_paired) else []
    log:
        logdir + "/cnvkit/cnvkit_reference.log"
    threads: rule_res("cnvkit_reference").get("threads", 8)
    resources:
        mem_mb = rule_res("cnvkit_reference").get("mem_mb", 16000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_reference")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs",
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("cnvkit_reference", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit reference build at {current_time}")
            sample_outdir = os.path.dirname(str(output.ref))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"cnvkit_reference_{current_time}.sh")
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
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit reference build: {e}")
            raise RuntimeError(f"Error occurred during cnvkit reference build: {e}, you can check the log file {log_path} for more details.")

rule cnvkit_batch:
    """Run CNVkit batch analysis on all samples."""
    input:
        treat_bams = expand(indir + "/{sid}/{sid}" + bam_suffix + ".bam", sid=samples),
        ref = outdir + "/reference/reference.cnn" if control_samples else "",
        fasta =  _require(config.get("genome", {}).get("fasta"), "genome.fasta", "cnvkit"),
        access = access if access else "",
        targets = targets_bed if (_method == "hybrid" and targets_bed) else "",
    output:
        cnr = expand(outdir + "/{sample_id}/{sample_id}.cnr", sample_id=samples) if not use_paired else [],
        cns = expand(outdir + "/{sample_id}/{sample_id}.cns", sample_id=samples) if not use_paired else []
    log:
        logdir + "/cnvkit/cnvkit_batch.log"
    threads: rule_res("cnvkit_batch").get("threads", 8)
    resources:
        mem_mb = rule_res("cnvkit_batch").get("mem_mb", 32000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_batch")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs",
        staging = outdir + "/batch",
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("cnvkit_batch", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit batch at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"cnvkit_batch_{current_time}.sh")
            plain_script = os.path.join(str(outdir), f"cnvkit_batch.sh")
            os.makedirs(str(outdir), exist_ok=True)
            os.makedirs(params.staging, exist_ok=True)
            cmd = [params.cnvkit, "batch"]
            cmd += [str(bam) for bam in input.treat_bams]
            if input.ref:
                cmd += ["--reference", str(input.ref)]
            else:
                # No control samples: pass a bare -n so CNVkit builds a FLAT
                # reference (uniform expected coverage). Without it, batch
                # refuses to run: "Option -n/--normal must be given to build
                # a new reference if -r/--reference is not used."
                cmd += ["-n"]
                cmd += ["--method", params.method, "--fasta", str(input.fasta)]
                if input.access:
                    cmd += ["--access", str(input.access)]
                if input.targets:
                    cmd += ["--targets", str(input.targets)]
            cmd += ["--processes", str(threads), "--output-dir", params.staging, "--scatter", "--diagram"]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                # Relocate each sample's flat artifacts from the staging dir
                # into the declared per-sample dirs (outdir/<sample_id>/) and
                # STRIP bam_suffix: CNVkit names every artifact after the
                # input BAM basename (<sid><bam_suffix>.cnr, plots
                # <sid><bam_suffix>-scatter.png), while the declared outputs
                # are <sid>.cnr / <sid>.cns.
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
                # Drop staging dir remnants now that artifacts are relocated.
                f.write(f'rm -rf "{params.staging}"\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
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
        targets = lambda wildcards: [targets_bed] if (_method == "hybrid" and targets_bed) else [],
    output:
        # [] when disabled (pooled/flat mode) — empty STRING is job-build-fatal.
        cnr = [outdir + "/{sample_id}/{sample_id}.cnr"] if use_paired else [],
        cns = [outdir + "/{sample_id}/{sample_id}.cns"] if use_paired else [],
    log:
        logdir + "/{sample_id}/cnvkit_paired.log"
    threads: rule_res("cnvkit_paired").get("threads", 8)
    resources:
        mem_mb = rule_res("cnvkit_paired").get("mem_mb", 32000),
    conda:
        "cnvkit.yaml"
    container:
        sif("cnvkit.yaml", "cnvkit_paired")
    params:
        cnvkit = config.get("Procedure", {}).get("cnvkit") or "cnvkit.py",
        method = config.get("Params", {}).get("cnvkit", {}).get("method") or "wgs",
        staging = outdir + "/batch_paired/{sample_id}",
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("cnvkit_paired", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cnvkit batch (paired) for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"cnvkit_paired_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(str(outdir), wildcards.sample_id, f"cnvkit_paired.sh")
            os.makedirs(os.path.join(str(outdir), wildcards.sample_id), exist_ok=True)
            os.makedirs(params.staging, exist_ok=True)
            cmd = [
                params.cnvkit, "batch",
                str(input.tumor_bam),
                # Matched normal: single -n file, CNVkit builds this pair's
                # own reference from it (pooled reference must NOT be used
                # here: -r and -n are mutually exclusive in batch).
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
                # Relocate the pair's flat artifacts into outdir/<sample_id>/
                # AND strip bam_suffix from every name: CNVkit names all
                # artifacts after the input BAM basename
                # (<sid><bam_suffix>.cnr / <sid><bam_suffix>-scatter.png),
                # while the declared and node.py-registered outputs are
                # <sid>.cnr / <sid>.cns. A plain mv keeps the suffix and the
                # job then dies on "Missing output files" after a batch run
                # that logged nothing but success.
                sid = wildcards.sample_id
                f.write(f'for f in "{params.staging}/{sid}{bam_suffix}."* "{params.staging}/{sid}{bam_suffix}-"*; do\n')
                f.write(f'    if [ -e "$f" ]; then\n')
                f.write(f'        base="${{f##*/}}"\n')
                f.write(f'        mv "$f" "{outdir}/{sid}/{sid}${{base#{sid}{bam_suffix}}}"\n')
                f.write(f'    fi\n')
                f.write(f'done\n')
                f.write(f'rm -rf "{params.staging}"\n')
                # The batch exit code alone does not prove the declared
                # outputs exist; fail loudly here instead of leaving it to
                # Snakemake's cause-less "Missing output files".
                f.write(f'if [ ! -s "{outdir}/{sid}/{sid}.cnr" ] || [ ! -s "{outdir}/{sid}/{sid}.cns" ]; then\n')
                f.write(f'    echo "ERROR: declared output {outdir}/{sid}/{sid}.cnr or .cns missing/empty after cnvkit batch" >&2; exit 1\n')
                f.write(f'fi\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during cnvkit batch (paired) for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during cnvkit batch (paired) for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")

rule cnvkit_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        cnr = expand(outdir + "/{sample_id}/{sample_id}.cnr", sample_id=cnv_samples),
        cns = expand(outdir + "/{sample_id}/{sample_id}.cns", sample_id=cnv_samples)
    resources:
        mem_mb = rule_res("cnvkit_result").get("mem_mb", 8000),
