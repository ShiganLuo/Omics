include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
bam_suffix = config.get("bam_suffix", "")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
_derive = config.get("Params", {}).get("derive", {}) or {}

rule derive_targets:
    """Derive the capture-region BED from aligned BAM coverage.

    Upstream rule: runs right after markdup so exome/hybrid assays get one
    interval BED consumed by both GATK -L and CNVkit --targets. Vendor
    panels are often unavailable or on an older assembly build; the region
    is reconstructed from actual aligned reads on the CURRENT reference,
    so build mismatch cannot occur.
    """
    input:
        bams = expand(indir + "/{sid}/{sid}" + bam_suffix + ".bam", sid=samples),
    output:
        bed = outdir + "/targets.bed",
    log:
        logdir + "/derive_targets.log",
    threads: rule_res("derive_targets").get("threads", 4),
    resources:
        mem_mb = rule_res("derive_targets").get("mem_mb", 16000),
    conda:
        "../samtools/samtools.yaml"
    container:
        sif("../samtools/samtools.yaml", "derive_targets")
    params:
        script = os.path.join(ROOT_DIR, "modules/targets/bin/derive_targets.py"),
        min_depth = _derive.get("min_depth", 10),
        min_bq = _derive.get("min_bq", 20),
        min_mq = _derive.get("min_mq", 20),
        min_gap = _derive.get("min_gap", 50),
        min_region = _derive.get("min_region", 50),
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("derive_targets", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start derive_targets at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"derive_targets_{current_time}.sh")
            plain_script = os.path.join(outdir, "derive_targets.sh")
            os.makedirs(outdir, exist_ok=True)
            cmd = [
                "python", params.script,
                "-b", *[str(bam) for bam in input.bams],
                "-o", str(output.bed),
                "--min-depth", str(params.min_depth),
                "--min-bq", str(params.min_bq),
                "--min-mq", str(params.min_mq),
                "--min-gap", str(params.min_gap),
                "--min-region", str(params.min_region),
            ]
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during derive_targets: {e}")
            raise RuntimeError(f"Error occurred during derive_targets: {e}, you can check the log file {log_path} for more details.")
