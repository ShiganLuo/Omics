include: "../../common/common.smk"

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
_inputs = config.get("inputs", {}) or {}
_tel = config.get("Params", {}).get("telomere", {})
_reg = _tel.get("regions", {})


def get_scan_tsv(wildcards):
    """scan_assembly_telomere TSV (dict-mode inputs.regions.scan_tsv[sample],
    else config Params.telomere.regions.scan_tsv).

    Explicit config only — unset config fails fast naming the key to fill
    (never pass None into Snakemake's input: block).
    """
    per = (_inputs.get("regions", {}) or {}).get("scan_tsv", {}) or {}
    return _require(per.get(wildcards.sample_id) or _reg.get("scan_tsv"),
                    "inputs.regions.scan_tsv (per-sample dict) / "
                    "Params.telomere.regions.scan_tsv", "telomere regions")


rule regions_telomere:
    """Build telomere_core/subtelomere BED from an assembly telomere scan TSV."""
    input:
        scan_tsv = get_scan_tsv
    output:
        bed = os.path.join(outdir, "{sample_id}", "{sample_id}.telomere.bed"),
    log:
        logdir + "/{sample_id}/regions_telomere.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/regions/bin/region_build.py"),
        out = os.path.join(outdir, "{sample_id}", "{sample_id}.telomere.bed"),
        subtelomere = _reg.get("subtelomere") or 100000,
        min_len = _reg.get("min_len") or 1000,
    threads: rule_res("regions_telomere").get("threads", 1)
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "regions_telomere")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("regions_telomere", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Build telomere regions for {wildcards.sample_id} at {current_time}")
            script = receipt_script(f"regions_telomere_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/regions_telomere_{wildcards.sample_id}.sh")
            cmd = ["python", params.script, "telomere",
                   "--scan-tsv", input.scan_tsv,
                   "--subtelomere", str(params.subtelomere),
                   "--min-len", str(params.min_len),
                   "-o", params.out]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "regions_telomere completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during regions_telomere for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during regions_telomere: {e}, check {log_path}\n")


rule regions_result:
    input:
        bed = expand(os.path.join(outdir, "{sample_id}", "{sample_id}.telomere.bed"),
                     sample_id=samples),
