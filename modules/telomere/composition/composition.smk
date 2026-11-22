include: "../../common/common.smk"
import json  # noqa: E402

outdir = config.get("outdir", "output")
logdir = config.get("logdir_combine", config.get("logdir", "log"))
samples = config.get("samples", [])
_inputs = config.get("inputs", {}) or {}
_tel = config.get("Params", {}).get("telomere", {})
_comp = _tel.get("composition", {})
_tel_formats = _tel.get("formats") or ["png"]
if isinstance(_tel_formats, str):
    _tel_formats = [_tel_formats]


def get_composition_bams(wildcards):
    """All samples' fire.sorted.bam (字典模式: inputs.composition.bam[sample]).

    Every configured sample must have an entry — never silently drop a
    sample from the multi-sample analysis; missing entries fail fast naming
    the config key and sample(s).
    """
    per = (_inputs.get("composition", {}) or {}).get("bam", {}) or {}
    missing = [s for s in samples if s not in per]
    if missing:
        raise ValueError(
            f"Missing reference file for telomere composition: set "
            f"'inputs.composition.bam' for sample(s) {missing} in the "
            f"workflow config JSON, then rerun.")
    return [per[s] for s in samples]


rule composition_run:
    """Telomere length + repeat composition (multi-sample, telomere_profile)."""
    input:
        bams = get_composition_bams
    output:
        per_read = os.path.join(outdir, "telomere_per_read.tsv"),
        summary = os.path.join(outdir, "telomere_summary.tsv"),
        composition = os.path.join(outdir, "telomere_composition.tsv"),
        violin = [os.path.join(outdir, "telomere_length_violin." + fmt)
                  for fmt in _tel_formats],
    log:
        logdir + "/telomere/composition.log"
    params:
        script = os.path.join(ROOT_DIR, "modules/telomere/composition/bin/telomere_profile.py"),
        outdir = outdir,
        samples = samples,
        formats = _tel_formats,
        min_read_len = _comp.get("min_read_len") or 5000,
        max_reads = _comp.get("max_reads") or 2_000_000,
    threads: rule_res("composition_run").get("threads", 1)
    conda: "../telomere.yaml"
    container:
        sif("../telomere.yaml", "composition_run")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("composition_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start telomere composition at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"composition_run_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "composition_run.sh")
            per = (_inputs.get("composition", {}) or {}).get("bam", {}) or {}
            cmd = ["python", params.script, "-o", params.outdir,
                   "--min-read-len", str(params.min_read_len),
                   "--max-reads", str(params.max_reads)]
            for s in params.samples:
                if s in per:
                    cmd += ["-s", f"{s}={per[s]}"]
            for _f in params.formats:
                cmd += ["-f", _f]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "composition_run completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during composition_run: {e}\n")
            raise RuntimeError(f"Error occurred during composition_run: {e}, check {log_path}\n")


rule composition_result:
    input:
        per_read = os.path.join(outdir, "telomere_per_read.tsv"),
        summary = os.path.join(outdir, "telomere_summary.tsv"),
        composition = os.path.join(outdir, "telomere_composition.tsv"),
        violin = [os.path.join(outdir, "telomere_length_violin." + fmt)
                  for fmt in _tel_formats],
