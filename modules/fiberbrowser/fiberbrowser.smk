include: "../common/common.smk"
import json  # noqa: E402  (common.smk 未导入而本模块需要的包)

outdir = config.get("outdir", "output/fiberbrowser")
logdir = config.get("logdir", "log/group")
samples = config.get("samples", [])
sample_sams = config.get("sample_sams", {}) or {}
ROOT_DIR = config.get("ROOT_DIR", ".")

loci_map = config.get("Params", {}).get("fiberbrowser", {}).get("fiberbrowser", {}).get("loci_map", {}) or {}
sample_group_dict = config.get("Params", {}).get("fiberbrowser", {}).get("fiberbrowser", {}).get("sample_group_dict", {}) or {}
variants = config.get("Params", {}).get("fiberbrowser", {}).get("fiberbrowser", {}).get("variants", None) or ["plain", "zoom"]
formats = config.get("Params", {}).get("fiberbrowser", {}).get("fiberbrowser", {}).get("formats", None) or ["png"]
if isinstance(formats, str):
    formats = [formats]
loci_map_result = config.get("Params", {}).get("fiberbrowser", {}).get("fiberbrowser", {}).get("loci_map", {}) or {}
variants_result = config.get("Params", {}).get("fiberbrowser", {}).get("fiberbrowser", {}).get("variants", None) or ["plain", "zoom"]
formats_result = config.get("Params", {}).get("fiberbrowser", {}).get("fiberbrowser", {}).get("formats", None) or ["png"]
if isinstance(formats_result, str):
    formats_result = [formats_result]


rule fiberbrowser:
    """Paper-style fiber-seq locus browser figures (Fig 2H scheme)."""
    input:
        sams = sorted(set(sample_sams.values())),
    output:
        plots = [outdir + "/{locus}/{locus}" + ("" if v == "plain" else "_zoom")
                 + "." + fmt for v in variants for fmt in formats],
    wildcard_constraints:
        locus = r"[^/]+"
    log:
        logdir + "/fiberbrowser/fiberbrowser_{locus}.log"
    threads: rule_res("fiberbrowser").get("threads", 1)
    benchmark: logdir + "/fiberbrowser/fiberbrowser_{locus}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberbrowser", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 4000),
    conda:
        "fiberbrowser.yaml"
    container:
        sif("fiberbrowser.yaml", "fiberbrowser")
    params:
        script = ROOT_DIR + "/modules/fiberbrowser/bin/fiberbrowser.py",
        region = lambda wildcards: loci_map[wildcards.locus],
        sample_sams = sample_sams,
        sample_group_dict = sample_group_dict,
        variants = variants,
        formats = formats,
        max_fibers = config.get("Params", {}).get("fiberbrowser", {}).get("max_fibers", 12),
        zoom_bp = config.get("Params", {}).get("fiberbrowser", {}).get("zoom_bp", 2100),
        gtf = config.get("Params", {}).get("fiberbrowser", {}).get("gtf"),
        zoom_center = config.get("Params", {}).get("fiberbrowser", {}).get("zoom_center"),
        zoom_center_mode = config.get("Params", {}).get("fiberbrowser", {}).get("zoom_center_mode"),
        sample_order = config.get("Params", {}).get("fiberbrowser", {}).get("sample_order"),
        dpi = config.get("Params", {}).get("fiberbrowser", {}).get("dpi", 200),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberbrowser", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start fiberbrowser for locus {wildcards.locus} at {current_time}")
            locus_outdir = os.path.dirname(str(output.plots[0]))
            os.makedirs(locus_outdir, exist_ok=True)
            script = receipt_script(f"fiberbrowser_{wildcards.locus}_{current_time}.sh")
            plain_script = os.path.join(locus_outdir, "fiberbrowser.sh")
            cmd = [
                "python", params.script,
                "--locus", f"{params.region}|{wildcards.locus}",
                "--out-dir", locus_outdir,
                "--sample-sams", json.dumps(params.sample_sams),
                "--sample-group-dict", json.dumps(params.sample_group_dict),
                "--max-fibers", str(params.max_fibers),
                "--zoom-bp", str(params.zoom_bp),
                "--dpi", str(params.dpi),
            ]
            if params.gtf:
                cmd += ["--gtf", str(params.gtf)]
            if params.zoom_center is not None:
                cmd += ["--zoom-center", str(params.zoom_center)]
            if params.zoom_center_mode:
                cmd += ["--zoom-center-mode", str(params.zoom_center_mode)]
            if params.sample_order:
                cmd += ["--sample-order", ",".join(params.sample_order)]
            for fmt in params.formats:
                cmd += ["--format", fmt]
            for v in params.variants:
                cmd += ["--variant", v]
            quoted = " ".join(shlex.quote(str(c)) for c in cmd)
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(quoted + "\n")
                f.write(f'echo "fiberbrowser for {wildcards.locus} at '
                        f'{current_time} completed successfully"\n')
            rule_logger.info("Executing: " + quoted)
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(
                f"Error occurred during fiberbrowser for {wildcards.locus}: {e}\n")
            raise RuntimeError(
                f"Error occurred during fiberbrowser for {wildcards.locus}: "
                f"{e}, you can check the log file {log_path} for more details.\n")


rule fiberbrowser_result:
    """Aggregate rule: every locus x variant x format figure."""

    input:
        plots = expand(
            outdir + "/{locus}/{locus}{suffix}.{fmt}",
            locus=list(loci_map_result),
            suffix=["" if v == "plain" else "_zoom" for v in variants_result],
            fmt=formats_result),
