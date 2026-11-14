include: "../common/common.smk"
import json  # noqa: E402  (common.smk 未导入而本模块需要的包)

outdir = config.get("outdir", "output/fiberbrowser")
logdir = config.get("logdir", "log/group")
samples = config.get("samples", [])
sample_sams = config.get("sample_sams", {}) or {}
ROOT_DIR = config.get("ROOT_DIR", ".")

# 指令派生名(变量名 = config 原名)。loci_map / sample_group_dict 由
# node.py 解析 Params.fiberbrowser.loci 后写回 config;配置形状校验归
# schema 与 node.py,模块对配置只取值。
fiberbrowser_params = config.get("Params", {}).get("fiberbrowser", {})
loci_map = fiberbrowser_params.get("loci_map") or {}
sample_group_dict = fiberbrowser_params.get("sample_group_dict") or {}
variants = fiberbrowser_params.get("variants") or ["plain", "zoom"]
formats = fiberbrowser_params.get("formats") or ["png"]


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
    resources:
        mem_mb = rule_res("fiberbrowser").get("mem_mb", 4000),
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
        max_fibers = fiberbrowser_params.get("max_fibers", 12),
        zoom_bp = fiberbrowser_params.get("zoom_bp", 2100),
        dpi = fiberbrowser_params.get("dpi", 200),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberbrowser", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start fiberbrowser for locus {wildcards.locus} at {current_time}")
            locus_outdir = os.path.dirname(str(output.plots[0]))
            os.makedirs(locus_outdir, exist_ok=True)
            script = os.path.join(
                RECEIPT_DIR,
                f"fiberbrowser_{wildcards.locus}_{current_time}.sh")
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
            locus=list(loci_map),
            suffix=["" if v == "plain" else "_zoom" for v in variants],
            fmt=formats),
