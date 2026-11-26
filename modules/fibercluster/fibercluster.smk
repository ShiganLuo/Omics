include: "../common/common.smk"

# ---------------------------------------------------------------------------
# fibercluster: cluster single Fiber-seq fibers by accessibility pattern at
# target regions (one engine run per channel/k setting) and curate the runs
# into a reader-facing report package (canonical single-molecule figures +
# cross-locus cluster-proportion overview + sensitivity appendix).
#
# Group-level module: one run consumes all samples (no {sample_id} wildcard);
# the {setting} wildcard names one channel/k configuration (e.g. fire_k2).
# Region names live inside the regions BED (data, not wildcards), so the
# per-setting outputs are tracked by a done marker (modules.md 输出声明:
# 禁用 directory(); 标记文件推荐). Engine CLIs consume config-native shapes;
# this file is wiring only.
# ---------------------------------------------------------------------------
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "common/4_elements")
samples = config.get("samples", [])

fibercluster_params = config.get("Params", {}).get("fibercluster", {})
fc_params = fibercluster_params.get("fiber_cluster", {})
fr_params = fibercluster_params.get("fiber_report", {})
# Settings labels = the {setting} wildcard values (one engine run each).
# Zero-config default: no `settings` and no per-setting overrides -> one run
# labeled "auto" with engine defaults (-k auto: silhouette-chosen per region).
# Sensitivity presets stay available: list labels here AND give each label an
# entry in channel_by_setting/k_by_setting only when it deviates from
# Params.fibercluster.fiber_cluster.{channel,k} (the module-level defaults).
_settings = config.get("settings") or ["auto"]
if isinstance(_settings, str):
    _settings = [_settings]


def _fiber_cluster_inputs(wildcards):
    """Declare every file one fiber_cluster run consumes.

    Shape A per modules.md: per-sample element paths spliced from indir +
    sample_id + suffix. nuc/fire overlays join only when their suffix keys
    are configured (optional-file pattern #1), matching the engine's
    optional --nuc/--fire inputs.

    Returns a dict for ``unpack``.
    """
    files = {
        "regions": _require(config.get("regions_bed"), "regions_bed", tool="fibercluster"),
        "msp": [os.path.join(indir, s, s + (fc_params.get("msp_suffix") or ".msp.bed.gz"))
                for s in samples],
    }
    nuc_suffix = fc_params.get("nuc_suffix")
    if nuc_suffix:
        files["nuc"] = [os.path.join(indir, s, s + nuc_suffix) for s in samples]
    fire_suffix = fc_params.get("fire_suffix")
    if fire_suffix:
        files["fire"] = [os.path.join(indir, s, s + fire_suffix) for s in samples]
    return files


def _fiber_report_inputs(wildcards):
    """Report consumes every setting's marker plus the same regions BED."""
    return {
        "regions": _require(config.get("regions_bed"), "regions_bed", tool="fibercluster"),
        "markers": [outdir + f"/{label}/fiber_cluster.done" for label in _settings],
    }


rule fiber_cluster:
    """Cluster single fibers per region for one channel/k setting.

    Per-region outputs (tables + single-molecule figure) live beside the
    done marker inside {outdir}/{setting}/<region>/; region names come from
    the regions BED at run time, so the marker is the tracked output.
    """
    input:
        unpack(_fiber_cluster_inputs),
    output:
        done = outdir + "/{setting}/fiber_cluster.done",
    log:
        logdir + "/fibercluster/{setting}.log"
    threads: rule_res("fiber_cluster").get("threads", 4)
    benchmark: logdir + "/fibercluster/{setting}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiber_cluster", input_size=input_size_named(input, 'msp'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq/fiberseq.yaml"
    container:
        sif("../fiberseq/fiberseq.yaml", "fiber_cluster")
    params:
        script = os.path.join(ROOT_DIR, "modules", "fibercluster", "bin", "fiber_cluster.py"),
        p = fc_params,
        samples = samples,
        sif_path = sif("../fiberseq/fiberseq.yaml", "fiber_cluster"),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiber_cluster", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start fiber_cluster for setting {wildcards.setting} at {current_time}")
            setting_outdir = os.path.dirname(str(output.done))
            # 规则重跑 = 全量重算;先清上一轮产物,避免脏文件掩蔽失败检测
            # (如 regions_failed.txt 残留会把 rc=1 误判为「部分区域失败」)
            shutil.rmtree(setting_outdir, ignore_errors=True)
            os.makedirs(setting_outdir, exist_ok=True)
            script = receipt_script(f"fiber_cluster_{wildcards.setting}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # 设置名可能含路径分隔符
            plain_script = os.path.join(setting_outdir, f"fiber_cluster_{wildcards.setting}.sh")
            s_label = str(wildcards.setting)
            # channel/k fall back to the module-level defaults; per-setting
            # dicts override ONLY the labels they name (sensitivity presets).
            s_channel = (params.p.get("channel_by_setting") or {}).get(s_label) \
                or params.p.get("channel") or "msp"
            s_k = (params.p.get("k_by_setting") or {}).get(s_label, params.p.get("k", 0))
            s_k = int(s_k) if s_k else 0  # 0 = engine auto-k
            # 回执脚本经 subprocess.call 直跑(不走 Snakemake shell(),否则会被
            # 再包一层容器执行且无法捕获退出码),故命令须显式 apptainer exec 经
            # SIF 执行:引擎依赖 pysam(仅 fiberseq.sif 内有)。--bind /data/16T
            # 覆盖脚本/元素/输出全部路径。
            cmd = [
                "apptainer", "exec", "--bind", "/data/16T", params.sif_path,
                "python", params.script,
                "-r", str(input.regions),
                "-o", setting_outdir,
                "-C", str(s_channel),
                "-k", str(s_k),
                "-S", ",".join(params.samples),
                "-w",  # 规则重跑即全量重算(mtime 由 done 标记跟踪)
            ]
            # optional knobs: forward only what the config actually sets, so
            # the engine's own defaults stay the single source of truth
            _opt_flags = {
                "bin_size": "-b", "bin_mode": "-B", "min_msp_bp": "-P",
                "method": "-e", "min_fibers": "-n", "min_span_frac": "-s",
                "bin_threshold": "-t", "max_fibers": "-a",
                "right_tracks": "-T", "fire_fdr": "-q", "seed": "-x",
            }
            for key, flag in _opt_flags.items():
                val = params.p.get(key)
                if val is not None and val != "":
                    cmd += [flag, str(val)]
            cache_dir = params.p.get("cache_dir")
            if cache_dir:
                cmd += ["-y", str(cache_dir)]
            for s_id, msp_path in zip(params.samples, input.msp):
                cmd += ["-m", f"{s_id}={msp_path}"]
            nuc_inputs = getattr(input, "nuc", [])
            for s_id, nuc_path in zip(params.samples, nuc_inputs):
                cmd += ["-u", f"{s_id}={nuc_path}"]
            fire_inputs = getattr(input, "fire", [])
            for s_id, fire_path in zip(params.samples, fire_inputs):
                cmd += ["-f", f"{s_id}={fire_path}"]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "fiber_cluster for {wildcards.setting} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(str(x) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                import glob as _glob
                import subprocess as _sp
                rc = _sp.call(f"bash {plain_script} >> {log_path} 2>&1", shell=True)
                failed_list = os.path.join(setting_outdir, "regions_failed.txt")
                if rc == 1 and os.path.exists(failed_list):
                    # 引擎契约(README §1.1):区域级失败(如 n_fibers<10)收进
                    # regions_failed.txt 且退出码 1 = 部分区域失败,不阻断 DAG;
                    # 失败清单随逐区域产物一并被跟踪。
                    rule_logger.warning(
                        f"fiber_cluster {wildcards.setting}: 部分区域失败(rc=1),"
                        f"清单见 {failed_list}")
                elif rc != 0:
                    raise RuntimeError(
                        f"fiber_cluster {wildcards.setting} failed with rc={rc}, "
                        f"see {log_path}")
                if not _glob.glob(os.path.join(setting_outdir, "*", "cluster_assignments.tsv")):
                    raise RuntimeError(
                        f"fiber_cluster {wildcards.setting}: 无任何区域产出(全部失败?), "
                        f"see {log_path}")
                # 标记文件: 区域名来自 regions BED,逐区域产物随标记被 DAG 跟踪
                open(str(output.done), "w").close()
            rule_logger.info(f"fiber_cluster completed for setting {wildcards.setting}")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiber_cluster for {wildcards.setting}: {e}\n")
            raise RuntimeError(f"Error occurred during fiber_cluster for {wildcards.setting}: {e}, you can check the log file {log_path} for more details.\n")


rule fiber_report:
    """Curate all setting runs into figures/ + summary/ + sensitivity/."""
    input:
        unpack(_fiber_report_inputs),
    output:
        loci_summary = outdir + "/summary/loci_summary.tsv",
        proportions = outdir + "/summary/cluster_proportions.png",
        sensitivity = outdir + "/sensitivity/sensitivity_table.tsv",
        figures = outdir + "/figures.done",
    log:
        logdir + "/fibercluster/fiber_report.log"
    threads: 1
    benchmark: logdir + "/fibercluster/fiber_report.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiber_report", attempt=attempt).get("mem_mb", 4000),
    conda:
        "../fiberseq/fiberseq.yaml"
    container:
        sif("../fiberseq/fiberseq.yaml", "fiber_report")
    params:
        script = os.path.join(ROOT_DIR, "modules", "fibercluster", "bin", "fiber_report.py"),
        p = fr_params,
        settings = list(_settings),
        samples = samples,
        sif_path = sif("../fiberseq/fiberseq.yaml", "fiber_report"),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiber_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start fiber_report at {current_time}")
            report_outdir = os.path.dirname(str(output.loci_summary))
            os.makedirs(report_outdir, exist_ok=True)
            script = receipt_script(f"fiber_report_{current_time}.sh")
            plain_script = os.path.join(report_outdir, "fiber_report.sh")
            canonical = params.p.get("canonical")
            if not canonical:
                # default: the FIRST setting is the deliverable run; a
                # single-run setup needs no canonical key at all
                canonical = params.settings[0]
                rule_logger.info(f"fiber_report.canonical not set -> using first setting {canonical!r}")
            cmd = ["apptainer", "exec", "--bind", "/data/16T", params.sif_path,
                   "python", params.script, "-c", str(canonical), "-o", str(outdir)]
            for label in params.settings:
                cmd += ["-r", f"{label}={outdir}/{label}"]
            cmd += ["-S", ",".join(params.samples)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "fiber_report was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(str(x) for x in cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                import subprocess as _sp
                # shell() 会被 Snakemake 包进容器执行,但无法捕获退出码(rc=1 契约);
                # 这里用 subprocess.call 直跑回执脚本,故命令显式 apptainer exec 经
                # SIF 执行(见 cmd 构造处)。
                rc = _sp.call(f"bash {plain_script} >> {log_path} 2>&1", shell=True)
                if rc != 0:
                    raise RuntimeError(f"fiber_report failed with rc={rc}, see {log_path}")
                # figures/ 目录名随 regions BED 变化,figures.done 兜底跟踪
                open(str(output.figures), "w").close()
            rule_logger.info("fiber_report completed")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiber_report: {e}\n")
            raise RuntimeError(f"Error occurred during fiber_report: {e}, you can check the log file {log_path} for more details.\n")


rule fibercluster_result:
    """Panel rule: aggregate this module's deliverables as the node.py
    registration point (modules.md 收尾规则)."""

    input:
        markers = expand(outdir + "/{setting}/fiber_cluster.done", setting=_settings),
        loci_summary = outdir + "/summary/loci_summary.tsv",
        proportions = outdir + "/summary/cluster_proportions.png",
        sensitivity = outdir + "/sensitivity/sensitivity_table.tsv",
        figures = outdir + "/figures.done",
