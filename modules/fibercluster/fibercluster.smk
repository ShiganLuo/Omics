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
# Settings labels = the {setting} wildcard values (one engine run each; the
# sigtest `contrasts` shape). Per-label channel/k overrides live with their
# consumer in Params.fibercluster.fiber_cluster.{channel,k}_by_setting
# (the feature_col_by_contrast shape); a plain string is accepted and
# expanded to a one-element list.
_settings = config.get("settings") or ["default"]
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
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiber_cluster", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start fiber_cluster for setting {wildcards.setting} at {current_time}")
            setting_outdir = os.path.dirname(str(output.done))
            os.makedirs(setting_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiber_cluster_{wildcards.setting}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # 设置名可能含路径分隔符
            plain_script = os.path.join(setting_outdir, f"fiber_cluster_{wildcards.setting}.sh")
            s_label = str(wildcards.setting)
            # channel/k are properties OF the setting: each label must carry
            # an entry in channel_by_setting/k_by_setting (no shared global —
            # a missing entry fails naming the JSON key to fill)
            s_channel = (params.p.get("channel_by_setting") or {}).get(s_label)
            if s_channel is None:
                raise RuntimeError(
                    f"Params.fibercluster.fiber_cluster.channel_by_setting has no entry "
                    f"for setting '{s_label}' — add it (msp|fire).")
            s_k = (params.p.get("k_by_setting") or {}).get(s_label)
            if s_k is None:
                raise RuntimeError(
                    f"Params.fibercluster.fiber_cluster.k_by_setting has no entry "
                    f"for setting '{s_label}' — add it (int).")
            cmd = [
                "python", params.script,
                "-r", str(input.regions),
                "-o", setting_outdir,
                "-b", str(100 if params.p.get("bin_size") is None else params.p.get("bin_size")),
                "-B", str(params.p.get("bin_mode") or "any"),
                "-P", str(10 if params.p.get("min_msp_bp") is None else params.p.get("min_msp_bp")),
                "-e", str(params.p.get("method") or "hierarchical"),
                "-C", str(s_channel),
                "-k", str(s_k),
                "-n", str(10 if params.p.get("min_fibers") is None else params.p.get("min_fibers")),
                "-s", str(0.8 if params.p.get("min_span_frac") is None else params.p.get("min_span_frac")),
                "-t", str(0.5 if params.p.get("bin_threshold") is None else params.p.get("bin_threshold")),
                "-a", str(300 if params.p.get("max_fibers") is None else params.p.get("max_fibers")),
                "-T", str(params.p.get("right_tracks") or "fire,lnk,nuc"),
                "-q", str(0.1 if params.p.get("fire_fdr") is None else params.p.get("fire_fdr")),
                "-x", str(0 if params.p.get("seed") is None else params.p.get("seed")),
                # display order = the samples list order (its own ordering is
                # the single source; no separate sample_order key)
                "-S", ",".join(params.samples),
                "-w",  # 规则重跑即全量重算(mtime 由 done 标记跟踪)
            ]
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
                shell(f"bash {plain_script} >> {log_path} 2>&1")
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
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiber_report", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start fiber_report at {current_time}")
            report_outdir = os.path.dirname(str(output.loci_summary))
            os.makedirs(report_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiber_report_{current_time}.sh")
            plain_script = os.path.join(report_outdir, "fiber_report.sh")
            canonical = params.p.get("canonical")
            if not canonical:
                raise RuntimeError(
                    "Params.fibercluster.fiber_report.canonical is empty — name the "
                    "setting label whose figures/summary become the deliverable.")
            cmd = ["python", params.script, "-c", str(canonical), "-o", str(outdir)]
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
                shell(f"bash {plain_script} >> {log_path} 2>&1")
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
