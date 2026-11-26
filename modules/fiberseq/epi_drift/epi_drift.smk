include: "../../common/common.smk"


outdir = config.get("outdir", "fiberseq_analysis/epi_drift")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
_inputs = config.get("inputs", {}).get("epi_drift", {})
_bams = _inputs.get("bam") or {}
_p = config.get("Params", {}).get("fiberseq_epi_drift", {})
_groups = _p.get("groups") or {}
_glabels = list(_groups)
_run_p = _p.get("fiberseq_epi_drift_run", {})
_formats_fig = _p.get("fiberseq_epi_drift_fig", {}).get("formats") or ["png"]
if isinstance(_formats_fig, str):
    _formats_fig = [_formats_fig]


def _group_bam_list(wildcards):
    out = []
    for label in _glabels:
        for s in _groups.get(label) or []:
            bam = _bams.get(s)
            if bam and bam not in out:
                out.append(bam)
    return out


def _group_args():
    args = []
    for label in _glabels:
        args += ["--group", f"{label}={','.join(_groups.get(label) or [])}"]
    for label in _glabels:
        for s in _groups.get(label) or []:
            if s in _bams:
                args += ["--bam", f"{s}={_bams[s]}"]
    return args


rule fiberseq_epi_drift_run:
    """Per-fiber m5CpG fraction variance over random windows for two groups
    (fig9 data): epigenetic_drift_windows.tsv + epigenetic_drift.json."""
    input:
        fai = lambda wildcards: _require(genome_cfg.get("fai_index"),
                                         "genome.fai_index",
                                         "fiberseq_epi_drift"),
        bams = _group_bam_list,
    output:
        windows_tsv = outdir + "/epigenetic_drift_windows.tsv",
        summary_json = outdir + "/epigenetic_drift.json",
    log:
        logdir + "/epi_drift_run.log"
    threads: rule_res("fiberseq_epi_drift").get("threads", 4)
    benchmark: logdir + "/epi_drift_run.benchmark.tsv"
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_epi_drift", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_epi_drift_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/epi_drift/bin/epi_drift_run.py",
        group_args = lambda wildcards: _group_args(),
        win = _run_p.get("win", 5000),
        min_fibers = _run_p.get("min_fibers", 20),
        subsample = _run_p.get("subsample", 20),
        n_windows = _run_p.get("n_windows", 600),
        seed = _run_p.get("seed", 42),
        min_chrom_len = _run_p.get("min_chrom_len", 30000000),
        chrom_prefix = _run_p.get("chrom_prefix", ""),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_epi_drift_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start epi_drift run at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            if not _glabels or len(_glabels) != 2:
                raise ValueError(
                    "Params.fiberseq_epi_drift.groups must define exactly "
                    f"2 groups (label -> [samples]), got: {_groups!r}")
            for label in _glabels:
                missing = [s for s in (_groups.get(label) or [])
                           if s not in _bams]
                if missing:
                    raise ValueError(
                        f"group {label!r} references samples without a BAM "
                        f"entry in inputs.epi_drift.bam: {missing}")
            script = receipt_script(f"fiberseq_epi_drift_run_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "fiberseq_epi_drift_run.sh")
            cmd = [
                "python", str(params.script),
                "--fai", str(input.fai),
                "--windows-out", str(output.windows_tsv),
                "--summary-out", str(output.summary_json),
                "--win", str(params.win),
                "--min-fibers", str(params.min_fibers),
                "--subsample", str(params.subsample),
                "--n-windows", str(params.n_windows),
                "--seed", str(params.seed),
                "--min-chrom-len", str(params.min_chrom_len),
                "--chrom-prefix", str(params.chrom_prefix),
            ] + list(params.group_args)
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_epi_drift_run at {current_time} '
                        f'completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_epi_drift_run: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_epi_drift_run: {e}\n")
        finally:
            rule_logger.info(f"fiberseq_epi_drift_run finished at "
                             f"{time.strftime('%Y%m%d_%H%M%S', time.localtime())}")


rule fiberseq_epi_drift_fig:
    """fig9-style drift plot (histogram + drift-vs-mean scatter)."""
    input:
        windows_tsv = outdir + "/epigenetic_drift_windows.tsv",
        summary_json = outdir + "/epigenetic_drift.json",
    output:
        plots = [outdir + "/epigenetic_drift." + fmt for fmt in _formats_fig],
    log:
        logdir + "/epi_drift_fig.log"
    threads: rule_res("fiberseq_epi_drift").get("threads", 2)
    benchmark: logdir + "/epi_drift_fig.benchmark.tsv"
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_epi_drift", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 4000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_epi_drift_fig")
    params:
        script = ROOT_DIR + "/modules/fiberseq/epi_drift/bin/epi_drift_plot.py",
        formats = _formats_fig,
        dpi = _p.get("fiberseq_epi_drift_fig", {}).get("dpi", 200),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_epi_drift_fig", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start epi_drift fig at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = receipt_script(f"fiberseq_epi_drift_fig_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "fiberseq_epi_drift_fig.sh")
            cmd = [
                "python", str(params.script),
                "--windows", str(input.windows_tsv),
                "--summary", str(input.summary_json),
                "--out-base", outdir + "/epigenetic_drift",
                "--dpi", str(params.dpi),
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_epi_drift_fig at {current_time} '
                        f'completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_epi_drift_fig: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_epi_drift_fig: {e}\n")
        finally:
            rule_logger.info(f"fiberseq_epi_drift_fig finished at "
                             f"{time.strftime('%Y%m%d_%H%M%S', time.localtime())}")
