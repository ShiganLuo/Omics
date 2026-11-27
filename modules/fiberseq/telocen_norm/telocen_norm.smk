include: "../../common/common.smk"


outdir = config.get("outdir", "fiberseq_analysis/telocen_norm")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
_inputs = config.get("inputs", {}).get("telocen_norm", {})
_bams = _inputs.get("bam") or {}
_p = config.get("Params", {}).get("fiberseq_telocen_norm", {})
_groups = _p.get("groups") or {}
_glabels = list(_groups)
_formats_fig = _p.get("fiberseq_telocen_norm_fig", {}).get("formats") or ["png"]
if isinstance(_formats_fig, str):
    _formats_fig = [_formats_fig]


def _group_bam_list(wildcards):
    """BAM paths of all samples in the two comparison groups (deduped)."""
    out = []
    for label in _glabels:
        for s in _groups.get(label) or []:
            bam = _bams.get(s)
            if bam and bam not in out:
                out.append(bam)
    return out


def _group_cli():
    """['--group', 'youth=y', ...] — both rules need the group mapping."""
    args = []
    for label in _glabels:
        args += ["--group", f"{label}={','.join(_groups.get(label) or [])}"]
    return args


def _bam_cli():
    """['--bam', 'S=path', ...] — only the compute rule scans BAMs."""
    args = []
    for label in _glabels:
        for s in _groups.get(label) or []:
            if s in _bams:
                args += ["--bam", f"{s}={_bams[s]}"]
    return args


rule fiberseq_telocen_norm_run:
    """Depth-normalized telomere/centromere region metrics for two groups
    (fig2 data): region_normalized.tsv + class_diff_normalized.tsv."""
    input:
        regions = lambda wildcards: _require(genome_cfg.get("telocen_bed"),
                                             "genome.telocen_bed",
                                             "fiberseq_telocen_norm"),
        bams = _group_bam_list,
    output:
        region_tsv = outdir + "/region_normalized.tsv",
        class_diff = outdir + "/class_diff_normalized.tsv",
    log:
        logdir + "/telocen_norm_run.log"
    threads: rule_res("fiberseq_telocen_norm").get("threads", 16)
    benchmark: logdir + "/telocen_norm_run.benchmark.tsv"
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_telocen_norm", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_telocen_norm_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/telocen_norm/bin/telocen_normalized.py",
        group_args = lambda wildcards: _group_cli() + _bam_cli(),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_telocen_norm_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start telocen_norm run at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            if not _glabels or len(_glabels) != 2:
                raise ValueError(
                    "Params.fiberseq_telocen_norm.groups must define exactly "
                    f"2 groups (label -> [samples]), got: {_groups!r}")
            for label in _glabels:
                missing = [s for s in (_groups.get(label) or [])
                           if s not in _bams]
                if missing:
                    raise ValueError(
                        f"group {label!r} references samples without a BAM "
                        f"entry in inputs.telocen_norm.bam: {missing}")
            script = receipt_script(f"fiberseq_telocen_norm_run_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "fiberseq_telocen_norm_run.sh")
            cmd = [
                "python", str(params.script),
                "--regions", str(input.regions),
                "--outdir", outdir,
                "--threads", str(threads),
            ] + list(params.group_args)
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_telocen_norm_run at {current_time} '
                        f'completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_telocen_norm_run: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_telocen_norm_run: {e}\n")
        finally:
            rule_logger.info(f"fiberseq_telocen_norm_run finished at "
                             f"{time.strftime('%Y%m%d_%H%M%S', time.localtime())}")


rule fiberseq_telocen_norm_fig:
    """fig2-style paired plot over region classes (depth-normalized)."""
    input:
        region_tsv = outdir + "/region_normalized.tsv",
    output:
        plots = [outdir + "/telocen_normalized." + fmt for fmt in _formats_fig],
    log:
        logdir + "/telocen_norm_fig.log"
    threads: rule_res("fiberseq_telocen_norm").get("threads", 2)
    benchmark: logdir + "/telocen_norm_fig.benchmark.tsv"
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_telocen_norm", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 4000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_telocen_norm_fig")
    params:
        script = ROOT_DIR + "/modules/fiberseq/telocen_norm/bin/plot_telocen_norm.py",
        group_args = lambda wildcards: _group_cli(),
        formats = _formats_fig,
        dpi = _p.get("fiberseq_telocen_norm_fig", {}).get("dpi", 200),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_telocen_norm_fig", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start telocen_norm fig at {current_time}")
            os.makedirs(outdir, exist_ok=True)
            script = receipt_script(f"fiberseq_telocen_norm_fig_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(outdir, "fiberseq_telocen_norm_fig.sh")
            cmd = [
                "python", str(params.script),
                "--input", str(input.region_tsv),
                "--out-base", outdir + "/telocen_normalized",
                "--dpi", str(params.dpi),
            ] + list(params.group_args)
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_telocen_norm_fig at {current_time} '
                        f'completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_telocen_norm_fig: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_telocen_norm_fig: {e}\n")
        finally:
            rule_logger.info(f"fiberseq_telocen_norm_fig finished at "
                             f"{time.strftime('%Y%m%d_%H%M%S', time.localtime())}")
