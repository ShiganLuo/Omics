include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/motif")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
_motif_params = config.get("Params", {}).get("fiberseq_motif", {})
# 参数精确到 rule:formats 按 rule 独立配置(Params.fiberseq_motif.<rule_name>.formats,
# 无 module 级共享回退;rule_param 见 common.smk),每 rule 一个 parse 期变量。


def _rule_formats(rule_name):
    """Per-rule figure formats (parse 期解析,output 展开与 CLI 转发共用)。"""
    fmts = rule_param(_motif_params, rule_name, "formats", ["pdf"])
    return [fmts] if isinstance(fmts, str) else (fmts or ["pdf"])


_formats_profile = _rule_formats("fiberseq_motif_profile")
_formats_result = _rule_formats("fiberseq_motif_result")

rule fiberseq_motif_scan:
    """Scan a MEME motif file in FIRE peak sequences (any species)."""
    input:
        peaks = indir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        fasta = lambda wildcards: _require(genome_cfg.get("fasta"), "genome.fasta", "fiberseq_motif"),
        motif_file = lambda wildcards: _require(genome_cfg.get("motif_file"), "genome.motif_file", "fiberseq_motif"),
    output:
        sites = outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
    log:
        logdir + "/{sample_id}/fiberseq_motif_scan.log"
    threads: rule_res("fiberseq_motif").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_motif_scan.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_motif", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 4000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_motif_scan")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/motif_scan.py",
        pwm_threshold = rule_param(_motif_params, "fiberseq_motif_scan", "pwm_threshold", 0.8),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_motif_scan", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start motif scan for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.sites))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_motif_scan_{wildcards.sample_id}_{current_time}.sh")
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_motif_scan_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.motif_file),
                "-r", str(input.peaks),
                "-f", str(input.fasta),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-t", str(params.pwm_threshold),
                "--threads", str(threads),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_motif_scan for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_motif_scan for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_motif_scan for {wildcards.sample_id}: {e}\n")

rule fiberseq_motif_profile:
    """Aggregate m6A/accessibility/nucleosome profiles around motif sites."""
    input:
        sites = outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
    output:
        profile = outdir + "/{sample_id}/{sample_id}.motif_profile.tsv",
        scores = outdir + "/{sample_id}/{sample_id}.motif_footprint_scores.tsv",
        # ② 按格式展开：每种格式一个输出文件（base 不带扩展名）
        plots = [outdir + "/{sample_id}/{sample_id}.motif_profile." + fmt
                 for fmt in _formats_profile],
    log:
        logdir + "/{sample_id}/fiberseq_motif_profile.log"
    threads: rule_res("fiberseq_motif").get("threads", 1)
    benchmark: logdir + "/{sample_id}/fiberseq_motif_profile.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_motif", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_motif_profile")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/motif_profile.py",
        flank = rule_param(_motif_params, "fiberseq_motif_profile", "flank", 250),
        min_fibers = rule_param(_motif_params, "fiberseq_motif_profile", "min_fibers", 10),
        formats = _formats_profile,   # ③ 转发
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("fiberseq_motif_profile", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start motif profile for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.profile))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_motif_profile_{wildcards.sample_id}_{current_time}.sh")
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_motif_profile_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-m", str(input.sites),
                "-a", str(input.m6a),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-f", str(params.flank),
                "-i", str(params.min_fibers),
                "-l", log_path,
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_motif_profile for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_motif_profile for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_motif_profile for {wildcards.sample_id}: {e}\n")

rule fiberseq_motif_result:
    input:
        sites = expand(outdir + "/{sample_id}/{sample_id}.motif_sites.bed",
                       sample_id=samples),
        profile = expand(outdir + "/{sample_id}/{sample_id}.motif_profile.tsv",
                         sample_id=samples),
        plots = expand(outdir + "/{sample_id}/{sample_id}.motif_profile.{fmt}",
                       sample_id=samples, fmt=_formats_result),

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_motif_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
