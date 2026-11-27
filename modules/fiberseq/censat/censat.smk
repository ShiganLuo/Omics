include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/censat")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
_formats_run = config.get("Params", {}).get("fiberseq_censat", {}).get("fiberseq_censat_run", {}).get("formats") or ["pdf"]
if isinstance(_formats_run, str):
    _formats_run = [_formats_run]
_formats_result = config.get("Params", {}).get("fiberseq_censat", {}).get("fiberseq_censat_result", {}).get("formats") or ["pdf"]
if isinstance(_formats_result, str):
    _formats_result = [_formats_result]

rule fiberseq_censat_run:
    """Region-class accessibility and methylation profiles (CDR/HOR/flank)."""
    input:
        regions = lambda wildcards: _require(genome_cfg.get("censat_bed"), "genome.censat_bed", "fiberseq_censat"),
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
        cpg = indir + "/{sample_id}/{sample_id}.cpg.bed.gz",
    output:
        regions_tsv = outdir + "/{sample_id}/{sample_id}.censat_regions.tsv",
        class_summary = outdir + "/{sample_id}/{sample_id}.censat_class_summary.tsv",
        profile = outdir + "/{sample_id}/{sample_id}.censat_profile.tsv",
        plots = [outdir + "/{sample_id}/{sample_id}.censat_profile." + fmt
                 for fmt in _formats_run],
    log:
        logdir + "/{sample_id}/fiberseq_censat_run.log"
    threads: rule_res("fiberseq_censat").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_censat_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_censat", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_censat_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/censat_profile.py",
        pileup = genome_cfg.get("mcpg_pileup"),
        msp_min_len = config.get("Params", {}).get("fiberseq_censat", {}).get("fiberseq_censat_run", {}).get("msp_min_len", 50),
        flank = config.get("Params", {}).get("fiberseq_censat", {}).get("fiberseq_censat_run", {}).get("flank", 10000),
        bin = config.get("Params", {}).get("fiberseq_censat", {}).get("fiberseq_censat_run", {}).get("bin", 500),
        formats = _formats_run,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_censat_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start censat profile for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.regions_tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"fiberseq_censat_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)
            plain_script = os.path.join(sample_outdir, f"fiberseq_censat_run_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-r", str(input.regions),
                "-s", str(input.msp),
                "-u", str(input.nuc),
                "-c", str(input.cpg),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-m", str(params.msp_min_len),
                "-f", str(params.flank),
                "-b", str(params.bin),
                "-l", log_path,
            ]
            if params.pileup:
                cmd += ["-p", str(params.pileup)]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_censat_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_censat_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_censat_run for {wildcards.sample_id}: {e}\n")

rule fiberseq_censat_result:
    input:
        class_summary = expand(outdir + "/{sample_id}/{sample_id}.censat_class_summary.tsv",
                               sample_id=samples),
        plots = expand(outdir + "/{sample_id}/{sample_id}.censat_profile.{fmt}",
                       sample_id=samples, fmt=_formats_result),
