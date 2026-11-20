include: "../../common/common.smk"


indir = config.get("indir", "common/4_elements")
outdir = config.get("outdir", "fiberseq_analysis/track")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})

rule fiberseq_track_run:
    """bedGraph/bigWig accessibility and m6A density tracks + trackDb."""
    input:
        msp = indir + "/{sample_id}/{sample_id}.msp.bed.gz",
        m6a = indir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        nuc = indir + "/{sample_id}/{sample_id}.nuc.bed.gz",
        chrom_sizes = config.get("genome",{}).get("fasta", "") + ".fai",
    output:
        acc_bw = outdir + "/{sample_id}/{sample_id}.accessibility.bw",
        m6a_bw = outdir + "/{sample_id}/{sample_id}.m6a_density.bw",
        trackdb = outdir + "/{sample_id}/trackDb.txt",
    log:
        logdir + "/{sample_id}/fiberseq_track_run.log"
    threads: rule_res("fiberseq_track").get("threads", 4)
    benchmark: logdir + "/{sample_id}/fiberseq_track_run.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "fiberseq_track", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 65000),
    conda:
        "../fiberseq.yaml"
    container:
        sif("../fiberseq.yaml", "fiberseq_track_run")
    params:
        script = ROOT_DIR + "/modules/fiberseq/bin/fiberseq_track.py",
        # 参数精确到 rule:Params.fiberseq_track.fiberseq_track_run.<key>(config 全量路径读取)
        bin = config.get("Params", {}).get("fiberseq_track", {}).get("fiberseq_track_run", {}).get("bin", 10),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("fiberseq_track_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start track generation for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.acc_bw))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"fiberseq_track_run_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            # timestamped copy stays in log/cmd as history; the outdir copy
            # (no timestamp, overwritten each run) documents the current command
            # next to the results it produced
            plain_script = os.path.join(sample_outdir, f"fiberseq_track_run_{wildcards.sample_id}.sh")
            cmd = [
                "python", params.script,
                "-s", str(input.msp),
                "-a", str(input.m6a),
                "-u", str(input.nuc),
                "-c", str(input.chrom_sizes),
                "-o", sample_outdir,
                "--sample", wildcards.sample_id,
                "-b", str(params.bin),
                "-l", log_path,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "fiberseq_track_run for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during fiberseq_track_run for {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during fiberseq_track_run for {wildcards.sample_id}: {e}\n")

rule fiberseq_track_result:
    input:
        tracks = expand(outdir + "/{sample_id}/{sample_id}.accessibility.bw",
                        sample_id=samples),
