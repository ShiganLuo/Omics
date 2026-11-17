include: "../../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
# gtf resolved per-genome in rules
# 参数精确到 rule:Params.tailer.<rule_name>.<key>(rule_param 见 common.smk;
# 无 module 级共享回退,tailer_global 独立键;rev_comp/threshold 与
# star_3pg_gene_specific 共享,各 rule 各自独立键)。
_tailer_params = config.get("Params", {}).get("tailer", {})
_read_num_tailer_global = rule_param(_tailer_params, "tailer_global", "read", 2)
_threshold_tailer_global = rule_param(_tailer_params, "tailer_global", "threshold", 100)
_rev_comp_tailer_global = rule_param(_tailer_params, "tailer_global", "rev_comp", False)

rule tailer_global:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam",
        bai = indir + "/{sample_id}/{sample_id}.bam.bai",
    output:
        tail = outdir + "/{genome}/{sample_id}/{sample_id}_tail.csv",
    log:
        logdir + "/{genome}/{sample_id}/tailer.log"
    threads: rule_res("tailer_global").get("threads", 1)
    benchmark: logdir + "/{genome}/{sample_id}/tailer.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "tailer_global", input_size=input_size_named(input, 'bam'),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "../tailer.yaml"
    container:
        sif("../tailer.yaml", "tailer_global")
    params:
        gtf = lambda wildcards: config.get("genome", {}).get(wildcards.genome, {}).get("gtf"),
        read_num = _read_num_tailer_global,
        threshold = _threshold_tailer_global,
        rev_comp = "--rev_comp" if _rev_comp_tailer_global else "",
        Tailer = config.get("Procedure", {}).get("tailer") or config.get("Procedure", {}).get("Tailer") or "Tailer"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("tailer_global", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start tailer_global for sample {wildcards.sample_id} at {current_time}")

            outdir_tail = os.path.dirname(str(output.tail))
            os.makedirs(outdir_tail, exist_ok=True)
            generated_tail = os.path.join(os.path.dirname(str(input.bam)), f"{wildcards.sample_id}_tail.csv")
            if os.path.exists(generated_tail):
                os.remove(generated_tail)

            cmd = f"{params.Tailer} -a {params.gtf} -read {params.read_num} -t {params.threshold} {params.rev_comp} {input.bam} >> {log_path} 2>&1"
            rule_logger.info(f"Running: {cmd}")
            shell(cmd)

            if not os.path.exists(generated_tail):
                raise FileNotFoundError(
                    f"Tailer finished but did not create expected file: {generated_tail}"
                )
            shutil.copy2(generated_tail, output.tail)

            rule_logger.info(f"tailer_global completed for sample {wildcards.sample_id}")
        except Exception as e:
            logger.error(f"Error occurred during tailer_global for sample {wildcards.sample_id}: {e}")
            raise e