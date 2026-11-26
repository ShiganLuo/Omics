include: "../../common/common.smk"
epifany = config.get("Procedure", {}).get("epifany") or "Epifany"

indir = config.get("indir", "data/psm_fdr")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
# 参数精确到 rule(2026-11 规范):config 全量路径字面链
# Params.protein_inference.protein_inference.<key>,无 module 级共享回退。

rule protein_inference:
    input:
        idxml = indir + "/{sample_id}/{sample_id}_filtered.idXML"
    output:
        protein_idxml = outdir + "/{sample_id}/{sample_id}_protein.idXML"
    log:
        logdir + "/{sample_id}/protein_inference.log"
    benchmark: logdir + "/{sample_id}/protein_inference.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "protein_inference", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml", "protein_inference")
    params:
        epifany = epifany,
        method = config.get("Params", {}).get("protein_inference", {}).get("protein_inference", {}).get("method", "epifany"),
        greedy_group_resolution = config.get("Params", {}).get("protein_inference", {}).get("protein_inference", {}).get("greedy_group_resolution", "none")
    run:
        log_path = str(log)
        open(log_path, 'a').close()
        rule_logger = setup_logger("protein_inference", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start protein inference for sample {wildcards.sample_id} at {current_time}")
            # Epifany segfaults on idXML without peptide hits (e.g. when no PSM
            # survived FDR); fail early with an actionable message instead.
            # 测试模式输入是空占位,跳过内容校验(工具不会真正执行)。
            if not SCRIPT_ONLY:
                n_hits = 0
                with open(input.idxml) as fh:
                    for line in fh:
                        if "<PeptideHit" in line:
                            n_hits += 1
                            break
                if n_hits == 0:
                    raise RuntimeError(
                        f"no peptide hits in {input.idxml}: zero PSMs survived FDR "
                        f"for sample {wildcards.sample_id}; check the psm_fdr step "
                        "(search/identification quality) before protein inference"
                    )
            script = receipt_script(f"protein_inference_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/protein_inference.sh")
            os.makedirs(os.path.dirname(plain_script), exist_ok=True)  # 目录靠 snakemake 建,回执复制前先建
            # Epifany (OpenMS 3.1) has no -algorithm:method: it only runs the
            # Bayesian inference; greedy_group_resolution is a TOPP-level
            # option (-greedy_group_resolution), not an algorithm subsection one.
            if params.method not in ("", "epifany"):
                raise ValueError(
                    f"Params.protein_inference.protein_inference.method={params.method!r} is not supported: "
                    "this Epifany build only provides Bayesian 'epifany' inference"
                )
            cmd = [
                params.epifany,
                "-in", input.idxml,
                "-out", output.protein_idxml,
                "-greedy_group_resolution", params.greedy_group_resolution,
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f"echo 'Protein inference for {wildcards.sample_id} was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule protein_inference was call failed,error: {e}")
            raise RuntimeError(f"rule protein_inference was call failed,error: {e}")

rule protein_inference_result:
    input:
        protein_idxml = outdir + "/{sample_id}/{sample_id}_protein.idXML"
