include: "../../common/common.smk"
percolator = config.get("Procedure", {}).get("percolator") or "PercolatorAdapter"

shell.prefix("set -x; set -e;")


indir = config.get("indir", "data/search_engine")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
# 参数精确到 rule(2026-11 规范):config 全量路径字面链,无 module 级共享回退:
#   Params.psm_rescoring.search_engine                    # 引擎选择,每 rule 独立键
#   Params.psm_rescoring.percolator.psm_rescoring.<key>   # 参数组单消费,挂 rule 名下
#   Params.psm_rescoring.psm_feature_extractor.psm_rescoring.<key>

def get_input_for_psm_rescoring(wildcards):
    """Get input files for PSM rescoring."""
    search_engine = config.get("Params", {}).get("psm_rescoring", {}).get("search_engine", None)
    if search_engine == "comet":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_comet.idXML"}
    elif search_engine == "msgf":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_msgf.idXML"}
    elif search_engine == "sage":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_sage.idXML"}
    elif search_engine == "msfragger":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_msfragger.idXML"}
    else:
        raise ValueError("No search engine specified in the configuration.")

rule psm_rescoring:
    input:
        unpack(get_input_for_psm_rescoring)
    output:
        scored_idxml = outdir + "/{sample_id}/{sample_id}_scored.idXML"
    log:
        logdir + "/{sample_id}/psm_rescoring.log"
    benchmark: logdir + "/{sample_id}/psm_rescoring.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "psm_rescoring", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 32000),
    conda:
        "../openms.yaml"
    container:
        sif("../openms.yaml", "psm_rescoring")
    params:
        percolator = percolator,
        feature_extractor = config.get("Procedure", {}).get("psm_feature_extractor") or "PSMFeatureExtractor",
        extra = config.get("Params", {}).get("psm_rescoring", {}).get("psm_feature_extractor", {}).get("psm_rescoring", {}).get("extra", ""),
        train_FDR = config.get("Params", {}).get("psm_rescoring", {}).get("percolator", {}).get("psm_rescoring", {}).get("train_FDR", 0.05),
        test_FDR = config.get("Params", {}).get("psm_rescoring", {}).get("percolator", {}).get("psm_rescoring", {}).get("test_FDR", 0.05),
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("psm_rescoring", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start PSM rescoring for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"psm_rescoring_{wildcards.sample_id}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(outdir, f"{wildcards.sample_id}/psm_rescoring.sh")
            feat_idxml = os.path.join(outdir, f"{wildcards.sample_id}/{wildcards.sample_id}_features.idXML")
            # PercolatorAdapter needs search-engine-specific features, which
            # PSMFeatureExtractor adds between the search and the rescoring.
            feat_cmd = [
                params.feature_extractor,
                "-in", input.idxml,
                "-out", feat_idxml,
            ]
            if params.extra:
                feat_cmd += str(params.extra).split()
            # OpenMS 3.x PercolatorAdapter uses -trainFDR/-testFDR (no
            # underscores); the legacy 2.x -feature option no longer exists.
            perc_cmd = [
                params.percolator,
                "-in", feat_idxml,
                "-out", output.scored_idxml,
                "-trainFDR", str(params.train_FDR),
                "-testFDR", str(params.test_FDR),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(feat_cmd) + "\n")
                f.write(" ".join(perc_cmd) + "\n")
                f.write(f"echo 'PSM rescoring for {wildcards.sample_id} was completed successfully'\n")
            shutil.copyfile(script, plain_script)
            if not SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule psm_rescoring was call failed,error: {e}")
            raise RuntimeError(f"rule psm_rescoring was call failed,error: {e}")

rule psm_rescoring_result:
    input:
        expand(outdir + "/{sid}/{sid}_scored.idXML", sid=samples)
