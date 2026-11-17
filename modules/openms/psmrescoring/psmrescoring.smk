include: "../../common/common.smk"

shell.prefix("set -x; set -e;")


indir = config.get("indir", "data/search_engine")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

# Get parameters
# 参数精确到 rule(2026-11 规范;rule_param 见 common.smk),无 module 级共享回退:
#   Params.psm_rescoring.search_engine                    # 引擎选择,每 rule 独立键
#   Params.psm_rescoring.percolator.psm_rescoring.<key>   # 参数组单消费,挂 rule 名下
#   Params.psm_rescoring.psm_feature_extractor.psm_rescoring.<key>
percolator_params = config.get("Params", {}).get("psm_rescoring", {}).get("percolator", {})
_psm_feature_params = config.get("Params", {}).get("psm_rescoring", {}).get("psm_feature_extractor", {})
_params_root = config.get("Params", {})

# search engine selector:原 module 级共享读取 Params.search_engine →
# rule 级键 Params.psm_rescoring.search_engine(parse 期变量)。
_search_engine_psm_rescoring = rule_param(_params_root, "psm_rescoring", "search_engine", None)
# Get executables
percolator = config.get("Procedure", {}).get("percolator") or "PercolatorAdapter"
def get_input_for_psm_rescoring(wildcards):
    """Get input files for PSM rescoring."""
    if _search_engine_psm_rescoring == "comet":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_comet.idXML"}
    elif _search_engine_psm_rescoring == "msgf":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_msgf.idXML"}
    elif _search_engine_psm_rescoring == "sage":
        return {"idxml": indir + f"/{wildcards.sample_id}/{wildcards.sample_id}_sage.idXML"}
    elif _search_engine_psm_rescoring == "msfragger":
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
        extra = rule_param(_psm_feature_params, "psm_rescoring", "extra", ""),
        train_FDR = rule_param(percolator_params, "psm_rescoring", "train_FDR", 0.05),
        test_FDR = rule_param(percolator_params, "psm_rescoring", "test_FDR", 0.05),
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("psm_rescoring", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start PSM rescoring for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(RECEIPT_DIR, f"psm_rescoring_{wildcards.sample_id}_{current_time}.sh")
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
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"rule psm_rescoring was call failed,error: {e}")
            raise RuntimeError(f"rule psm_rescoring was call failed,error: {e}")

rule psm_rescoring_result:
    input:
        expand(outdir + "/{sid}/{sid}_scored.idXML", sid=samples)

    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "psm_rescoring_result", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
