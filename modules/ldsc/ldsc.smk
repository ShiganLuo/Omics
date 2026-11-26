include: "../common/common.smk"


outdir = config.get("outdir", "results/ldsc")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})

gwas = config.get("Params", {}).get("ldsc", {}).get("gwas", {}) or {}
traits = sorted(gwas.keys())
BIN = ROOT_DIR + "/modules/ldsc/bin"

# 参数精确到 rule:annot_prefix 按 rule 各自配置
# (Params.ldsc.ldsc_annot.annot_prefix / Params.ldsc.ldsc_h2_cts.annot_prefix;
# 无 module 级共享回退)。
_annot_prefix_ldsc_annot = config.get("Params", {}).get("ldsc", {}).get("ldsc_annot", {}).get("annot_prefix", "peaks") or "peaks"
_annot_prefix_ldsc_h2_cts = config.get("Params", {}).get("ldsc", {}).get("ldsc_h2_cts", {}).get("annot_prefix", "peaks") or "peaks"

wildcard_constraints:
    trait = ".+",

rule ldsc_annot:
    """Split FIRE peak clusters into LDSC annotation BEDs."""
    input:
        peaks = lambda wildcards: _require(config.get("Params", {}).get("ldsc", {}).get("annot_bed"),
                                           "Params.ldsc.annot_bed", "LDSC"),
    output:
        manifest = outdir + "/annot/" + _annot_prefix_ldsc_annot + ".annot_manifest.tsv",
    log:
        logdir + "/ldsc_annot.log"
    threads: rule_res("ldsc_annot").get("threads", 1)
    benchmark: logdir + "/ldsc_annot.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "ldsc_annot", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "ldsc.yaml"
    container:
        sif("ldsc.yaml", "ldsc_annot")
    params:
        script = BIN + "/ldsc_annot.py",
        prefix = _annot_prefix_ldsc_annot,
        annot_cluster_from_name = config.get("Params", {}).get("ldsc", {}).get("annot_cluster_from_name"),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("ldsc_annot", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ldsc_annot at {current_time}")
            group_outdir = outdir + "/annot"
            os.makedirs(group_outdir, exist_ok=True)
            script = receipt_script(f"ldsc_annot_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(group_outdir, f"ldsc_annot.sh")
            cmd = [
                "python", params.script,
                "-p", str(input.peaks),
                "-o", group_outdir,
                "--prefix", params.prefix,
                "-l", log_path,
            ]
            if params.annot_cluster_from_name:
                cmd.append("-c")
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "ldsc_annot at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during ldsc_annot: {e}\n")
            raise RuntimeError(f"Error occurred during ldsc_annot: {e}\n")

rule ldsc_munge:
    """Munge one GWAS sumstats file with munge_sumstats.py."""
    input:
        sumstats = lambda wildcards: _require(
            gwas.get(wildcards.trait), f"Params.ldsc.gwas.{wildcards.trait}", "LDSC"),
    output:
        munged = outdir + "/{trait}/{trait}.sumstats.gz",
    log:
        logdir + "/ldsc_munge_{trait}.log"
    threads: rule_res("ldsc_munge").get("threads", 1)
    benchmark: logdir + "/ldsc_munge_{trait}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "ldsc_munge", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "ldsc.yaml"
    container:
        sif("ldsc.yaml", "ldsc_munge")
    params:
        munge = config.get("Procedure", {}).get("munge_sumstats") or "munge_sumstats.py",
        extra = config.get("Params", {}).get("ldsc", {}).get("munge_args", ""),
        merge_alleles = config.get("Params", {}).get("ldsc", {}).get("merge_alleles"),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("ldsc_munge", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ldsc_munge for {wildcards.trait} at {current_time}")
            trait_outdir = outdir + f"/{wildcards.trait}"
            os.makedirs(trait_outdir, exist_ok=True)
            script = receipt_script(f"ldsc_munge_{wildcards.trait}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(trait_outdir, f"ldsc_munge.sh")
            cmd = [
                params.munge,
                "--sumstats", str(input.sumstats),
                "--out", trait_outdir + f"/{wildcards.trait}",
                "--merge-alleles", _require(params.merge_alleles,
                                            "Params.ldsc.merge_alleles", "LDSC"),
            ]
            if params.extra:
                cmd += str(params.extra).split()
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "ldsc_munge for {wildcards.trait} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during ldsc_munge for {wildcards.trait}: {e}\n")
            raise RuntimeError(f"Error occurred during ldsc_munge for {wildcards.trait}: {e}\n")

rule ldsc_h2_cts:
    """Partitioned heritability (ldsc.py --h2-cts) for one GWAS trait."""
    input:
        munged = outdir + "/{trait}/{trait}.sumstats.gz",
        manifest = outdir + "/annot/" + _annot_prefix_ldsc_h2_cts + ".annot_manifest.tsv",
    output:
        results = outdir + "/{trait}/{trait}.cts.results",
    log:
        logdir + "/ldsc_h2_cts_{trait}.log"
    threads: rule_res("ldsc_h2_cts").get("threads", 1)
    benchmark: logdir + "/ldsc_h2_cts_{trait}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "ldsc_h2_cts", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "ldsc.yaml"
    container:
        sif("ldsc.yaml", "ldsc_h2_cts")
    params:
        ldsc = config.get("Procedure", {}).get("ldsc") or "ldsc.py",
        ref_ld = config.get("Params", {}).get("ldsc", {}).get("ref_ld_chr"),
        w_ld = config.get("Params", {}).get("ldsc", {}).get("w_ld_chr"),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("ldsc_h2_cts", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start ldsc_h2_cts for {wildcards.trait} at {current_time}")
            trait_outdir = outdir + f"/{wildcards.trait}"
            os.makedirs(trait_outdir, exist_ok=True)
            script = receipt_script(f"ldsc_h2_cts_{wildcards.trait}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(trait_outdir, f"ldsc_h2_cts.sh")
            cmd = [
                params.ldsc,
                "--h2-cts", str(input.munged),
                "--ref-ld-chr", str(_require(params.ref_ld, "Params.ldsc.ref_ld_chr", "LDSC")),
                "--w-ld-chr", str(_require(params.w_ld, "Params.ldsc.w_ld_chr", "LDSC")),
                "--out", trait_outdir + f"/{wildcards.trait}",
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "ldsc_h2_cts for {wildcards.trait} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during ldsc_h2_cts for {wildcards.trait}: {e}\n")
            raise RuntimeError(f"Error occurred during ldsc_h2_cts for {wildcards.trait}: {e}\n")

rule ldsc_result:
    """Panel of per-trait LDSC results."""
    input:
        results = expand(outdir + "/{trait}/{trait}.cts.results", trait=traits),
