
include: "../common/common.smk"


outdir=config.get("outdir","output"); logdir=config.get("logdir","logs")
plink2=config.get("Procedure",{}).get("plink2") or "plink2"
vcf=config.get("input_vcf") or outdir+"/variants/filtered/population_snps.vcf.gz"
prefix=outdir+"/analysis/plink/population_snps"

def run_plink(rule_name, log_path, output_path, commands):
    os.makedirs(os.path.dirname(log_path),exist_ok=True); os.makedirs(os.path.dirname(output_path),exist_ok=True)
    script = receipt_script(f"{rule_name}_{time.strftime('%Y%m%d_%H%M%S')}.sh")
    os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
    plain_script=os.path.join(os.path.dirname(output_path), f"{rule_name}.sh")
    with open(log_path, "a") as log_handle: log_handle.write("")
    with open(script,"w") as handle:
        handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        for command in commands: handle.write(" ".join(str(x) for x in command)+"\n")
    shutil.copyfile(script, plain_script)
    if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
        script_only_touch(output)
    else:
        shell(f"bash {plain_script} >> {log_path} 2>&1")
rule plink2_population_pgen:
    input: vcf=vcf, tbi=vcf+".tbi"
    output:
        # intermediate: consumed by plink2_population_bed/plink2_population_pca
        # (registered outputs are .eigenvec/.eigenval), Params.plink2_population.keep_intermediates
        pgen=maybe_temp(prefix+".pgen", "plink2_population"), pvar=maybe_temp(prefix+".pvar", "plink2_population"), psam=maybe_temp(prefix+".psam", "plink2_population")
    log: logdir+"/analysis/plink2_population_pgen.log"
    benchmark: logdir+"/analysis/plink2_population_pgen.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "plink2_population_pgen", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda: "plink2_population.yaml"
    container:
        sif("plink2_population.yaml", "plink2_population_pgen")
    params: plink2=plink2, prefix=prefix
    run:
        run_plink("plink2_population_pgen",str(log),output.pgen,[[params.plink2,"--vcf",input.vcf,"dosage=DS","--make-pgen","--out",params.prefix]])

rule plink2_population_bed:
    input: pgen=prefix+".pgen", pvar=prefix+".pvar", psam=prefix+".psam"
    output: bed=maybe_temp(prefix+".bed", "plink2_population"), bim=maybe_temp(prefix+".bim", "plink2_population"), fam=maybe_temp(prefix+".fam", "plink2_population")
    log: logdir+"/analysis/plink2_population_bed.log"
    benchmark: logdir+"/analysis/plink2_population_bed.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "plink2_population_bed", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda: "plink2_population.yaml"
    container:
        sif("plink2_population.yaml", "plink2_population_bed")
    params: plink2=plink2, prefix=prefix
    run:
        run_plink("plink2_population_bed",str(log),output.bed,[[params.plink2,"--pfile",params.prefix,"--make-bed","--out",params.prefix]])

if config.get("Params", {}).get("analysis", {}).get("pca",True):
    rule plink2_population_pca:
        input: pgen=[prefix+".pgen",prefix+".pvar",prefix+".psam"]
        output: eigenvec=prefix+".eigenvec", eigenval=prefix+".eigenval"
        log: logdir+"/analysis/plink2_population_pca.log"
        conda: "plink2_population.yaml"
        container:
            sif("plink2_population.yaml", "plink2_population_pca")
        params: plink2=plink2, prefix=prefix, components=config.get("Params", {}).get("analysis", {}).get("pca_components",10)
        run:
            run_plink("plink2_population_pca",str(log),output.eigenvec,[[params.plink2,"--pfile",params.prefix,"--pca",str(params.components),"approx","--out",params.prefix]])

if config.get("Params", {}).get("analysis", {}).get("gwas",False):
    phenotype=config.get("Params", {}).get("analysis", {}).get("phenotype")
    if not phenotype: raise ValueError("Params.analysis.gwas=true requires phenotype")
    rule plink2_population_gwas:
        input:
            bed=[prefix+".bed",prefix+".bim",prefix+".fam"], phenotype=phenotype,
            covariates=config.get("Params", {}).get("analysis", {}).get("covariates",[])
        output: done=outdir+"/analysis/gwas/gwas.done"
        log: logdir+"/analysis/plink2_population_gwas.log"
        conda: "plink2_population.yaml"
        container:
            sif("plink2_population.yaml", "plink2_population_gwas")
        params: plink2=plink2, prefix=prefix, pheno_name=config.get("Params", {}).get("analysis", {}).get("phenotype_name"), covar_names=config.get("Params", {}).get("analysis", {}).get("covariate_names"), gwas_prefix=outdir+"/analysis/gwas/population"
        run:
            os.makedirs(os.path.dirname(str(output.done)),exist_ok=True)
            cmd=[params.plink2,"--bfile",params.prefix,"--pheno",input.phenotype]
            if params.pheno_name: cmd += ["--pheno-name",params.pheno_name]
            if input.covariates: cmd += ["--covar",input.covariates[0]]
            if params.covar_names: cmd += ["--covar-name",params.covar_names]
            cmd += ["--glm","hide-covar","allow-no-covars","--out",params.gwas_prefix]
            script = receipt_script(f"plink2_population_gwas_{time.strftime('%Y%m%d_%H%M%S')}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script=os.path.join(os.path.dirname(str(output.done)), "plink2_population_gwas.sh")
            with open(str(log),"w") as h: h.write("")
            with open(script,"w") as h: h.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(str(x) for x in cmd)+"\ntouch "+str(output.done)+"\necho \"plink2 completed successfully\"\n")
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log} 2>&1")