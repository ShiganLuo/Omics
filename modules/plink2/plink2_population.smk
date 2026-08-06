
include: "../common/common.smk"
import os, time

outdir=config.get("outdir","output"); logdir=config.get("logdir","logs")
plink2=config.get("Procedure",{}).get("plink2") or "plink2"
vcf=config.get("input_vcf") or outdir+"/variants/filtered/population_snps.vcf.gz"
prefix=outdir+"/analysis/plink/population_snps"
analysis=config.get("Params",{}).get("analysis",{})

def run_plink(rule_name, log_path, output_path, commands):
    os.makedirs(os.path.dirname(log_path),exist_ok=True); os.makedirs(os.path.dirname(output_path),exist_ok=True)
    script=os.path.join(os.path.dirname(output_path),f"{rule_name}_{time.strftime('%Y%m%d_%H%M%S')}.sh")
    with open(log_path,"w") as log_handle: log_handle.write("")
    with open(script,"w") as handle:
        handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        for command in commands: handle.write(" ".join(str(x) for x in command)+"\n")
    shell(f"bash {script} >> {log_path} 2>&1")

rule plink2_population_pgen:
    input: vcf=vcf, tbi=vcf+".tbi"
    output:
        pgen=prefix+".pgen", pvar=prefix+".pvar", psam=prefix+".psam"
    log: logdir+"/analysis/plink2_population_pgen.log"
    conda: "plink2_population.yaml"
    container:
        sif("plink2_population.yaml")
    params: plink2=plink2, prefix=prefix
    run:
        run_plink("plink2_population_pgen",str(log),output.pgen,[[params.plink2,"--vcf",input.vcf,"dosage=DS","--make-pgen","--out",params.prefix]])

rule plink2_population_bed:
    input: pgen=prefix+".pgen", pvar=prefix+".pvar", psam=prefix+".psam"
    output: bed=prefix+".bed", bim=prefix+".bim", fam=prefix+".fam"
    log: logdir+"/analysis/plink2_population_bed.log"
    conda: "plink2_population.yaml"
    container:
        sif("plink2_population.yaml")
    params: plink2=plink2, prefix=prefix
    run:
        run_plink("plink2_population_bed",str(log),output.bed,[[params.plink2,"--pfile",params.prefix,"--make-bed","--out",params.prefix]])

if analysis.get("pca",True):
    rule plink2_population_pca:
        input: pgen=[prefix+".pgen",prefix+".pvar",prefix+".psam"]
        output: eigenvec=prefix+".eigenvec", eigenval=prefix+".eigenval"
        log: logdir+"/analysis/plink2_population_pca.log"
        conda: "plink2_population.yaml"
        container:
            sif("plink2_population.yaml")
        params: plink2=plink2, prefix=prefix, components=analysis.get("pca_components",10)
        run:
            run_plink("plink2_population_pca",str(log),output.eigenvec,[[params.plink2,"--pfile",params.prefix,"--pca",str(params.components),"approx","--out",params.prefix]])

if analysis.get("gwas",False):
    phenotype=analysis.get("phenotype")
    if not phenotype: raise ValueError("Params.analysis.gwas=true requires phenotype")
    rule plink2_population_gwas:
        input:
            bed=[prefix+".bed",prefix+".bim",prefix+".fam"], phenotype=phenotype,
            covariates=analysis.get("covariates",[])
        output: done=outdir+"/analysis/gwas/gwas.done"
        log: logdir+"/analysis/plink2_population_gwas.log"
        conda: "plink2_population.yaml"
        container:
            sif("plink2_population.yaml")
        params: plink2=plink2, prefix=prefix, pheno_name=analysis.get("phenotype_name"), covar_names=analysis.get("covariate_names"), gwas_prefix=outdir+"/analysis/gwas/population"
        run:
            os.makedirs(os.path.dirname(str(output.done)),exist_ok=True)
            cmd=[params.plink2,"--bfile",params.prefix,"--pheno",input.phenotype]
            if params.pheno_name: cmd += ["--pheno-name",params.pheno_name]
            if input.covariates: cmd += ["--covar",input.covariates[0]]
            if params.covar_names: cmd += ["--covar-name",params.covar_names]
            cmd += ["--glm","hide-covar","allow-no-covars","--out",params.gwas_prefix]
            script=os.path.join(os.path.dirname(str(output.done)),f"plink2_population_gwas_{time.strftime('%Y%m%d_%H%M%S')}.sh")
            with open(str(log),"w") as h: h.write("")
            with open(script,"w") as h: h.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(str(x) for x in cmd)+"\ntouch "+str(output.done)+"\n")
            shell(f"bash {script} >> {log} 2>&1")
