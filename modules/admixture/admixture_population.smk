
include: "../common/common.smk"
import os,time
outdir=config.get("outdir","output"); logdir=config.get("logdir","logs"); admixture=config.get("Procedure",{}).get("admixture") or "admixture"; prefix=outdir+"/analysis/plink/population_snps"; ks=config.get("Params",{}).get("analysis",{}).get("admixture_k",[2,3,4,5,6])
rule admixture_population:
 input: bed=prefix+".bed",bim=prefix+".bim",fam=prefix+".fam"
 output: q=outdir+"/analysis/structure/admixture.K{k}.Q",p=outdir+"/analysis/structure/admixture.K{k}.P"
 log: logdir+"/analysis/admixture.K{k}.log"
 conda: "admixture_population.yaml"
 container:
     sif("admixture_population.yaml")
 params: admixture=admixture,bed=prefix+".bed",seed=config.get("Params",{}).get("analysis",{}).get("admixture_seed",2026)
 threads: 4
 run:
  os.makedirs(os.path.dirname(str(log)),exist_ok=True);os.makedirs(os.path.dirname(str(output.q)),exist_ok=True)
  script=os.path.join(os.path.dirname(str(output.q)),f"admixture_population_{wildcards.k}_{time.strftime('%Y%m%d_%H%M%S')}.sh")
  cmd=[params.admixture,"--seed",str(params.seed),f"-j{threads}",params.bed,str(wildcards.k)]
  with open(str(log),"w") as h:h.write("")
  with open(script,"w") as h:h.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(str(x) for x in cmd)+"\nmv "+params.bed+"."+str(wildcards.k)+".Q "+str(output.q)+"\nmv "+params.bed+"."+str(wildcards.k)+".P "+str(output.p)+"\n")
  shell(f"bash {script} >> {log} 2>&1")
