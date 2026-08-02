
include: "../common/common.smk"
import os,time
from itertools import combinations
outdir=config.get("outdir","output"); logdir=config.get("logdir","logs")
vcftools=config.get("Procedure",{}).get("vcftools") or "vcftools"; vcf=config.get("input_vcf") or outdir+"/variants/filtered/population_snps.vcf.gz"
populations=config.get("populations",{}); pops=sorted(populations); analysis=config.get("Params",{}).get("analysis",{})
if analysis.get("diversity",True):
 rule vcftools_population_diversity:
  input: vcf=vcf,tbi=vcf+".tbi",popfile=lambda wc:outdir+"/metadata/populations/"+wc.population+".txt"
  output: pi=outdir+"/analysis/diversity/{population}.windowed.pi",tajima=outdir+"/analysis/diversity/{population}.Tajima.D"
  log: logdir+"/analysis/vcftools_diversity_{population}.log"
  conda: "vcftools_population.yaml"
  params: vcftools=vcftools,window=analysis.get("window_size",100000),step=analysis.get("window_step",50000),prefix=lambda wc:outdir+"/analysis/diversity/"+wc.population
  run:
   os.makedirs(os.path.dirname(str(log)),exist_ok=True);os.makedirs(os.path.dirname(str(output.pi)),exist_ok=True)
   script=os.path.join(os.path.dirname(str(output.pi)),f"vcftools_diversity_{wildcards.population}_{time.strftime('%Y%m%d_%H%M%S')}.sh")
   cmds=[[params.vcftools,"--gzvcf",input.vcf,"--keep",input.popfile,"--window-pi",str(params.window),"--window-pi-step",str(params.step),"--out",params.prefix],[params.vcftools,"--gzvcf",input.vcf,"--keep",input.popfile,"--TajimaD",str(params.window),"--out",params.prefix]]
   with open(str(log),"w") as h:h.write("")
   with open(script,"w") as h:h.write("#!/usr/bin/env bash\nset -euo pipefail\n"+"\n".join(" ".join(str(x) for x in c) for c in cmds)+"\n")
   shell(f"bash {script} >> {log} 2>&1")
if analysis.get("fst",True) and len(pops)>1:
 rule vcftools_population_fst:
  input: vcf=vcf,tbi=vcf+".tbi",popfiles=[outdir+"/metadata/populations/"+p+".txt" for p in pops]
  output: outdir+"/analysis/fst/{population_a}_vs_{population_b}.weir.fst"
  log: logdir+"/analysis/vcftools_fst_{population_a}_vs_{population_b}.log"
  conda: "vcftools_population.yaml"
  params: vcftools=vcftools,pop_a=lambda wc:outdir+"/metadata/populations/"+wc.population_a+".txt",pop_b=lambda wc:outdir+"/metadata/populations/"+wc.population_b+".txt",prefix=lambda wc:outdir+"/analysis/fst/"+wc.population_a+"_vs_"+wc.population_b
  run:
   os.makedirs(os.path.dirname(str(log)),exist_ok=True);os.makedirs(os.path.dirname(str(output[0])),exist_ok=True)
   cmd=[params.vcftools,"--gzvcf",input.vcf,"--weir-fst-pop",params.pop_a,"--weir-fst-pop",params.pop_b,"--out",params.prefix]
   script=os.path.join(os.path.dirname(str(output[0])),f"vcftools_fst_{wildcards.population_a}_vs_{wildcards.population_b}_{time.strftime('%Y%m%d_%H%M%S')}.sh")
   with open(str(log),"w") as h:h.write("")
   with open(script,"w") as h:h.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(str(x) for x in cmd)+"\n")
   shell(f"bash {script} >> {log} 2>&1")
