
include: "../common/common.smk"
import os,time
outdir=config.get("outdir","output"); logdir=config.get("logdir","logs"); easy_sfs=config.get("Procedure",{}).get("easySFS") or "easySFS.py"; analysis=config.get("Params",{}).get("analysis",{})
if analysis.get("history",False):
 if not config.get("populations"): raise ValueError("easySFS requires populations")
 rule easySFS_population_history:
  input: vcf=config.get("input_vcf") or outdir+"/variants/filtered/population_snps.vcf.gz",tbi=(config.get("input_vcf") or outdir+"/variants/filtered/population_snps.vcf.gz")+".tbi",popmap=outdir+"/metadata/population.tsv"
  output: done=outdir+"/analysis/history/easySFS.done"
  log: logdir+"/analysis/easySFS_population.log"
  conda: "easySFS_population.yaml"
  params: tool=easy_sfs,outdir=outdir+"/analysis/history/easySFS",extra=analysis.get("history_config","")
  run:
   os.makedirs(os.path.dirname(str(log)),exist_ok=True);os.makedirs(params.outdir,exist_ok=True)
   script=os.path.join(params.outdir,f"easySFS_population_{time.strftime('%Y%m%d_%H%M%S')}.sh")
   cmd=[params.tool,"-i",input.vcf,"-p",input.popmap,"-o",params.outdir]
   if params.extra: cmd += str(params.extra).split()
   with open(str(log),"w") as h:h.write("")
   with open(script,"w") as h:h.write("#!/usr/bin/env bash\nset -euo pipefail\n"+" ".join(str(x) for x in cmd)+"\ntouch "+str(output.done)+"\n")
   shell(f"bash {script} >> {log} 2>&1")
