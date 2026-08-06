
include: "../common/common.smk"
import os
import time

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
tool = config.get("Procedure", {}).get("PopLDdecay") or "PopLDdecay"
vcf = config.get("input_vcf") or outdir + "/variants/filtered/population_snps.vcf.gz"
populations = sorted(config.get("populations", {}))
analysis = config.get("Params", {}).get("analysis", {})

if analysis.get("ld_decay", False):
    rule PopLDdecay_population:
        input:
            vcf=vcf,
            tbi=vcf + ".tbi",
            popfile=lambda wc: outdir + "/metadata/populations/" + wc.population + ".txt"
        output:
            stat=outdir + "/analysis/ld_decay/{population}.stat.gz"
        log:
            logdir + "/analysis/PopLDdecay_{population}.log"
        conda:
            "PopLDdecay_population.yaml"
        container:
            sif("PopLDdecay_population.yaml")
        params:
            tool=tool,
            prefix=lambda wc: outdir + "/analysis/ld_decay/" + wc.population,
            max_dist=analysis.get("ld_max_distance", 500000)
        run:
            log_path = str(log)
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(str(output.stat)), exist_ok=True)
            script = os.path.join(os.path.dirname(str(output.stat)), f"PopLDdecay_{wildcards.population}_{time.strftime('%Y%m%d_%H%M%S')}.sh")
            cmd = [params.tool, "-InVCF", input.vcf, "-OutStat", params.prefix, "-SubPop", input.popfile, "-MaxDist", str(params.max_dist)]
            with open(log_path, "w") as handle:
                handle.write("")
            with open(script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write(" ".join(str(item) for item in cmd) + "\n")
            shell(f"bash {script} >> {log_path} 2>&1")
