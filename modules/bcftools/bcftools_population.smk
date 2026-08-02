
include: "../common/common.smk"
import os
import time

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
bcftools = config.get("Procedure", {}).get("bcftools") or "bcftools"
joint_vcf = config.get("input_vcf") or outdir + "/variants/joint/joint.vcf.gz"
filters = config.get("Params", {}).get("filtering", {})

rule bcftools_population_filter:
    input:
        vcf=joint_vcf,
        tbi=joint_vcf + ".tbi"
    output:
        vcf=outdir + "/variants/filtered/population_snps.vcf.gz",
        tbi=outdir + "/variants/filtered/population_snps.vcf.gz.tbi"
    log:
        logdir + "/variants/bcftools_population_filter.log"
    threads: 4
    conda:
        "bcftools_population.yaml"
    params:
        bcftools=bcftools,
        maf=filters.get("maf", 0.05),
        missing=1 - filters.get("max_missing", 0.9),
        qual=filters.get("min_qual", 30),
        min_dp=filters.get("min_dp", 8),
        max_dp=filters.get("max_dp", 200)
    run:
        log_path = str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(output.vcf), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            script = os.path.join(os.path.dirname(output.vcf), f"bcftools_population_filter_{current_time}.sh")
            cmd1 = [params.bcftools, "view", "-m2", "-M2", "-v", "snps", input.vcf]
            cmd2 = [params.bcftools, "filter", "-e", f"QUAL<{params.qual} || DP<{params.min_dp} || DP>{params.max_dp}"]
            cmd3 = [params.bcftools, "view", "-i", f"F_MISSING<={params.missing} && MAF>={params.maf}", "-Oz", "-o", output.vcf]
            with open(script, "w") as handle:
                handle.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                handle.write(" ".join(str(x) for x in cmd1) + " | " + " ".join(str(x) for x in cmd2) + " | " + " ".join(str(x) for x in cmd3) + "\n")
                handle.write(f"{params.bcftools} index -t {output.vcf}\n")
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as exc:
            with open(log_path, "a") as handle:
                handle.write(f"bcftools filtering failed: {exc}\n")
            raise
