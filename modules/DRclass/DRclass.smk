include: "../common/common.smk"
# DRclass: DNA(WES somatic) x RNA-seq 联合证据分类模块.
# 所有输入必须为具体文件路径(config: samples.<id>.xxx 显式声明, 无目录扫描).

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", {})
Procedure = config.get("Procedure", {})

# 指令派生名(变量名 = config 原名): 兼容嵌套 {default:..., references:{<g>:{...}}}
# 与扁平 {fasta:...} 两种 genome 形态; or 回退,无 if/else 分支.
genome = config.get("genome", {})
genome = genome.get("references", {}).get(genome.get("default", ""), {}) or genome


rule drclass_anchors:
    """M4 锚表: somatic PASS(DNA锚) ∪ RNA VCF 检出(RNA锚, 排除DNA锚)."""
    input:
        somatic=lambda w: samples[w.sample_id]["somatic_vcf"],
        rna=lambda w: samples[w.sample_id]["rna_vcf"],
    output:
        tsv=outdir + "/{sample_id}/anchors.tsv",
        bed=outdir + "/{sample_id}/anchors.bed",
    log:
        logdir + "/{sample_id}/drclass_anchors.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=Procedure.get("python3") or "python3",
    benchmark: logdir + "/{sample_id}/drclass_anchors.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_anchors").get("mem_mb", 2000),
    conda:
        "DRclass.yaml"
    container:
        sif("DRclass.yaml", "drclass_anchors")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("drclass_anchors", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start drclass anchors for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"drclass_anchors_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"drclass_anchors_{wildcards.sample_id}.sh")

            cmd = [
                params.py_proc, params.py, "anchors",
                "--somatic-vcf", str(input.somatic),
                "--rna-vcf", str(input.rna),
                "--out-tsv", str(output.tsv),
                "--out-bed", str(output.bed),
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "drclass anchors for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during drclass anchors for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during drclass anchors for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule drclass_pileup:
    """M3 RNA 定量: bcftools mpileup 对 RNA BAM 在锚位点计数 (-L anchors.bed)."""
    input:
        bam=lambda w: samples[w.sample_id]["rna_bam"],
        bai=lambda w: samples[w.sample_id]["rna_bam"].replace(".bam", ".bai"),
        fasta=genome["fasta"],
        fai=genome["fai_index"],
        bed=outdir + "/{sample_id}/anchors.bed",
    output:
        tsv=outdir + "/{sample_id}/rna.pileup.tsv",
    log:
        logdir + "/{sample_id}/drclass_pileup.log"
    params:
        bcftools=Procedure.get("bcftools") or "bcftools",
        min_mq=config.get("Params", {}).get("DRclass", {}).get("drclass_pileup", {}).get("min_mq", 20),
        min_bq=config.get("Params", {}).get("DRclass", {}).get("drclass_pileup", {}).get("min_bq", 30),
    benchmark: logdir + "/{sample_id}/drclass_pileup.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_pileup").get("mem_mb", 8000),
    threads: rule_res("drclass_pileup").get("threads", 4)
    conda:
        "DRclass.yaml"
    container:
        sif("DRclass.yaml", "drclass_pileup")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("drclass_pileup", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start drclass pileup for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"drclass_pileup_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"drclass_pileup_{wildcards.sample_id}.sh")

            mpileup_args = [
                str(params.bcftools), "mpileup", "-f", str(input.fasta), "-R", str(input.bed), "-d", "1000",
                "--min-MQ", str(params.min_mq), "--min-BQ", str(params.min_bq),
                "-a", "FORMAT/AD,FORMAT/DP", "-Ou", "-q", "0", "-Q", "0", str(input.bam),
            ]
            query_args = [
                str(params.bcftools), "query", "-f", r"%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n",
            ]
            mpileup_cmd = " ".join(shlex.quote(str(x)) for x in mpileup_args)
            query_cmd = " ".join(shlex.quote(str(x)) for x in query_args) + " > " + shlex.quote(str(output.tsv))
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                # 多命令回执(管道两段)只落一份 .sh,一处 shell() 判断
                f.write(mpileup_cmd + " | " + query_cmd + "\n")
                f.write(f'echo "drclass pileup for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(mpileup_args) + " | " + " ".join(query_args))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during drclass pileup for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during drclass pileup for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


def drclass_classify_inputs(wildcards):
    """声明全部输入文件;可选文件配置了才进 input(建图期存在性校验)."""
    s = samples[wildcards.sample_id]
    files = {
        "anchors": outdir + f"/{wildcards.sample_id}/anchors.tsv",
        "pileup": outdir + f"/{wildcards.sample_id}/rna.pileup.tsv",
        "somatic": s["somatic_vcf"],
        "rna": s["rna_vcf"],
    }
    if s.get("germline_vcf"):
        files["germline"] = s["germline_vcf"]
    for name, key in (("strain", "strain_vcf"), ("splice", "splice_bed"), ("editing", "editing_bed")):
        v = s.get(key) or []
        if v:
            files[name] = v if isinstance(v, (list, tuple)) else [v]
    return files


rule drclass_classify:
    """状态机判定 -> O1 dna_anchored.tsv + O2 rna_anchored.tsv."""
    input:
        unpack(drclass_classify_inputs),
    output:
        dna=outdir + "/{sample_id}/{sample_id}.dna_anchored.tsv",
        rna=outdir + "/{sample_id}/{sample_id}.rna_anchored.tsv",
    log:
        logdir + "/{sample_id}/drclass_classify.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=Procedure.get("python3") or "python3",
        min_af=config.get("Params", {}).get("DRclass", {}).get("drclass_classify", {}).get("min_af", 0.02),
        min_ad=config.get("Params", {}).get("DRclass", {}).get("drclass_classify", {}).get("min_ad", 4),
        min_dp=config.get("Params", {}).get("DRclass", {}).get("drclass_classify", {}).get("min_dp", 20),
        germline_af=config.get("Params", {}).get("DRclass", {}).get("drclass_classify", {}).get("germline_af", 0.20),
        germline_dp=config.get("Params", {}).get("DRclass", {}).get("drclass_classify", {}).get("germline_dp", 10),
        splice_window=config.get("Params", {}).get("DRclass", {}).get("drclass_classify", {}).get("splice_window", 5),
    benchmark: logdir + "/{sample_id}/drclass_classify.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_classify").get("mem_mb", 4000),
    conda:
        "DRclass.yaml"
    container:
        sif("DRclass.yaml", "drclass_classify")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("drclass_classify", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start drclass classify for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.dna))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"drclass_classify_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"drclass_classify_{wildcards.sample_id}.sh")

            cmd = [
                params.py_proc, params.py, "classify",
                "--sample", str(wildcards.sample_id),
                "--anchors", str(input.anchors),
                "--pileup", str(input.pileup),
                "--somatic-vcf", str(input.somatic),
                "--rna-vcf", str(input.rna),
            ]
            germline = getattr(input, "germline", "")
            if germline:
                cmd += ["--germline-vcf", str(germline)]
            for f in getattr(input, "strain", []):
                cmd += ["--strain-vcf", str(f)]
            for f in getattr(input, "splice", []):
                cmd += ["--splice-bed", str(f)]
            for f in getattr(input, "editing", []):
                cmd += ["--editing-bed", str(f)]
            cmd += [
                "--min-af", str(params.min_af), "--min-ad", str(params.min_ad), "--min-dp", str(params.min_dp),
                "--germline-af", str(params.germline_af), "--germline-dp", str(params.germline_dp),
                "--splice-window", str(params.splice_window),
                "--out-dna", str(output.dna), "--out-rna", str(output.rna),
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "drclass classify for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during drclass classify for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during drclass classify for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")


rule drclass_summary:
    """O3: 每样本标签计数 + 穿透率 + R4 质控, 全样本汇总."""
    input:
        tables=expand(outdir + "/{s}/{s}.dna_anchored.tsv", s=samples.keys()) +
               expand(outdir + "/{s}/{s}.rna_anchored.tsv", s=samples.keys()),
    output:
        tsv=outdir + "/DRclass_summary.tsv",
    log:
        logdir + "/drclass_summary.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=Procedure.get("python3") or "python3",
        max_r4=config.get("Params", {}).get("DRclass", {}).get("drclass_summary", {}).get("max_r4_fraction", 0.30),
    benchmark: logdir + "/drclass_summary.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_summary").get("mem_mb", 2000),
    conda:
        "DRclass.yaml"
    container:
        sif("DRclass.yaml", "drclass_summary")
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("drclass_summary", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start drclass summary at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"drclass_summary_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "drclass_summary.sh")

            cmd = [
                params.py_proc, params.py, "summary",
                "--tables",
            ] + [str(t) for t in input.tables] + [
                "--max-r4", str(params.max_r4),
                "--out", str(output.tsv),
            ]
            with open(script_path, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "drclass summary was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during drclass summary: {e}\n")
            raise RuntimeError(f"Error occurred during drclass summary: {e}, you can check the log file {log_path} for more details.\n")


rule drclass_result:
    input:
        summary=outdir + "/DRclass_summary.tsv",
        dna_anchored=expand(outdir + "/{s}/{s}.dna_anchored.tsv", s=samples.keys()),
        rna_anchored=expand(outdir + "/{s}/{s}.rna_anchored.tsv", s=samples.keys()),
