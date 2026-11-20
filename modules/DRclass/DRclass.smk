include: "../common/common.smk"
# DRclass: DNA(WES somatic) x RNA-seq 联合证据分类模块.
# 所有输入必须为具体文件路径(config: samples.<id>.xxx 显式声明, 无目录扫描).

outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

SAMPLES = config.get("samples", {})
PROC = config.get("Procedure", {})

# genome 兼容两种结构: 扁平 {fasta:...} (模块直跑) 或嵌套
# {default:..., references:{<g>:{fasta:...}}} (subworkflow 传入)
GENOME = config.get("genome", {})
if "fasta" not in GENOME and "references" in GENOME:
    _g = GENOME.get("default") or sorted(GENOME["references"].keys())[0]
    GENOME = GENOME["references"].get(_g, {})


def _s(sid):
    return SAMPLES[sid]

rule drclass_anchors:
    """M4 锚表: somatic PASS(DNA锚) ∪ RNA VCF 检出(RNA锚, 排除DNA锚)."""
    input:
        somatic=lambda w: _s(w.sample_id)["somatic_vcf"],
        rna=lambda w: _s(w.sample_id)["rna_vcf"],
    output:
        tsv=outdir + "/{sample_id}/anchors.tsv",
        bed=outdir + "/{sample_id}/anchors.bed",
    log:
        logdir + "/{sample_id}/drclass_anchors.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=PROC.get("python3") or "python3",
    benchmark: logdir + "/{sample_id}/drclass_anchors.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_anchors").get("mem_mb", 2000),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"drclass_anchors_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "drclass_anchors.sh")

            cmd = [
                params.py_proc, params.py, "anchors",
                "--somatic-vcf", str(input.somatic),
                "--rna-vcf", str(input.rna),
                "--out-tsv", str(output.tsv),
                "--out-bed", str(output.bed),
            ]
            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(shlex.join([str(x) for x in cmd]) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule drclass_pileup:
    """M3 RNA 定量: bcftools mpileup 对 RNA BAM 在锚位点计数 (-L anchors.bed)."""
    input:
        bam=lambda w: _s(w.sample_id)["rna_bam"],
        bai=lambda w: _s(w.sample_id)["rna_bam"].replace(".bam", ".bai"),
        fasta=GENOME["fasta"],
        fai=GENOME["fai_index"],
        bed=outdir + "/{sample_id}/anchors.bed",
    output:
        tsv=outdir + "/{sample_id}/rna.pileup.tsv",
    log:
        logdir + "/{sample_id}/drclass_pileup.log"
    params:
        bcftools=PROC.get("bcftools") or "bcftools",
        min_mq=config.get("Params", {}).get("drclass", {}).get("drclass_pileup", {}).get("min_mq", 20),
        min_bq=config.get("Params", {}).get("drclass", {}).get("drclass_pileup", {}).get("min_bq", 30),
    benchmark: logdir + "/{sample_id}/drclass_pileup.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_pileup").get("mem_mb", 8000),
    threads: rule_res("drclass_pileup").get("threads", 4)
    run:
        log_path = str(log)
        open(log_path, "a").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"drclass_pileup_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "drclass_pileup.sh")

            mpileup_cmd = shlex.join([
                str(params.bcftools), "mpileup", "-f", str(input.fasta), "-R", str(input.bed), "-d", "1000",
                "--min-MQ", str(params.min_mq), "--min-BQ", str(params.min_bq),
                "-a", "FORMAT/AD,FORMAT/DP", "-Ou", "-q", "0", "-Q", "0", str(input.bam),
            ])
            query_cmd = shlex.join([
                str(params.bcftools), "query", "-f", r"%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n",
            ]) + " > " + shlex.quote(str(output.tsv))
            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                # 多命令回执(管道两段)只落一份 .sh,一处 shell() 判断
                f.write(mpileup_cmd + " | " + query_cmd + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule drclass_classify:
    """状态机判定 -> O1 dna_anchored.tsv + O2 rna_anchored.tsv."""
    input:
        anchors=outdir + "/{sample_id}/anchors.tsv",
        pileup=outdir + "/{sample_id}/rna.pileup.tsv",
        somatic=lambda w: _s(w.sample_id)["somatic_vcf"],
        rna=lambda w: _s(w.sample_id)["rna_vcf"],
        germline=lambda w: _s(w.sample_id).get("germline_vcf"),
        strain=lambda w: _s(w.sample_id).get("strain_vcf", []),
        splice=lambda w: _s(w.sample_id).get("splice_bed", []),
        editing=lambda w: _s(w.sample_id).get("editing_bed", []),
    output:
        dna=outdir + "/{sample_id}/{sample_id}.dna_anchored.tsv",
        rna=outdir + "/{sample_id}/{sample_id}.rna_anchored.tsv",
    log:
        logdir + "/{sample_id}/drclass_classify.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=PROC.get("python3") or "python3",
        min_af=config.get("Params", {}).get("drclass", {}).get("drclass_classify", {}).get("min_af", 0.02),
        min_ad=config.get("Params", {}).get("drclass", {}).get("drclass_classify", {}).get("min_ad", 4),
        min_dp=config.get("Params", {}).get("drclass", {}).get("drclass_classify", {}).get("min_dp", 20),
        germline_af=config.get("Params", {}).get("drclass", {}).get("drclass_classify", {}).get("germline_af", 0.20),
        germline_dp=config.get("Params", {}).get("drclass", {}).get("drclass_classify", {}).get("germline_dp", 10),
        splice_window=config.get("Params", {}).get("drclass", {}).get("drclass_classify", {}).get("splice_window", 5),
    benchmark: logdir + "/{sample_id}/drclass_classify.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_classify").get("mem_mb", 4000),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.dirname(str(output.dna))
            os.makedirs(sample_outdir, exist_ok=True)
            script_path = os.path.join(RECEIPT_DIR, f"drclass_classify_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "drclass_classify.sh")

            cmd = [
                params.py_proc, params.py, "classify",
                "--sample", str(wildcards.sample_id),
                "--anchors", str(input.anchors),
                "--pileup", str(input.pileup),
                "--somatic-vcf", str(input.somatic),
                "--rna-vcf", str(input.rna),
            ]
            if input.germline:
                cmd += ["--germline-vcf", str(input.germline)]
            strain_vcfs = input.strain if isinstance(input.strain, (list, tuple)) else ([input.strain] if input.strain else [])
            for f in strain_vcfs:
                cmd += ["--strain-vcf", str(f)]
            splice_beds = input.splice if isinstance(input.splice, (list, tuple)) else ([input.splice] if input.splice else [])
            for f in splice_beds:
                cmd += ["--splice-bed", str(f)]
            editing_beds = input.editing if isinstance(input.editing, (list, tuple)) else ([input.editing] if input.editing else [])
            for f in editing_beds:
                cmd += ["--editing-bed", str(f)]
            cmd += [
                "--min-af", str(params.min_af), "--min-ad", str(params.min_ad), "--min-dp", str(params.min_dp),
                "--germline-af", str(params.germline_af), "--germline-dp", str(params.germline_dp),
                "--splice-window", str(params.splice_window),
                "--out-dna", str(output.dna), "--out-rna", str(output.rna),
            ]
            with open(script_path, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(shlex.join([str(x) for x in cmd]) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise

rule drclass_summary:
    """O3: 每样本标签计数 + 穿透率 + R4 质控, 全样本汇总."""
    input:
        tables=expand(outdir + "/{s}/{s}.dna_anchored.tsv", s=SAMPLES.keys()) +
               expand(outdir + "/{s}/{s}.rna_anchored.tsv", s=SAMPLES.keys()),
    output:
        tsv=outdir + "/DRclass_summary.tsv",
    log:
        logdir + "/drclass_summary.log"
    params:
        py=os.path.join(ROOT_DIR, "modules", "DRclass", "bin", "drclass.py"),
        py_proc=PROC.get("python3") or "python3",
        max_r4=config.get("Params", {}).get("drclass", {}).get("drclass_summary", {}).get("max_r4_fraction", 0.30),
    benchmark: logdir + "/drclass_summary.benchmark.tsv",
    resources:
        mem_mb=rule_res("drclass_summary").get("mem_mb", 2000),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
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
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(shlex.join([str(x) for x in cmd]) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script_path, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            raise
