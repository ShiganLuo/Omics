"""HiFi genome assembly module (hifiasm).

Assembles PacBio HiFi BAM into primary-contig FASTA and runs assembly QC
(statistics + terminal telomere capture figures, rule hifiasm_qc). Downstream
modules (centromere, telomere assembly_scan/tidk, ...) consume the output file
at the layout contract {outdir}/{sample_id}/assembly/asm.bp.p_ctg.fa.
"""

include: "../common/common.smk"

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
bam_substring = config.get("bam_substring") or ""
formats = config.get("Params", {}).get("hifiasm", {}).get("qc", {}).get("formats") or ["png"]
if isinstance(formats, str):
    formats = [formats]

def get_input_for_hifiasm(wildcards):
    """Get input files for hifiasm assembly.

    Dict-mode (inputs.hifiasm.bam / inputs.hifiasm.pbi, per-sample) wins —
    e.g. Fiberseq's staged ``{raw}/{sample}/{sample}.bam`` layout; else the
    module's standalone flat shape A ``{indir}/{sample_id}[.bam_substring].
    bam`` (+ its colocated ``.bam.pbi`` index, the same naming contract
    shape A uses).
    """
    per = (config.get("inputs", {}) or {}).get("hifiasm", {}) or {}
    bam_map = per.get("bam", {}) or {}
    pbi_map = per.get("pbi", {}) or {}
    if wildcards.sample_id in bam_map:
        bam = bam_map[wildcards.sample_id]
        return {"bam": bam,
                "bai": pbi_map.get(wildcards.sample_id) or (bam + ".pbi")}
    in_dict = {}
    if bam_substring != "":
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}." + bam_substring + ".bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}." + bam_substring + ".bam.pbi")
    else:
        in_dict["bam"] = os.path.join(indir, f"{wildcards.sample_id}.bam")
        in_dict["bai"] = os.path.join(indir, f"{wildcards.sample_id}.bam.pbi")
    return in_dict

rule hifiasm_assemble:
    """Assemble HiFi reads with hifiasm."""
    input:
        unpack(get_input_for_hifiasm)
    output:
        fasta = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.fa",
        p_ctg_gfa = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.gfa",
        hap1_ctg_gfa = outdir + "/{sample_id}/assembly/asm.bp.hap1.p_ctg.gfa",
        hap2_ctg_gfa = outdir + "/{sample_id}/assembly/asm.bp.hap2.p_ctg.gfa",
        p_ctg_lowq = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.lowQ.bed",
        hap1_lowq = outdir + "/{sample_id}/assembly/asm.bp.hap1.p_ctg.lowQ.bed",
        hap2_lowq = outdir + "/{sample_id}/assembly/asm.bp.hap2.p_ctg.lowQ.bed",
        p_utg_lowq = outdir + "/{sample_id}/assembly/asm.bp.p_utg.lowQ.bed",
        r_utg_lowq = outdir + "/{sample_id}/assembly/asm.bp.r_utg.lowQ.bed",
        p_ctg_noseq = maybe_temp(outdir + "/{sample_id}/assembly/asm.bp.p_ctg.noseq.gfa", "hifiasm"),
        hap1_noseq = maybe_temp(outdir + "/{sample_id}/assembly/asm.bp.hap1.p_ctg.noseq.gfa", "hifiasm"),
        hap2_noseq = maybe_temp(outdir + "/{sample_id}/assembly/asm.bp.hap2.p_ctg.noseq.gfa", "hifiasm"),
        p_utg_gfa = maybe_temp(outdir + "/{sample_id}/assembly/asm.bp.p_utg.gfa", "hifiasm"),
        p_utg_noseq = maybe_temp(outdir + "/{sample_id}/assembly/asm.bp.p_utg.noseq.gfa", "hifiasm"),
        r_utg_gfa = maybe_temp(outdir + "/{sample_id}/assembly/asm.bp.r_utg.gfa", "hifiasm"),
        r_utg_noseq = maybe_temp(outdir + "/{sample_id}/assembly/asm.bp.r_utg.noseq.gfa", "hifiasm"),
        ec_bin = maybe_temp(outdir + "/{sample_id}/assembly/asm.ec.bin", "hifiasm"),
        ovlp_source = maybe_temp(outdir + "/{sample_id}/assembly/asm.ovlp.source.bin", "hifiasm"),
        ovlp_reverse = maybe_temp(outdir + "/{sample_id}/assembly/asm.ovlp.reverse.bin", "hifiasm"),
        hifi_fq = maybe_temp(outdir + "/{sample_id}/assembly/{sample_id}.hifi.fq.gz", "hifiasm"),
        hifi_fq_meta = maybe_temp(outdir + "/{sample_id}/assembly/{sample_id}.hifi.fq.gz.meta.json", "hifiasm"),
        stamp = maybe_temp(outdir + "/{sample_id}/assembly/asm.hifiasm.stamp", "hifiasm"),
    log:
        logdir + "/{sample_id}/hifiasm_assemble.log"
    benchmark:
        logdir + "/{sample_id}/hifiasm_assemble.benchmark.tsv"
    conda: "hifiasm.yaml"
    container:
        sif("hifiasm.yaml", "hifiasm_assemble")
    threads: rule_res("hifiasm_assemble").get("threads", 32)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "hifiasm_assemble", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 64000),
    params:
        script = os.path.join(ROOT_DIR, "modules/hifiasm/bin/hifiasm_assemble.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        hifiasm = config.get("Procedure", {}).get("hifiasm") or "hifiasm",
        samtools = config.get("Procedure", {}).get("samtools") or "samtools",
        prefix = outdir + "/{sample_id}/assembly/asm",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("hifiasm_assemble", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hifiasm assembly for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.fasta))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"hifiasm_assemble_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"hifiasm_assemble_{wildcards.sample_id}.sh")
            cmd = [
                params.python, params.script,
                "--bam", str(input.bam),
                "--out-prefix", params.prefix,
                "--out-fasta", str(output.fasta),
                "--out-fastq", str(output.hifi_fq),
                "--threads", str(threads),
                "--samtools", params.samtools,
                "--hifiasm", params.hifiasm,
            ]
            cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(cmd_str + "\n")
                f.write(f'echo "hifiasm assembly for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + cmd_str)
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during hifiasm assembly for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during hifiasm assembly for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

rule hifiasm_qc:
    """Assembly QC: statistics + terminal telomere capture visualization."""
    input:
        fasta = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.fa",
    output:
        stats_tsv = outdir + "/{sample_id}/assembly/qc/{sample_id}_assembly_stats.tsv",
        stats_txt = outdir + "/{sample_id}/assembly/qc/{sample_id}_assembly_stats.txt",
        tel_tsv = outdir + "/{sample_id}/assembly/qc/{sample_id}_telomere_ends.tsv",
        tel_txt = outdir + "/{sample_id}/assembly/qc/{sample_id}_telomere_qc.txt",
        overview_plots = [outdir + "/{sample_id}/assembly/qc/{sample_id}_assembly_overview." + fmt for fmt in formats],
        tel_plots = [outdir + "/{sample_id}/assembly/qc/{sample_id}_telomere_qc." + fmt for fmt in formats],
        map_plots = [outdir + "/{sample_id}/assembly/qc/{sample_id}_telomere_map." + fmt for fmt in formats],
    log:
        logdir + "/{sample_id}/hifiasm_qc.log"
    benchmark:
        logdir + "/{sample_id}/hifiasm_qc.benchmark.tsv"
    conda: "hifiasm.yaml"
    container:
        sif("hifiasm.yaml", "hifiasm_qc")
    threads: rule_res("hifiasm_qc").get("threads", 1)
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "hifiasm_qc", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    params:
        script = os.path.join(ROOT_DIR, "modules/hifiasm/bin/hifiasm_qc.py"),
        python = config.get("Procedure", {}).get("python") or "python",
        formats = formats,
        scan_length = (config.get("Params", {}).get("hifiasm", {})
                       .get("qc", {}).get("scan_length") or 100000),
        min_tel_len = (config.get("Params", {}).get("hifiasm", {})
                       .get("qc", {}).get("min_tel_len") or 1000),
        end_margin = (config.get("Params", {}).get("hifiasm", {})
                      .get("qc", {}).get("end_margin") or 1000),
        top_contigs = (config.get("Params", {}).get("hifiasm", {})
                       .get("qc", {}).get("top_contigs") or 30),
        expected_chrom_ends = (config.get("Params", {}).get("hifiasm", {})
                               .get("qc", {}).get("expected_chrom_ends") or 0),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("hifiasm_qc", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start hifiasm QC for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.stats_tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            script = receipt_script(f"hifiasm_qc_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"hifiasm_qc_{wildcards.sample_id}.sh")
            cmd = [
                params.python, params.script,
                "--fasta", input.fasta,
                "--sample_name", wildcards.sample_id,
                "--output_dir", sample_outdir,
                "--scan_length", str(params.scan_length),
                "--min_tel_len", str(params.min_tel_len),
                "--end_margin", str(params.end_margin),
                "--top-contigs", str(params.top_contigs),
            ]
            if params.expected_chrom_ends:
                cmd += ["--expected-chrom-ends", str(params.expected_chrom_ends)]
            for _f in params.formats:
                cmd += ["-F", _f]
            cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(cmd_str + "\n")
                f.write(f'echo "hifiasm qc for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + cmd_str)
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during hifiasm QC for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during hifiasm QC for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")

rule hifiasm_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        fasta = outdir + "/{sample_id}/assembly/asm.bp.p_ctg.fa"

rule hifiasm_qc_result:
    """Result aggregation for assembly QC outputs."""
    input:
        stats_txt = outdir + "/{sample_id}/assembly/qc/{sample_id}_assembly_stats.txt",
        tel_txt = outdir + "/{sample_id}/assembly/qc/{sample_id}_telomere_qc.txt",
