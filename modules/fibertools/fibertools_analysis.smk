"""Fiber-seq module for single-molecule chromatin accessibility analysis.

Provides rules for:
  - ft_predict_m6a: Predict m6A modifications from PacBio HiFi kinetics
  - ft_add_nucleosomes: Add nucleosome calls to Fiber-seq BAM
  - ft_fire: Call Fiber-seq Inferred Regulatory Elements (FIREs)
  - ft_extract: Extract Fiber-seq data to BED/bigBed format
  - ft_call_peaks: Call FIRE peaks using FDR-based peak calling

Requires fibertools-rs (ft) CLI tool.
See: https://fiberseq.github.io/
"""
include: "../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
samples = config.get("samples", [])
_bam_suffix = config.get("bam_suffix", ".fiberseq.fire.aligned.bam")

rule ft_extract:
    """Extract Fiber-seq data to BED format.

    Input: Fiber-seq BAM with m6A, nucleosome, and FIRE calls (aligned).
    Output: Separate compressed BED12 files for m6a, nuc, msp, cpg.
    """
    input:
        bam = indir + "/{sample_id}/{sample_id}" + _bam_suffix
    output:
        m6a = outdir + "/{sample_id}/{sample_id}.m6a.bed.gz",
        nuc = outdir + "/{sample_id}/{sample_id}.nuc.bed.gz",
        msp = outdir + "/{sample_id}/{sample_id}.msp.bed.gz",
        cpg = outdir + "/{sample_id}/{sample_id}.cpg.bed.gz",
    log:
        logdir + "/{sample_id}/ft_extract.log"
    threads: rule_res("ft_extract").get("threads", 8)
    resources:
        mem_mb = rule_res("ft_extract").get("mem_mb", 8000),
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml", "ft_extract")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft"
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("ft_extract", log_file=log_path)
        try:
            rule_logger.info(f"Extracting Fiber-seq data for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(RECEIPT_DIR, f"ft_extract_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "ft_extract.sh")
            cmd = [
                params.ft, "extract",
                "-r",
                "-t", str(threads),
                "--m6a", output.m6a,
                "--nuc", output.nuc,
                "--msp", output.msp,
                "-c", output.cpg,
                input.bam,
            ]
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Fiber-seq extraction completed for sample {wildcards.sample_id}"\n')
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_extract failed: {e}\n")
            raise RuntimeError(f"ft_extract failed: {e}\n")

def _call_peaks_mode():
    """Normalize Params.fibertools.call_peaks_mode ('fdr' | 'fire_frac').

    Single source shared by the input function (DAG-time gating) and the
    run block (command construction), per modules.md run 块禁止模块级全局变量.
    """
    return config.get("Params", {}).get("fibertools", {}).get("call_peaks_mode", "fdr")


def _ft_call_peaks_inputs(wildcards):
    """Declare every file ft_call_peaks consumes; missing ones fail at DAG
    build time instead of mid-job (modules.md: 消费的文件必须声明为 input).

    Returns a dict for ``unpack``: the aligned BAM plus, in fdr mode, the
    null-model references (chrom sizes + shuffle script) — optional-file
    pattern #1 from modules.md 可选文件的两种合规写法.
    """
    files = {
        "bam": indir + "/{sample_id}/{sample_id}" + _bam_suffix,
        "shuffle_script": ROOT_DIR + "/modules/fibertools/bin/shuffle_bed.py",
    }
    if _call_peaks_mode() == "fdr":
        # fdr mode shuffles fibers genome-wide; without chrom sizes the
        # null cannot be built — fail at DAG time naming the JSON key.
        files["fai"] = _require(
            config.get("genome", {}).get("fai"),
            "genome.fai", "ft_call_peaks",
        )
    return files


def _call_peaks_mode():
    """Normalize Params.fibertools.call_peaks_mode ('fdr' | 'fire_frac').

    Single source shared by the input function (DAG-time gating) and the
    run block (command construction), per modules.md run 块禁止模块级全局变量.
    """
    return config.get("Params", {}).get("fibertools", {}).get("call_peaks_mode", "fdr")


def _ft_call_peaks_inputs(wildcards):
    """Declare every file ft_call_peaks consumes; missing ones fail at DAG
    build time instead of mid-job (modules.md: 消费的文件必须声明为 input).

    Returns a dict for ``unpack``: the aligned BAM plus the fdr-mode null
    references (chrom sizes + engine script), added only when configured —
    optional-file pattern #1 from modules.md 可选文件的两种合规写法.
    """
    files = {
        "bam": indir + "/{sample_id}/{sample_id}" + _bam_suffix,
        "engine": ROOT_DIR + "/modules/fibertools/bin/call_peaks.py",
        "shuffle_script": ROOT_DIR + "/modules/fibertools/bin/shuffle_bed.py",
    }
    if _call_peaks_mode() == "fdr":
        # fdr mode shuffles fibers genome-wide; without chrom sizes the
        # null cannot be built — fail at DAG time naming the JSON key.
        files["fai"] = _require(
            config.get("genome", {}).get("fai"),
            "genome.fai", "ft_call_peaks",
        )
    return files


rule ft_call_peaks:
    """Call FIRE peaks using FDR-based peak calling.

    Input: Fiber-seq BAM with FIRE calls (+ genome.fai in fdr mode).
    Output: BED file with called FIRE peaks and the FDR threshold table.

    All support logic (aligned-read guard, fiber-span extraction, shuffled
    null, scratch cleanup) lives in bin/call_peaks.py — the rule only
    validates and invokes it (modules.md 核心原则).
    """
    input:
        unpack(_ft_call_peaks_inputs),
    output:
        peaks = outdir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        fdr_table = outdir + "/{sample_id}/{sample_id}.fire_fdr_table.tsv",
    log:
        logdir + "/{sample_id}/ft_call_peaks.log"
    threads: rule_res("ft_call_peaks").get("threads", 8)
    resources:
        mem_mb = rule_res("ft_call_peaks").get("mem_mb", 8000),
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml", "ft_call_peaks")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
        cfg = config.get("Params", {}).get("fibertools", {}),
        tmp_dir = config.get("tmp_dir") or os.path.join(outdir, "tmp"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("ft_call_peaks", log_file=log_path)
        try:
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            command_script = os.path.join(RECEIPT_DIR, f"ft_call_peaks_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"ft_call_peaks_{wildcards.sample_id}.sh")
            rule_logger.info(f"Calling FIRE peaks for sample {wildcards.sample_id}")
            cmd = [
                "python", str(input.engine),
                "--bam", str(input.bam),
                "--out-peaks", str(output.peaks),
                "--fdr-table", str(output.fdr_table),
                "--mode", params.cfg.get("call_peaks_mode", "fdr"),
                "--threads", str(threads),
                "--ft", params.ft,
                "--tmp-dir", params.tmp_dir,
                "--sample-id", wildcards.sample_id,
            ]
            if params.cfg.get("call_peaks_mode", "fdr") == "fdr":
                cmd += [
                    "--fai", str(input.fai),
                    "--shuffle-script", str(input.shuffle_script),
                    "--sd-cov", str(params.cfg.get("sd_cov", 5.0)),
                    "--max-fdr", str(params.cfg.get("max_fdr", 0.05)),
                ]
            else:
                cmd += ["--min-fire-frac", str(params.cfg.get("min_fire_frac", 0.1))]
            if params.cfg.get("haps", False):
                cmd.append("--haps")
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"export TMPDIR={params.tmp_dir}\n")
                f.write(" ".join(shlex.quote(str(x)) for x in cmd) + "\n")
                f.write(f'echo "FIRE peak calling completed for sample {wildcards.sample_id}"\n')
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_call_peaks failed: {e}\n")
            raise RuntimeError(f"ft_call_peaks failed: {e}\n")

rule ft_qc:
    """Collect QC metrics from a Fiber-seq BAM.

    Input: Fiber-seq BAM with m6A, nucleosome, MSP, and FIRE calls.
    Output: TSV file with QC metrics.
    """
    input:
        bam = indir + "/{sample_id}/{sample_id}" + _bam_suffix
    output:
        qc = outdir + "/{sample_id}/{sample_id}.qc.tsv",
    log:
        logdir + "/{sample_id}/ft_qc.log"
    threads: rule_res("ft_qc").get("threads", 8)
    resources:
        mem_mb = rule_res("ft_qc").get("mem_mb", 8000),
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml", "ft_qc")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
        use_acf = config.get("Params", {}).get("fibertools", {}).get("acf", False),
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("ft_qc", log_file=log_path)
        try:
            rule_logger.info(f"Collecting QC metrics for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(RECEIPT_DIR, f"ft_qc_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "ft_qc.sh")
            cmd = [
                params.ft, "qc",
                "-t", str(threads),
                input.bam,
                output.qc,
            ]
            if params.use_acf:
                cmd.append("--acf")
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "QC metrics collected for sample {wildcards.sample_id}"\n')
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_qc failed: {e}\n")
            raise RuntimeError(f"ft_qc failed: {e}\n")

wildcard_constraints:
    hap = "hap[12]",

rule ft_extract_fire:
    """Extract per-fiber FIRE elements to BED9 via `ft fire --extract`.

    Input: Fiber-seq BAM with FIRE calls (aligned).
    Output: BED of FIRE elements per fiber (consumed by downstream analysis).
    """
    input:
        bam = indir + "/{sample_id}/{sample_id}" + _bam_suffix
    output:
        fire = outdir + "/{sample_id}/{sample_id}.fire.bed.gz",
    log:
        logdir + "/{sample_id}/ft_extract_fire.log"
    threads: rule_res("ft_extract_fire").get("threads", 8)
    resources:
        mem_mb = rule_res("ft_extract_fire").get("mem_mb", 8000),
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml", "ft_extract_fire")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
        is_ont = config.get("Params", {}).get("fibertools", {}).get("ont", False),
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger("ft_extract_fire", log_file=log_path)
        try:
            rule_logger.info(f"Extracting FIRE elements for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(RECEIPT_DIR, f"ft_extract_fire_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "ft_extract_fire.sh")
            tmp_fire = str(output.fire) + ".unsorted.tmp.bed"
            cmd = [
                params.ft, "fire",
                "-t", str(threads),
                "--extract",
                input.bam,
                tmp_fire,
            ]
            if params.is_ont:
                cmd.extend(["--ont"])
            tmp_dir = config.get("tmp_dir") or os.path.join(outdir, "tmp")
            os.makedirs(tmp_dir, exist_ok=True)
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f'export TMPDIR="{tmp_dir}"\n')
                f.write(" ".join(cmd) + "\n")
                f.write(f"LC_ALL=C sort -k1,1 -k2,2n {tmp_fire} | bgzip -c > {output.fire}\n")
                f.write(f"tabix -f -p bed {output.fire}\n")
                f.write(f"rm -f {tmp_fire}\n")
                f.write(f'echo "FIRE extraction completed for sample {wildcards.sample_id}"\n')
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_extract_fire failed: {e}\n")
            raise RuntimeError(f"ft_extract_fire failed: {e}\n")

rule ft_split_hap:
    """Split an HP-tagged aligned BAM into haplotype BAMs with samtools."""
    input:
        bam = indir + "/{sample_id}/{sample_id}" + _bam_suffix
    output:
        bam = outdir + "/{sample_id}/{sample_id}.{hap}.bam",
        bai = outdir + "/{sample_id}/{sample_id}.{hap}.bam.bai",
    log:
        logdir + "/{sample_id}/ft_split_hap_{hap}.log"
    threads: rule_res("ft_split_hap").get("threads", 4)
    resources:
        mem_mb = rule_res("ft_split_hap").get("mem_mb", 8000),
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml", "ft_split_hap")
    params:
        hp_value = lambda wildcards: wildcards.hap[-1],
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger(f"ft_split_hap_{wildcards.hap}", log_file=log_path)
        try:
            rule_logger.info(f"Splitting {wildcards.hap} reads for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(RECEIPT_DIR, f"ft_split_hap_{wildcards.hap}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"ft_split_hap_{wildcards.hap}.sh")
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f'samtools view -b -@ {threads} -d HP:{params.hp_value} '
                        f'{input.bam} > {output.bam}\n')
                f.write(f'samtools index {output.bam}\n')
                f.write(f'echo "Haplotype split completed for sample {wildcards.sample_id}"\n')
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_split_hap failed: {e}\n")
            raise RuntimeError(f"ft_split_hap failed: {e}\n")

rule ft_extract_hap:
    """Extract m6A/nuc/MSP/CpG BEDs from a haplotype-split BAM."""
    input:
        bam = outdir + "/{sample_id}/{sample_id}.{hap}.bam",
        bai = outdir + "/{sample_id}/{sample_id}.{hap}.bam.bai",
    output:
        m6a = outdir + "/{sample_id}/{sample_id}.{hap}.m6a.bed.gz",
        nuc = outdir + "/{sample_id}/{sample_id}.{hap}.nuc.bed.gz",
        msp = outdir + "/{sample_id}/{sample_id}.{hap}.msp.bed.gz",
        cpg = outdir + "/{sample_id}/{sample_id}.{hap}.cpg.bed.gz",
    log:
        logdir + "/{sample_id}/ft_extract_{hap}.log"
    threads: rule_res("ft_extract_hap").get("threads", 8)
    resources:
        mem_mb = rule_res("ft_extract_hap").get("mem_mb", 8000),
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml", "ft_extract_hap")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
    run:
        log_path = str(log)
        open(log_path, 'w').close()
        rule_logger = setup_logger(f"ft_extract_{wildcards.hap}", log_file=log_path)
        try:
            rule_logger.info(f"Extracting {wildcards.hap} elements for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            command_script = os.path.join(RECEIPT_DIR, f"ft_extract_{wildcards.hap}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"ft_extract_{wildcards.hap}.sh")
            cmd = [
                params.ft, "extract",
                "-r",
                "-t", str(threads),
                "--m6a", output.m6a,
                "--nuc", output.nuc,
                "--msp", output.msp,
                "-c", output.cpg,
                input.bam,
            ]
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "Haplotype extraction completed for sample {wildcards.sample_id}"\n')
            shutil.copyfile(command_script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_extract_hap failed: {e}\n")
            raise RuntimeError(f"ft_extract_hap failed: {e}\n")


rule fibertools_analysis_result:
    """Panel rule: aggregate this module's main outputs as a dependency
    point for subworkflow `use rule` (modules.md 收尾规则)."""
    input:
        m6a = expand(outdir + "/{sample_id}/{sample_id}.m6a.bed.gz", sample_id=samples),
        nuc = expand(outdir + "/{sample_id}/{sample_id}.nuc.bed.gz", sample_id=samples),
        msp = expand(outdir + "/{sample_id}/{sample_id}.msp.bed.gz", sample_id=samples),
        cpg = expand(outdir + "/{sample_id}/{sample_id}.cpg.bed.gz", sample_id=samples),
        fire = expand(outdir + "/{sample_id}/{sample_id}.fire.bed.gz", sample_id=samples),
        peaks = expand(outdir + "/{sample_id}/{sample_id}.fire_peaks.bed", sample_id=samples),
        fdr_table = expand(outdir + "/{sample_id}/{sample_id}.fire_fdr_table.tsv", sample_id=samples),
        qc = expand(outdir + "/{sample_id}/{sample_id}.qc.tsv", sample_id=samples)
