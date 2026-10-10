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
    threads: 8
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
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
            command_script = os.path.join(sample_outdir, f"ft_extract_{current_time}.sh")
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
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_extract failed: {e}\n")
            raise RuntimeError(f"ft_extract failed: {e}\n")


rule ft_call_peaks:
    """Call FIRE peaks using FDR-based peak calling.

    Input: Fiber-seq BAM with FIRE calls.
    Output: BED file with called FIRE peaks.
    """
    input:
        bam = indir + "/{sample_id}/{sample_id}" + _bam_suffix
    output:
        peaks = outdir + "/{sample_id}/{sample_id}.fire_peaks.bed",
        fdr_table = outdir + "/{sample_id}/{sample_id}.fire_fdr_table.tsv",
    log:
        logdir + "/{sample_id}/ft_call_peaks.log"
    threads: 8
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
    params:
        ft = config.get("Procedure", {}).get("fibertools") or "ft",
        fdr_cfg = config.get("Params", {}).get("fibertools", {}),
        fai = config.get("genome", {}).get("fai", ""),
        shuffle_script = ROOT_DIR + "/modules/fibertools/bin/shuffle_bed.py",
    run:
        log_path = str(log)
        rule_logger = setup_logger("ft_call_peaks", log_file=log_path)
        try:
            sample_outdir = os.path.join(outdir, wildcards.sample_id)
            os.makedirs(sample_outdir, exist_ok=True)
            sorted_bam = os.path.join(sample_outdir, f"{wildcards.sample_id}.fire.sorted.bam")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            command_script = os.path.join(sample_outdir, f"ft_call_peaks_{current_time}.sh")
            rule_logger.info(f"Calling FIRE peaks for sample {wildcards.sample_id}")
            # Build ft call-peaks command (fdr mode needs a shuffled null BED)
            mode = params.fdr_cfg.get("call_peaks_mode", "fdr")
            if mode == "fdr" and not params.fai:
                raise RuntimeError(
                    "call_peaks_mode=fdr requires 'genome.fai' in the config JSON "
                    "(chrom sizes for the shuffled null)."
                )
            fibers_bed = os.path.join(sample_outdir, f"{wildcards.sample_id}.fibers.bed")
            shuffled_bed = os.path.join(sample_outdir, f"{wildcards.sample_id}.fibers.shuffled.bed")
            ft_cmd = [
                params.ft, "call-peaks",
                "-t", str(threads),
                "-o", str(output.peaks),
                "--sd-cov", str(params.fdr_cfg.get("sd_cov", 5.0)),
            ]
            if mode == "fdr":
                ft_cmd.extend([
                    "--max-fdr", str(params.fdr_cfg.get("max_fdr", 0.05)),
                    "--shuffled", shuffled_bed,
                    "--fdr-table-out", str(output.fdr_table),
                ])
            else:
                ft_cmd.extend(["--min-fire-frac", str(params.fdr_cfg.get("min_fire_frac", 0.1))])
            if params.fdr_cfg.get("haps", False):
                ft_cmd.append("--haps")
            ft_cmd.append(sorted_bam)
            ft_cmd_str = " ".join(ft_cmd)

            # Write shell script — all commands run inside conda/container env via shell()
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f'echo "Checking alignment for {wildcards.sample_id}"\n')
                f.write(f'aligned=$(samtools view -c -F 4 {input.bam})\n')
                f.write(f'echo "Aligned reads: $aligned"\n')
                f.write(f'if [ "$aligned" -eq 0 ]; then\n')
                f.write(f'    echo "BAM has no aligned reads — skipping ft call-peaks"\n')
                f.write(f'    printf "#chrom\\tstart\\tend\\tname\\tscore\\tstrand\\n" > {output.peaks}\n')
                f.write(f'    exit 0\n')
                f.write(f'fi\n')
                f.write(f'echo "Sorting BAM"\n')
                f.write(f'samtools sort -@ {threads} -o {sorted_bam} {input.bam}\n')
                f.write(f'samtools index {sorted_bam}\n')
                if mode == "fdr":
                    f.write(f'echo "Building shuffled null"\n')
                    f.write(f'samtools view -F 0x904 {sorted_bam} | '
                            f'awk \'BEGIN{{OFS="\\t"}}{{print $3, $4-1, $4-1+length($10), $1}}\' > {fibers_bed}\n')
                    f.write(f'python {params.shuffle_script} -i {fibers_bed} '
                            f'-g {params.fai} -o {shuffled_bed} --seed 42\n')
                f.write(f'echo "Calling FIRE peaks ({mode} mode)"\n')
                f.write(ft_cmd_str + "\n")
                if mode != "fdr":
                    f.write(f'printf "#FDR table not computed (fire_frac mode)\\n" > {output.fdr_table}\n')
                f.write(f'echo "FIRE peak calling completed for sample {wildcards.sample_id}"\n')

            shell(f"bash {command_script} >> {log_path} 2>&1")
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
    threads: 8
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
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
            command_script = os.path.join(sample_outdir, f"ft_qc_{current_time}.sh")
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
            shell(f"bash {command_script} >> {log_path} 2>&1")
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
        fire = outdir + "/{sample_id}/{sample_id}.fire.bed",
    log:
        logdir + "/{sample_id}/ft_extract_fire.log"
    threads: 8
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
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
            command_script = os.path.join(sample_outdir, f"ft_extract_fire_{current_time}.sh")
            cmd = [
                params.ft, "fire",
                "-t", str(threads),
                "--extract",
                input.bam,
                output.fire,
            ]
            if params.is_ont:
                cmd.extend(["--ont"])
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "FIRE extraction completed for sample {wildcards.sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
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
    threads: 4
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
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
            command_script = os.path.join(sample_outdir, f"ft_split_hap_{wildcards.hap}_{current_time}.sh")
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f'samtools view -b -@ {threads} -d HP:{params.hp_value} '
                        f'{input.bam} > {output.bam}\n')
                f.write(f'samtools index {output.bam}\n')
                f.write(f'echo "Haplotype split completed for sample {wildcards.sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
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
    threads: 8
    conda:
        "fibertools.yaml"
    container:
        sif("fibertools.yaml")
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
            command_script = os.path.join(sample_outdir, f"ft_extract_{wildcards.hap}_{current_time}.sh")
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
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"ft_extract_hap failed: {e}\n")
            raise RuntimeError(f"ft_extract_hap failed: {e}\n")
