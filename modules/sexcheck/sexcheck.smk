include: "../common/common.smk"

import json

outdir = config.get("outdir", "output/QC/sexcheck")
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", logdir)
samples = config.get("samples", []) or {}
# sample_id -> BAM path (absolute or relative to cwd); the subworkflow
# resolves the alignment output for each sample (see subworkflow.md wiring)
sample_bams = config.get("sample_bams", {}) or {}

SEXCHECK_SCRIPT = os.path.join(ROOT_DIR, "modules", "sexcheck", "bin", "sexcheck.py")

# Sex-contig names from config genome.sex_chroms (single-genome wiring);
# absent/short lists yield None -> the script falls back to name sniffing
_sex_chroms = (config.get("genome", {}) or {}).get("sex_chroms") or []
X_CHROM = _sex_chroms[0] if len(_sex_chroms) > 0 else None
Y_CHROM = _sex_chroms[1] if len(_sex_chroms) > 1 else None


rule sexcheck_idxstats:
    """Per-sample samtools idxstats (chrom/length/mapped/unmapped) tagged
    with the sample id — single-sample rule."""
    input:
        bam=lambda wildcards: sample_bams[wildcards.sample_id],
    output:
        tsv=outdir + "/idxstats/{sample_id}/{sample_id}.idxstats.tsv",
    log:
        logdir + "/{sample_id}/sexcheck_idxstats.log"
    threads: rule_res("sexcheck_idxstats").get("threads", 1)
    resources:
        mem_mb=rule_res("sexcheck_idxstats").get("mem_mb", 8000),
    conda:
        "sexcheck.yaml"
    container:
        sif("sexcheck.yaml", "sexcheck_idxstats")
    params:
        samtools=config.get("Procedure", {}).get("samtools") or "samtools",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sexcheck_idxstats", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sexcheck idxstats for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.tsv))
            os.makedirs(sample_outdir, exist_ok=True)
            cmd = [
                params.samtools, "idxstats",
                str(input.bam),
                "|", "awk", "-F", "'\\t'", f"-v", f"s={wildcards.sample_id}",
                "'{print s\"\\t\"$1\"\\t\"$2\"\\t\"$3\"\\t\"$4}'",
                ">", str(output.tsv),
            ]
            script = os.path.join(RECEIPT_DIR, f"sexcheck_idxstats_{current_time}.sh")
            # plain_script must carry sample_id: all samples share the
            # idxstats/ outdir, a fixed name races under parallel jobs.
            plain_script = os.path.join(sample_outdir, f"sexcheck_idxstats.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(str(c) for c in cmd) + "\n")
                f.write(f'echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during sexcheck idxstats for sample {wildcards.sample_id}: {e}")
            raise RuntimeError(f"Error occurred during sexcheck idxstats for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.")


rule sexcheck:
    """Infer sample sex from sex-chromosome coverage (group rule: consumes
    every sample's idxstats table and calls XX/XY/unknown)."""
    input:
        idxstats=expand(outdir + "/idxstats/{sample_id}/{sample_id}.idxstats.tsv",
                        sample_id=list(sample_bams)),
    output:
        tsv=outdir + "/sexcheck.tsv",
        done=outdir + "/sexcheck.done",
    log:
        logdir_combine + "/sexcheck/sexcheck.log"
    threads: rule_res("sexcheck").get("threads", 1)
    resources:
        mem_mb=rule_res("sexcheck").get("mem_mb", 8000),
    conda:
        "sexcheck.yaml"
    container:
        sif("sexcheck.yaml", "sexcheck")
    params:
        script=SEXCHECK_SCRIPT,
        # sex-contig names from genome.references.<g>.sex_chroms (auto-filled
        # by node.py GENOME_SEX_CHROMS for NCBI-accession assemblies); null
        # falls back to the script's chrX/X + chrY/Y name sniffing.
        # Multi-genome workflows pass genome.sex_chroms_map + genome.samples
        # instead (single-genome workflows omit them).
        x_chrom=X_CHROM,
        y_chrom=Y_CHROM,
        sex_chroms_map=(config.get("genome", {}) or {}).get("sex_chroms_map") or {},
        sample_genomes=(config.get("genome", {}) or {}).get("samples") or {},
        # numeric-only passthrough: null placeholders (and test-mode path
        # injections) fall back to the script's own defaults
        ratio=(lambda v: v if isinstance(v, (int, float)) else None)(
            config.get("Params", {}).get("sexcheck", {}).get("ratio_xy_threshold")),
        y_min=(lambda v: v if isinstance(v, (int, float)) else None)(
            config.get("Params", {}).get("sexcheck", {}).get("y_min_frac")),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sexcheck", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start sexcheck at {current_time}")
            summary_dir = os.path.dirname(str(output.tsv))
            os.makedirs(summary_dir, exist_ok=True)
            merged = os.path.join(summary_dir, "idxstats_all.tsv")
            with open(merged, "w") as out:
                for f in input.idxstats:
                    with open(f) as fh:
                        out.write(fh.read())
            cmd = [
                "python3", params.script,
                "--tsv", merged,
                "--outdir", summary_dir,
                "--log", log_path,
            ]
            if params.ratio:
                cmd += ["--ratio-xy-threshold", str(params.ratio)]
            if params.y_min:
                cmd += ["--y-min-frac", str(params.y_min)]
            if params.x_chrom:
                cmd += ["--x-chrom", str(params.x_chrom)]
            if params.y_chrom:
                cmd += ["--y-chrom", str(params.y_chrom)]
            if params.sex_chroms_map:
                cmd += ["--sex-chroms-map", json.dumps(params.sex_chroms_map)]
            for s, g in (params.sample_genomes or {}).items():
                cmd += ["--sample-genome", f"{s}={g}"]
            script = os.path.join(RECEIPT_DIR, f"sexcheck_{current_time}.sh")
            plain_script = os.path.join(summary_dir, "sexcheck.sh")
            with open(script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(shlex.quote(str(c)) for c in cmd) + "\n")
                f.write('echo "script completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during sexcheck: {e}")
            raise RuntimeError(f"Error occurred during sexcheck: {e}, you can check the log file {log_path} for more details.")


rule sexcheck_result:
    """Result aggregation rule for subworkflow use rule import."""
    input:
        tsv=outdir + "/sexcheck.tsv",
        done=outdir + "/sexcheck.done",
    resources:
        mem_mb = rule_res("sexcheck_result").get("mem_mb", 8000),
