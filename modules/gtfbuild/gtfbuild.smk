include: "../common/common.smk"

outdir = config.get("outdir", "results/gtfbuild")
logdir = config.get("logdir", "log")
ROOT_DIR = config.get("ROOT_DIR", ".")
env = config.get("env", {})
genome_cfg = config.get("genome", {})
genome_name = config.get("genome_name") or "genome"
p = config.get("Params", {}).get("gtfbuild", {})

GTF_OUT = outdir + f"/{genome_name}.gtf"


rule gtfbuild_run:
    """Aggregate target: projected gene annotation GTF + QC report."""
    input:
        gtf = GTF_OUT,
        qc = outdir + f"/{genome_name}.gtf_qc.tsv",


rule gtfbuild_liftoff:
    """Project a source GTF (e.g. GENCODE) onto a target assembly with LiftOff.

    Same-species projection (e.g. GRCm39 -> T2T_mhaESC_v1.5) is the standard
    way to obtain a gene annotation for an assembly without an official GTF.
    Output GTF keeps gene_id/gene_name, so downstream peak-annotation and
    TSS-based tools work unchanged. Mark the provenance in the file name.
    """
    input:
        target_fasta = lambda wildcards: _require(genome_cfg.get("fasta"), "genome.fasta", "gtfbuild"),
        source_fasta = lambda wildcards: _require(genome_cfg.get("source_fasta"), "genome.source_fasta", "gtfbuild"),
        source_gtf = lambda wildcards: _require(genome_cfg.get("source_gtf"), "genome.source_gtf", "gtfbuild"),
    output:
        gff3 = outdir + f"/{genome_name}.lifted.gff3",
        unmapped = outdir + f"/{genome_name}.unmapped.gff3",
        gtf = GTF_OUT,
    log:
        logdir + "/gtfbuild_liftoff.log"
    threads: p.get("threads", 32)
    conda:
        "gtfbuild.yaml"
    container:
        sif("gtfbuild.yaml")
    params:
        coverage = p.get("coverage", 0.5),
        copies = "-copies" if p.get("copies", False) else "",
        extra = p.get("extra_args", ""),
        inter_dir = outdir + f"/{genome_name}.liftoff_intermediate",
        add_genes = ROOT_DIR + "/modules/gtfbuild/bin/gtfbuild_add_gene_features.py",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("gtfbuild_liftoff", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start LiftOff projection for {genome_name} at {current_time}")
            sample_outdir = os.path.dirname(str(output.gff3))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"gtfbuild_liftoff_{genome_name}_{current_time}.sh")
            lift_cmd = [
                "liftoff",
                "-g", str(input.source_gtf),
                "-o", str(output.gff3),
                "-u", str(output.unmapped),
                "-p", str(threads),
                "-a", str(params.coverage),
                "-dir", str(params.inter_dir),
                str(params.copies), str(params.extra),
                str(input.target_fasta), str(input.source_fasta),
            ]
            gffread_cmd = [
                "gffread", str(output.gff3), "-T", "-o", str(output.gtf),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(c for c in lift_cmd if c) + "\n")
                f.write(" ".join(gffread_cmd) + "\n")
                f.write(f'python {params.add_genes} --gtf {output.gtf} --out {output.gtf}.tmp\n')
                f.write(f"mv {output.gtf}.tmp {output.gtf}\n")
                f.write(f'echo "gtfbuild_liftoff for {genome_name} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during gtfbuild_liftoff for {genome_name}: {e}\n")
            raise RuntimeError(f"Error occurred during gtfbuild_liftoff for {genome_name}: {e}\n")


rule gtfbuild_qc:
    """Projection QC: mapping rate / complete models / extra copies."""
    input:
        gff3 = outdir + f"/{genome_name}.lifted.gff3",
        gtf = GTF_OUT,
        unmapped = outdir + f"/{genome_name}.unmapped.gff3",
        source_gtf = lambda wildcards: _require(genome_cfg.get("source_gtf"), "genome.source_gtf", "gtfbuild"),
    output:
        qc = outdir + f"/{genome_name}.gtf_qc.tsv",
    log:
        logdir + "/gtfbuild_qc.log"
    threads: 1
    conda:
        "gtfbuild.yaml"
    container:
        sif("gtfbuild.yaml")
    params:
        script = ROOT_DIR + "/modules/gtfbuild/bin/gtfbuild_qc.py",
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("gtfbuild_qc", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start gtfbuild QC for {genome_name} at {current_time}")
            sample_outdir = os.path.dirname(str(output.qc))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"gtfbuild_qc_{genome_name}_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--src-gtf", str(input.source_gtf),
                "--lifted-gff3", str(input.gff3),
                "--lifted-gtf", str(input.gtf),
                "--unmapped", str(input.unmapped),
                "--genome", genome_name,
                "--out", str(output.qc),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "gtfbuild_qc for {genome_name} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during gtfbuild_qc for {genome_name}: {e}\n")
            raise RuntimeError(f"Error occurred during gtfbuild_qc for {genome_name}: {e}\n")


# ---- Coordinate liftover of interval files (e.g. RepeatMasker BED) ----
# Repeat annotations must be projected through whole-genome alignment
# (synteny blocks), not per-sequence mapping: TE copies are near-identical
# and would land on paralogous loci. paftools.js liftover carries intervals
# inside long collinear alignments to their orthologous positions.

BED_SRC = genome_cfg.get("source_bed") or ""
BED_STEM = os.path.splitext(os.path.basename(BED_SRC))[0] if BED_SRC else "bed"
PAF_OUT = outdir + f"/{genome_name}.wga.paf"


rule gtfbuild_bed_run:
    """Aggregate target: lifted BED + projection stats."""
    input:
        bed = outdir + f"/{BED_STEM}.lifted.bed",
        stats = outdir + f"/{BED_STEM}.lifted.stats.tsv",


rule gtfbuild_wga:
    """Whole-genome alignment (PAF) between target and source assembly."""
    input:
        target_fasta = lambda wildcards: _require(genome_cfg.get("fasta"), "genome.fasta", "gtfbuild"),
        source_fasta = lambda wildcards: _require(genome_cfg.get("source_fasta"), "genome.source_fasta", "gtfbuild"),
    output:
        paf = PAF_OUT,
    log:
        logdir + "/gtfbuild_wga.log"
    threads: p.get("threads", 32)
    conda:
        "gtfbuild.yaml"
    container:
        sif("gtfbuild.yaml")
    params:
        preset = p.get("wga_preset", "asm5"),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("gtfbuild_wga", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start whole-genome alignment for {genome_name} at {current_time}")
            sample_outdir = os.path.dirname(str(output.paf))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"gtfbuild_wga_{genome_name}_{current_time}.sh")
            cmd = [
                "minimap2",
                "-cx", str(params.preset),
                "-t", str(threads),
                "--secondary=no",
                str(input.target_fasta), str(input.source_fasta),
                ">", str(output.paf),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "gtfbuild_wga for {genome_name} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during gtfbuild_wga for {genome_name}: {e}\n")
            raise RuntimeError(f"Error occurred during gtfbuild_wga for {genome_name}: {e}\n")


rule gtfbuild_liftover_bed:
    """Project a BED (e.g. RepeatMasker TEs) to the target assembly via WGA."""
    input:
        paf = PAF_OUT,
        source_bed = lambda wildcards: _require(genome_cfg.get("source_bed"), "genome.source_bed", "gtfbuild"),
    output:
        bed = outdir + f"/{BED_STEM}.lifted.bed",
        stats = outdir + f"/{BED_STEM}.lifted.stats.tsv",
    log:
        logdir + "/gtfbuild_liftover_bed.log"
    threads: 1
    conda:
        "gtfbuild.yaml"
    container:
        sif("gtfbuild.yaml")
    params:
        min_mapq = p.get("liftover_min_mapq", 5),
        min_aln = p.get("liftover_min_aln", 50000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("gtfbuild_liftover_bed", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start BED liftover for {genome_name} at {current_time}")
            sample_outdir = os.path.dirname(str(output.bed))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"gtfbuild_liftover_bed_{BED_STEM}_{current_time}.sh")
            cmd = [
                "paftools.js", "liftover",
                "-q", str(params.min_mapq),
                "-l", str(params.min_aln),
                str(input.paf), str(input.source_bed),
                ">", str(output.bed),
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(
                    f'in_n=$(grep -vc "^#" {input.source_bed} || true)\n'
                    f'out_n=$(grep -vc "^#" {output.bed} || true)\n'
                    f'awk -v s="{BED_STEM}" -v i="$in_n" -v o="$out_n" \'BEGIN{{'
                    f'printf "metric\\tvalue\\nsource_bed\\t%s\\nin_intervals\\t%s\\n'
                    f'out_intervals\\t%s\\nlift_rate\\t%.4f\\n", s, i, o, (i>0? o/i : 0)}}\' '
                    f'> {output.stats}\n'
                )
                f.write(f'echo "gtfbuild_liftover_bed for {BED_STEM} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during gtfbuild_liftover_bed for {BED_STEM}: {e}\n")
            raise RuntimeError(f"Error occurred during gtfbuild_liftover_bed for {BED_STEM}: {e}\n")
