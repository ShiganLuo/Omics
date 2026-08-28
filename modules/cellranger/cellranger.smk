"""Cell Ranger module for single-cell RNA-seq analysis.

Provides three rules:
  - cellranger_ref:   Build reference genome from FASTA + GTF
  - cellranger_count: Align reads and generate gene-cell count matrix
  - cellranger_to_h5ad: Convert filtered matrix to h5ad for Scanpy

Cell Ranger must be available as a binary (not a conda package).
Set ``config["Procedure"]["cellranger"]`` to the cellranger executable path.
"""
include: "../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])
logdir_ref = config.get("logdir_ref", "logs/ref")
h5ad_outdir = config.get("h5ad_outdir", f"{outdir}/3_h5ad")
cellranger_ref_dir_provided = config.get("genome", {}).get("cellranger_ref_dir")
fasta = config.get("genome", {}).get("fasta") or ""
gtf = config.get("genome", {}).get("gtf") or ""
cellranger_ref_dir = outdir + "/cellranger_ref"
genome_name = config.get("Params", {}).get("cellranger", {}).get("mkref", {}).get("genome_name", "GRCh38")
cellranger_transcriptome_dir = cellranger_ref_dir + "/" + genome_name
cellranger_input_dict = config.get("cellranger_input_dict", {})


rule cellranger_ref:
    """Build Cell Ranger reference from Ensembl FASTA + GENCODE GTF."""
    input:
        fasta = fasta,
        gtf = gtf
    output:
        ref_dir = directory(cellranger_ref_dir)
    log:
        logdir_ref + "/cellranger_ref/cellranger_ref.log"
    threads: 16
    conda:
        "cellranger.yaml"
    container:
        sif("cellranger.yaml")
    params:
        script = os.path.join(ROOT_DIR, "modules", "cellranger", "bin", "cellranger_ref.py"),
        cellranger = config.get("Procedure", {}).get("cellranger") or "cellranger",
        genome_name = config.get("Params", {}).get("cellranger", {}).get("mkref", {}).get("genome_name", "GRCh38"),
        version = config.get("Params", {}).get("cellranger", {}).get("mkref", {}).get("version", "2024-A"),
        nthreads = 16,
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger(logger_name="cellranger_ref", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cellranger_ref at {current_time}")
            os.makedirs(str(output.ref_dir), exist_ok=True)
            command_script = os.path.join(str(output.ref_dir), f"cellranger_ref_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--fasta", input.fasta,
                "--gtf", input.gtf,
                "--output", str(output.ref_dir),
                "--genome", params.genome_name,
                "--version", params.version,
                "--cellranger", params.cellranger,
                "--nthreads", str(params.nthreads),
            ]
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(str(item) for item in cmd) + "\n")
                f.write(f'echo "Cell Ranger reference built at {str(output.ref_dir)} successfully"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"cellranger_ref failed: {e}\n")
            raise RuntimeError(f"cellranger_ref failed: {e}\n")

def get_input_for_cellranger_count(wildcards):
    """Get input files for cellranger_count rule."""
    in_dict = {}
    fastq_dir = cellranger_input_dict.get(wildcards.sample_id, {}).get("fastq_dir")    
    if not os.path.exists(fastq_dir):
        raise FileNotFoundError(f"FASTQ directory {fastq_dir} does not exist for sample {wildcards.sample_id}")
    if cellranger_ref_dir_provided and os.path.exists(cellranger_ref_dir_provided):
        in_dict["cellranger_ref_dir"] = cellranger_ref_dir_provided
    else:
        if fasta and gtf and os.path.exists(fasta) and os.path.exists(gtf):
            in_dict["cellranger_ref_dir"] = cellranger_transcriptome_dir
        else:
            raise ValueError(f"Neither transcriptome reference nor FASTA/GTF provided for sample {wildcards.sample_id}")
    in_dict["fastq_dir"] = fastq_dir
    return in_dict

rule cellranger_count:
    """Run cellranger count on a single sample."""
    input:
        unpack(get_input_for_cellranger_count)
    output:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        raw_matrix = directory(outdir + "/{sample_id}/raw_feature_bc_matrix"),
        filtered_matrix = directory(outdir + "/{sample_id}/filtered_feature_bc_matrix"),
        h5 = outdir + "/{sample_id}/filtered_feature_bc_matrix.h5"
    log:
        logdir + "/{sample_id}/cellranger_count.log"
    threads: 8
    conda:
        "cellranger.yaml"
    container:
        sif("cellranger.yaml")
    params:
        cellranger = config.get("Procedure", {}).get("cellranger") or "cellranger",
        create_bam = config.get("Params", {}).get("cellranger", {}).get("count", {}).get("create-bam") or False,
        nosecondary = config.get("Params", {}).get("cellranger", {}).get("count", {}).get("nosecondary") or True,
        sample_prefix = lambda wildcards:cellranger_input_dict.get(wildcards.sample_id, {}).get("sample_prefix")
    run:
        log_path = str(log)
        try:
            open(log_path, 'w').close()
            rule_logger = setup_logger(logger_name="cellranger_count", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start cellranger_count for sample {wildcards.sample_id} at {current_time}")
            sample_id = wildcards.sample_id
            sample_outdir = os.path.join(outdir, sample_id)
            # Clean up stale sample directory (cellranger fails if dir exists but isn't a valid pipestance)
            if os.path.exists(sample_outdir):
                import shutil
                rule_logger.warning(f"Removing stale sample directory: {sample_outdir}")
                shutil.rmtree(sample_outdir)
            # Don't create sample_outdir - cellranger will create it
            command_script = os.path.join(outdir, f"cellranger_count_{sample_id}_{current_time}.sh")
            cmd = [
                params.cellranger, "count",
                "--id", sample_id,
                "--transcriptome", input.cellranger_ref_dir,
                "--fastqs", input.fastq_dir,
                "--sample", params.sample_prefix,
                "--localcores", str(threads),
                "--create-bam", str(params.create_bam).lower(),
            ]
            if params.nosecondary:
                cmd.append("--nosecondary")
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(f"cd {outdir}\n")
                f.write(" ".join(cmd) + "\n")
                # Move outputs from outs/ to sample dir to match declared output paths
                f.write(f"outs=\"{sample_outdir}/outs\"\n")
                f.write(f"mv \"$outs/possorted_genome_bam.bam\" \"{output.bam}\"\n")
                f.write(f"mv \"$outs/raw_feature_bc_matrix\" \"{sample_outdir}/raw_feature_bc_matrix\"\n")
                f.write(f"mv \"$outs/filtered_feature_bc_matrix\" \"{sample_outdir}/filtered_feature_bc_matrix\"\n")
                f.write(f"mv \"$outs/filtered_feature_bc_matrix.h5\" \"{output.h5}\"\n")
                f.write(f'echo "Cell Ranger count completed for sample {sample_id}"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"cellranger_count failed: {e}\n")
            raise RuntimeError(f"cellranger_count failed: {e}\n")

rule cellranger_to_h5ad:
    """Convert Cell Ranger filtered matrix to h5ad."""
    input:
        bam = outdir + "/{sample_id}/{sample_id}.bam",
        filtered_matrix = directory(outdir + "/{sample_id}/filtered_feature_bc_matrix"),
    output:
        h5ad = h5ad_outdir + "/{sample_id}/{sample_id}_cellranger.h5ad"
    log:
        logdir + "/{sample_id}/cellranger_to_h5ad.log"
    threads: 1
    conda:
        "cellranger.yaml"
    container:
        sif("cellranger.yaml")
    params:
        script = os.path.join(ROOT_DIR, "modules", "cellranger", "bin", "cellranger_to_h5ad.py")
    run:
        log_path = str(log)
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            os.makedirs(os.path.dirname(output.h5ad), exist_ok=True)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.dirname(output.h5ad)
            command_script = os.path.join(sample_outdir, f"cellranger_to_h5ad_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--input", input.filtered_matrix,
                "--output", output.h5ad,
                "--sample-id", wildcards.sample_id
            ]
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(str(item) for item in cmd) + "\n")
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"cellranger_to_h5ad failed: {e}\n")
            raise RuntimeError(f"cellranger_to_h5ad failed: {e}\n")
