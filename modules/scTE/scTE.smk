include: "../common/common.smk"
ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
samples = config.get("samples", [])

_scTE_genome = config.get("Params", {}).get("scTE", {}).get("genome", "hg38")
_scTE_mode = config.get("Params", {}).get("scTE", {}).get("mode", "exclusive")
_scTE_index = f"{outdir}/index/{_scTE_genome}.{_scTE_mode}.idx"


rule scTE_build_index:
    """Build scTE genome index.

    Two modes:
      1. Resource-based: provide gene_gtf + te_bed → scTE_build -gene -te -g other
      2. Download-based: only genome name → scTE_build -g <genome>
    """
    output:
        index = _scTE_index
    log:
        logdir + "/scTE/scTE_build_index.log"
    threads: 1
    conda:
        "scTE.yaml"
    container:
        sif("scTE.yaml")
    params:
        scte_build = config.get("Procedure", {}).get("scTE_build") or "scTE_build",
        genome = _scTE_genome,
        mode = _scTE_mode,
        gene_gtf = config.get("genome", {}).get("gtf"),
        te_bed = config.get("genome", {}).get("te_bed"),
        index_dir = outdir + "/index",
    run:
        log_path = str(log)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger = setup_logger("scTE_build_index",log_file=log_path)
            out_prefix = os.path.join(params.index_dir, params.genome)
            command_script = os.path.join(
                params.index_dir, f"scTE_build_index_{current_time}.sh"
            )
            rule_logger.info(f"Building scTE index for genome {params.genome} with mode {params.mode}")
            if params.gene_gtf and params.te_bed:
                # Resource-based: build from local gene GTF + TE BED
                cmd = [
                    params.scte_build,
                    "-gene", params.gene_gtf,
                    "-te", params.te_bed,
                    "-o", out_prefix,
                    "-m", params.mode,
                    "-g", "other",
                ]
            else:
                # Download-based: fetch by genome name
                cmd = [
                    params.scte_build,
                    "-g", params.genome,
                    "-o", out_prefix,
                    "-m", params.mode,
                ]

            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(str(item) for item in cmd) + "\n")
                f.write(f'echo "scTE index built at {output.index} successfully"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"scTE_build_index failed: {e}\n")
            raise RuntimeError(f"scTE_build_index failed: {e}\n")

def get_input_for_scTE_quantify(wildcards):
    """Get input files for scTE_quantify rule."""
    in_dict = {}
    sample_id = wildcards.sample_id
    bam_path = f"{indir}/{sample_id}/{sample_id}.bam"
    in_dict["bam"] = bam_path
    scTE_index = config.get("genome", {}).get("scTE_index")
    if scTE_index and os.path.exists(scTE_index):
        in_dict["scTE_index"] = scTE_index
    else:
        in_dict["scTE_index"] = _scTE_index
    return in_dict
# ---------------------------------------------------------------------------
# Per-sample quantify
# ---------------------------------------------------------------------------
rule scTE_quantify:
    """Quantify TE expression from Cell Ranger BAM for one sample."""
    input:
        unpack(get_input_for_scTE_quantify)
    output:
        h5ad = outdir + "/{sample_id}/{sample_id}_scTE.h5ad"
    log:
        logdir + "/{sample_id}/scTE_quantify.log"
    threads: 8
    conda:
        "scTE.yaml"
    container:
        sif("scTE.yaml")
    params:
        script = os.path.join(ROOT_DIR, "modules", "scTE", "bin", "scTE_quantify.py"),
        cb_tag = config.get("Params", {}).get("scTE", {}).get("cb_tag", "CB"),
        umi_tag = config.get("Params", {}).get("scTE", {}).get("umi_tag", "UB"),
        scte_bin = config.get("Procedure", {}).get("scTE") or "scTE",
    run:
        log_path = str(log)
        try:
            rule_logger = setup_logger("scTE_quantify", log_file=log_path)
            rule_logger.info(f"Quantifying TE expression for sample {wildcards.sample_id}")
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            sample_outdir = os.path.dirname(output.h5ad)
            command_script = os.path.join(sample_outdir, f"scTE_quantify_{current_time}.sh")
            cmd = [
                "python", params.script,
                "--input", input.bam,
                "--output", output.h5ad,
                "--index", input.scTE_index,
                "--cb-tag", params.cb_tag,
                "--umi-tag", params.umi_tag,
                "--threads", str(threads),
                "--scte-bin", params.scte_bin,
                "--sample-id", wildcards.sample_id,
            ]
            with open(command_script, "w") as f:
                f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
                f.write(" ".join(str(item) for item in cmd) + "\n")
                f.write(f'echo "scTE quantification for sample {wildcards.sample_id} completed successfully"\n')
            shell(f"bash {command_script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"scTE_quantify failed: {e}\n")
            raise RuntimeError(f"scTE_quantify failed: {e}\n")

