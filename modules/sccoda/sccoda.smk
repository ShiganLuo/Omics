
include: "../common/common.smk"

ROOT_DIR = config.get("ROOT_DIR", ".")
indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
# reannotate: {counter: {tissue: override_json}} — forwarded by the
# subworkflow from scanpy's auto.cluster_annotations. Cross-module
# ordering guard only; see _reannotate_marker below.
reannotate = config.get("reannotate", {}) or {}

params = config.get("Params", {}).get("sccoda", {})
# tissue_samples: {tissue: [sample_id, ...]} — from node.py
tissue_samples = params.get("tissue_samples", {})
tissues = sorted(tissue_samples.keys())
counters = sorted(k for k, v in params.items()
                  if isinstance(v, dict) and k != "tissue_samples")
python = config.get("Procedure", {}).get("python") or "python"
script = os.path.join(ROOT_DIR, "modules", "sccoda", "bin", "sccoda.py")

def _sccoda_override(wildcards):
    """Path of the tissue's reannotate override JSON, [] when none."""
    entry = reannotate.get(wildcards.counter, {}) or {}
    if not isinstance(entry, dict):
        return []
    return entry.get(wildcards.tissue, "") or []


def _reannotate_marker(wildcards):
    """Marker file input enforcing auto -> reannotate -> sccoda ordering.

    The reannotate step overwrites the auto h5ad in place; sccoda depends on
    this marker so it never reads a pre-override auto.h5ad. Empty when no
    overrides are configured (the marker is then never produced)."""
    return (indir + f"/{wildcards.tissue}/{wildcards.tissue}_{wildcards.counter}_auto.h5ad.reannotated"
            if _sccoda_override(wildcards) else [])


# ---------------------------------------------------------------------------
# scCODA differential abundance over annotated cell types
# ---------------------------------------------------------------------------
rule sccoda_run:
    """scCODA differential abundance (pertpy): one job per tissue x counter,
    running every Params.sccoda.<counter>.compare entry scoped to its tissue.
    Consumes the scanpy_auto h5ad (obs: sample_id + cell_type); when cluster
    annotation overrides are configured for the tissue, the reannotate marker
    orders this rule after the in-place label fix."""
    input:
        h5ad = indir + "/{tissue}/{tissue}_{counter}_auto.h5ad",
        reannotated = _reannotate_marker
    output:
        report = outdir + "/{tissue}/sccoda/{counter}/sccoda_report.html",
        da_dir = directory(outdir + "/{tissue}/sccoda/{counter}/da")
    log:
        logdir + "/sccoda/{tissue}/sccoda_run_{counter}.log"
    threads: rule_res("sccoda_run").get("threads", 4)
    resources:
        mem_mb=rule_res("sccoda_run").get("mem_mb", 32000),
    conda:
        "sccoda.yaml"
    container:
        sif("sccoda.yaml")
    params:
        python = python,
        script = script,
        sccoda_cfg = lambda wildcards: params.get(wildcards.counter, {}),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("sccoda_run", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start scCODA differential abundance for tissue {wildcards.tissue} at {current_time}")
            out_root = os.path.dirname(str(output.report))
            os.makedirs(out_root, exist_ok=True)
            # dual-write receipts: timestamped copy in log/cmd, plain copy in outdir
            script = os.path.join(RECEIPT_DIR, f"sccoda_run_{wildcards.counter}_{wildcards.tissue}_{current_time}.sh")
            plain_script = os.path.join(out_root, f"sccoda_run_{wildcards.counter}_{wildcards.tissue}.sh")
            cfg = params.sccoda_cfg
            cmd = [params.python, params.script,
                   "--h5ad", f"{wildcards.tissue}={input.h5ad}",
                   "--out-dir", out_root,
                   "--ref-celltype", str(cfg.get("ref_celltype") or "automatic"),
                   "--formula", str(cfg.get("formula") or "condition"),
                   "--credible-threshold", str(cfg.get("credible_threshold", 0.90)),
                   "--num-samples", str(cfg.get("num_samples", 2000)),
                   "--num-warmup", str(cfg.get("num_warmup", 500))]
            # entries: tissue:ref_sample:treat_sample[=label] — every entry
            # is scoped to ONE tissue (sccoda.py rejects a compare whose
            # tissue has no --h5ad definition), so this per-tissue job only
            # takes the entries prefixed with its own wildcard. SAME filter
            # expression as node.py's registration gate.
            tissue_compares = [
                str(c) for c in (cfg.get("compare") or [])
                if str(c).split(":", 1)[0] == wildcards.tissue
            ]
            if not tissue_compares:
                raise RuntimeError(
                    f"Params.sccoda.{wildcards.counter}.compare "
                    f"has no entry for tissue '{wildcards.tissue}' — expected "
                    f"entries like '{wildcards.tissue}:ref_sample:"
                    f"treat_sample[=label]'")
            for cmp_ in tissue_compares:
                cmd += ["--compare", cmp_]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                # Cap numeric-library threads to the rule's thread budget
                # (config resources -> rule_res -> threads) so the process
                # honors what snakemake allocated instead of grabbing every
                # core; the value is injected at script generation, never
                # hardcoded here.
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "scCODA differential abundance for {wildcards.tissue} at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during scCODA differential abundance for tissue {wildcards.tissue}: {e}")
            raise RuntimeError(f"Error occurred during scCODA differential abundance for tissue {wildcards.tissue}: {e}")

rule sccoda_result:
    """Panel rule aggregating the module's main outputs."""
    input:
        expand(outdir + "/{tissue}/sccoda/{counter}/sccoda_report.html",
               tissue=tissues, counter=counters)
