include: "../../common/common.smk"
import json  # noqa: E402  (common.smk 未导入而本模块需要的包)

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"

pseudobulk_go_params = config.get("Params", {}).get("function", {}).get("pseudobulk_go", {})

root = pseudobulk_go_params.get("root") or indir
modes = pseudobulk_go_params.get("modes") or ["all", "gene", "te"]
formats = pseudobulk_go_params.get("formats") or ["png"]
if isinstance(formats, str):
    formats = [formats]

mode_dirs = {
    "all": "deg_pseudobulk",
    "gene": "deg_pseudobulk_gene",
    "te": "deg_pseudobulk_te",
}
go_unit_dir = outdir + "/function/go_kegg/{counter}/{mode}/{tissue}/{contrast}/{cell_type}"
prerequisites = config.get("prerequisites") or []


def deg_table(wildcards):
    """DEG table of one unit (pseudobulk 输出契约;缺文件建图期即失败)."""
    return os.path.join(
        root, wildcards.tissue, wildcards.counter,
        mode_dirs[wildcards.mode], f"deg_{wildcards.tissue}",
        wildcards.contrast, wildcards.cell_type,
        f"{wildcards.contrast}_{wildcards.cell_type}_DEG.csv")


wildcard_constraints:
    mode="(all|gene|te)"


checkpoint function_pseudobulk_go_index:
    """Discover enrichment units under the pseudobulk tree (after prereqs)."""
    input:
        prerequisites
    output:
        units = outdir + "/function/go_kegg/units.txt"
    log:
        logdir + "/function/function_pseudobulk_go_index.log"
    threads: rule_res("function_pseudobulk_go_index").get("threads", 1)
    resources:
        mem_mb=rule_res("function_pseudobulk_go_index").get("mem_mb", 8000),
    params:
        python=config.get("Procedure", {}).get("python") or "python",
        script=ROOT_DIR + "/modules/function/pseudobulk_go/bin/discover_units.py",
        root=root,
        modes=modes,
        mode_dirs=mode_dirs,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("function_pseudobulk_go_index", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start pseudobulk GO/KEGG unit discovery at {current_time}")
            out_units = str(output.units)
            os.makedirs(os.path.dirname(out_units), exist_ok=True)
            cmd = [
                params.python, params.script,
                "--root", str(params.root),
                "--modes", json.dumps(params.modes),
                "--mode-dirs", json.dumps(params.mode_dirs),
                "--out", out_units,
            ]
            script = receipt_script(f"function_pseudobulk_go_index_{current_time}.sh")
            plain_script = os.path.join(os.path.dirname(out_units), "function_pseudobulk_go_index.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write('echo "pseudobulk GO/KEGG unit discovery was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error during unit discovery: {e}")
            raise RuntimeError(
                f"Error during unit discovery: {e}, you can check the log "
                f"file {log_path} for more details.")


rule function_pseudobulk_go:
    """GO/KEGG enrichment for one pseudobulk DEG table (up and down split)."""
    input:
        deg=deg_table
    output:
        go_up=go_unit_dir + "/go_up.csv",
        go_down=go_unit_dir + "/go_down.csv",
        kegg_up=go_unit_dir + "/kegg_up.csv",
        kegg_down=go_unit_dir + "/kegg_down.csv",
        go_plot=[go_unit_dir + "/go_back_to_back." + fmt for fmt in formats],
        kegg_plot=[go_unit_dir + "/kegg_back_to_back." + fmt for fmt in formats],
        up_genes=go_unit_dir + "/up_genes.txt",
        down_genes=go_unit_dir + "/down_genes.txt",
    log:
        logdir + "/function/{counter}/{mode}/{tissue}/{contrast}/{cell_type}/function_pseudobulk_go.log"
    benchmark: logdir + "/function/{counter}/{mode}/{tissue}/{contrast}/{cell_type}/function_pseudobulk_go.benchmark.tsv"
    threads: rule_res("function_pseudobulk_go").get("threads", 1)
    resources:
        mem_mb=rule_res("function_pseudobulk_go").get("mem_mb", 8000),
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml", "function_pseudobulk_go")
    params:
        rscript=config.get("Procedure", {}).get("Rscript") or "Rscript",
        go_kegg_script=ROOT_DIR + "/modules/function/bin/go-kegg.r",
        species=pseudobulk_go_params.get("species", "") or "",
        lfc_cut=pseudobulk_go_params.get("lfc_cut", 1),
        p_cut=pseudobulk_go_params.get("p_cut", 0.05),
        top=pseudobulk_go_params.get("top", 10),
        formats=formats,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("function_pseudobulk_go", log_file=log_path)
        unit = (f"{wildcards.counter}/{wildcards.mode}/{wildcards.tissue}"
                f"/{wildcards.contrast}/{wildcards.cell_type}")
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start GO/KEGG enrichment for {unit} at {current_time}")
            unit_outdir = os.path.dirname(str(output.go_up))
            os.makedirs(unit_outdir, exist_ok=True)

            cmd = [
                params.rscript, params.go_kegg_script,
                "-i", input.deg,
                "-o", unit_outdir,
                "-s", params.species,
                "--gene-col", "gene",
                "--value-col", "log2FC",
                "--p-col", "pval_adj",
                "--lfc-cut", str(params.lfc_cut),
                "--p-cut", str(params.p_cut),
                "--top", str(params.top),
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            script = receipt_script(
                f"function_pseudobulk_go_{wildcards.counter}_{wildcards.mode}_"
                f"{wildcards.tissue}_{wildcards.contrast}_"
                f"{wildcards.cell_type}_{current_time}.sh")
            plain_script = os.path.join(
                unit_outdir,
                f"function_pseudobulk_go_{wildcards.counter}_{wildcards.mode}_"
                f"{wildcards.tissue}_{wildcards.contrast}_"
                f"{wildcards.cell_type}.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "GO/KEGG enrichment for {unit} '
                        f'was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during GO/KEGG enrichment for {unit}: {e}")
            raise RuntimeError(
                f"Error occurred during GO/KEGG enrichment for {unit}: "
                f"{e}, you can check the log file {log_path} for more details.")


def go_kegg_targets(wildcards):
    """go_up.csv of every discovered unit (resolved after the checkpoint)."""
    units_file = checkpoints.function_pseudobulk_go_index.get(**wildcards).output.units
    with open(units_file) as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    return [outdir + f"/function/go_kegg/{u}/go_up.csv" for u in lines]


rule function_pseudobulk_go_done:
    """Write the sentinel node.py registers once every unit finished."""
    input:
        go_kegg_targets
    output:
        done = outdir + "/function/go_kegg/enrichment.done"
    log:
        logdir + "/function/function_pseudobulk_go_done.log"
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("function_pseudobulk_go_done", log_file=log_path)
        current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        with open(str(output.done), "w") as fh:
            fh.write(current_time + "\n")
        rule_logger.info(f"pseudobulk GO/KEGG enrichment completed at {current_time}")


rule function_pseudobulk_go_result:
    """Panel rule: every discovered unit's enrichment table plus the sentinel."""
    input:
        units=go_kegg_targets,
        done=outdir + "/function/go_kegg/enrichment.done"
