include: "../../common/common.smk"

# GO/KEGG enrichment over the scRNAseq pseudobulk output tree
# (layout-C sub-tool of modules/function; engine: bin/go-kegg.r).
#
# Input contract (Omics pseudobulk output):
#   <root>/<tissue>/<counter>/deg_pseudobulk{,_gene,_te}/
#       deg_<tissue>/<contrast>/<cell_type>/<contrast>_<cell_type>_DEG.csv
#   (columns: gene, cell_type, comparison, baseMean, log2FC, lfcSE, stat,
#    pval, pval_adj, significant)
# Output (per unit): <outdir>/function/go_kegg/{counter}/{mode}/{tissue}/
#   {contrast}/{cell_type}/go_up.csv, go_down.csv, kegg_up.csv, kegg_down.csv,
#   go_back_to_back.<fmt>, kegg_back_to_back.<fmt>, up_genes.txt, down_genes.txt

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"

_pb = (config.get("Params", {}).get("function", {}) or {}).get("pseudobulk_go", {}) or {}
pb_root = _pb.get("root") or indir
pb_modes = _pb.get("modes") or ["gene"]
pb_species = _pb.get("species") or ""
pb_lfc_cut = _pb.get("lfc_cut", 1)
pb_p_cut = _pb.get("p_cut", 0.05)
pb_top = _pb.get("top", 10)
pb_format = _pb.get("format") or "png"

_MODE_DIRS = {
    "all": "deg_pseudobulk",
    "gene": "deg_pseudobulk_gene",
    "te": "deg_pseudobulk_te",
}

# Unit discovery is a RUNTIME checkpoint: in workflow runs the tree only
# exists after scanpy_pseudobulk, so the parse-time glob of a standalone run
# cannot be used. scan_unit_lines() re-implements the same match at
# checkpoint time. `prerequisites` = the per-(tissue, counter) pseudobulk
# report htmls (ordering edge; empty for standalone runs on an existing tree).
prerequisites = config.get("prerequisites") or []
import glob as _glob  # noqa: E402


def scan_unit_lines():
    """List 'counter/mode/tissue/contrast/cell_type' units under pb_root."""
    mode_of = {v: k for k, v in _MODE_DIRS.items()}
    pat = os.path.join(pb_root, "*", "*", "deg_pseudobulk*",
                       "deg_*", "*", "*", "*_DEG.csv")
    units = set()
    for path in _glob.glob(pat):
        rel = os.path.relpath(path, pb_root).split(os.sep)
        if len(rel) != 7:
            continue
        tissue, counter, mode_dir, deg_tissue, contrast, cell_type, fname = rel
        if deg_tissue != f"deg_{tissue}":
            continue
        if fname != f"{contrast}_{cell_type}_DEG.csv":
            continue
        mode = mode_of.get(mode_dir)
        if mode is None or mode not in pb_modes:
            continue
        units.add(f"{counter}/{mode}/{tissue}/{contrast}/{cell_type}")
    return sorted(units)


checkpoint function_pseudobulk_go_index:
    """Discover enrichment units under the pseudobulk tree (after prereqs)."""
    input:
        prerequisites
    output:
        units = outdir + "/function/go_kegg/units.txt"
    log:
        logdir + "/function/pseudobulk_go_index.log"
    threads: 1
    resources:
        mem_mb=rule_res("function_pseudobulk_go").get("mem_mb", 8000),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("function_pseudobulk_go_index", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            units = scan_unit_lines()
            os.makedirs(os.path.dirname(str(output.units)), exist_ok=True)
            with open(str(output.units), "w") as fh:
                fh.write("\n".join(units) + ("\n" if units else ""))
            rule_logger.info(
                f"Pseudobulk GO/KEGG unit discovery at {current_time}: "
                f"{len(units)} unit(s)")
        except Exception as e:
            rule_logger.error(f"Error during unit discovery: {e}")
            raise RuntimeError(
                f"Error during unit discovery: {e}, you can check the log "
                f"file {log_path} for more details.")


def _go_out(name):
    """Output path under the per-unit go_kegg dir (fixed layout)."""
    return outdir + ("/function/go_kegg/{counter}/{mode}/{tissue}"
                     "/{contrast}/{cell_type}/" + name)


wildcard_constraints:
    mode="(all|gene|te)"


rule function_pseudobulk_go:
    """GO/KEGG enrichment for one pseudobulk DEG table (up and down split)."""
    input:
        deg=lambda wildcards: os.path.join(
            pb_root, wildcards.tissue, wildcards.counter,
            _MODE_DIRS[wildcards.mode], f"deg_{wildcards.tissue}",
            wildcards.contrast, wildcards.cell_type,
            f"{wildcards.contrast}_{wildcards.cell_type}_DEG.csv")
    output:
        go_up=_go_out("go_up.csv"),
        go_down=_go_out("go_down.csv"),
        kegg_up=_go_out("kegg_up.csv"),
        kegg_down=_go_out("kegg_down.csv"),
        go_plot=_go_out(f"go_back_to_back.{pb_format}"),
        kegg_plot=_go_out(f"kegg_back_to_back.{pb_format}"),
        up_genes=_go_out("up_genes.txt"),
        down_genes=_go_out("down_genes.txt"),
    log:
        logdir + "/function/{counter}/{mode}/{tissue}/{contrast}/{cell_type}.go_kegg.log"
    threads: rule_res("function_pseudobulk_go").get("threads", 1)
    resources:
        mem_mb=rule_res("function_pseudobulk_go").get("mem_mb", 8000),
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml", "function_pseudobulk_go")
    params:
        go_kegg_script=ROOT_DIR + "/modules/function/bin/go-kegg.r",
        species=pb_species,
        lfc_cut=pb_lfc_cut,
        p_cut=pb_p_cut,
        top=pb_top,
        format=pb_format,
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("function_pseudobulk_go", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start GO/KEGG enrichment for {wildcards.mode}/{wildcards.tissue}/"
                f"{wildcards.contrast}/{wildcards.cell_type} at {current_time}")
            # Species is wired from config, never guessed (wrong species silently
            # produces plausible-looking nonsense).
            if params.species not in ("human", "mouse", "macaque"):
                raise RuntimeError(
                    "Params.function.pseudobulk_go.species must be one of "
                    "human/mouse/macaque; set it in the module config JSON.")
            if not os.path.isfile(input.deg):
                raise FileNotFoundError(
                    f"DEG table not found: {input.deg} — expected the pseudobulk "
                    "layout <root>/deg_pseudobulk*/deg_<tissue>/<contrast>/"
                    "<cell_type>/<contrast>_<cell_type>_DEG.csv")
            sample_outdir = os.path.dirname(str(output.go_up))
            os.makedirs(sample_outdir, exist_ok=True)

            cmd = [
                "Rscript", params.go_kegg_script,
                "-i", input.deg,
                "-o", sample_outdir,
                "-s", params.species,
                "--gene-col", "gene",
                "--value-col", "log2FC",
                "--p-col", "pval_adj",
                "--lfc-cut", str(params.lfc_cut),
                "--p-cut", str(params.p_cut),
                "--top", str(params.top),
                "--format", str(params.format),
            ]
            # dual-write receipts: timestamped copy in RECEIPT_DIR (flat pool),
            # plain copy beside the results
            script = os.path.join(
                RECEIPT_DIR,
                f"function_pseudobulk_go_{wildcards.counter}_{wildcards.mode}_"
                f"{wildcards.tissue}_{wildcards.contrast}_"
                f"{wildcards.cell_type}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "function_pseudobulk_go.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
                f.write(f"export MKL_NUM_THREADS={threads}\n")
                f.write(f"export NUMBA_NUM_THREADS={threads}\n")
                f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "function pseudobulk GO/KEGG for {wildcards.mode}/'
                        f'{wildcards.tissue}/{wildcards.contrast}/'
                        f'{wildcards.cell_type} at {current_time} completed '
                        f'successfully"\n')
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during GO/KEGG enrichment: {e}")
            raise RuntimeError(
                f"Error occurred during GO/KEGG enrichment for {wildcards.mode}/"
                f"{wildcards.tissue}/{wildcards.contrast}/{wildcards.cell_type}: "
                f"{e}, you can check the log file {log_path} for more details.")


def _go_targets(wildcards):
    """go_up.csv of every discovered unit (resolved after the checkpoint)."""
    units_file = checkpoints.function_pseudobulk_go_index.get(**wildcards).output.units
    lines = [ln.strip() for ln in open(units_file) if ln.strip()]
    return [outdir + f"/function/go_kegg/{u}/go_up.csv" for u in lines]


rule function_pseudobulk_go_result:
    """Aggregates every discovered unit; writes the sentinel node.py registers."""
    input:
        _go_targets
    output:
        outdir + "/function/go_kegg/enrichment.done"
    log:
        logdir + "/function/function_pseudobulk_go_result.log"
    run:
        with open(str(output), "w") as fh:
            fh.write(time.strftime("%Y%m%d_%H%M%S", time.localtime()) + "\n")
