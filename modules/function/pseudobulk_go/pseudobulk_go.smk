include: "../../common/common.smk"

# GO/KEGG enrichment over the scRNAseq pseudobulk output tree
# (layout-C sub-tool of modules/function; engine: bin/go-kegg.r).
#
# Input contract (Omics pseudobulk output):
#   <root>/<tissue>/pseudobulk/<counter>/deg_pseudobulk{,_gene,_te}/
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

# Parse-time discovery of (counter, mode, tissue, contrast, cell_type) units
# from the pseudobulk tree; all_pseudobulk_DEG.csv lives one level up and is
# not matched. The <tissue> dir and the deg_<tissue> dir must agree.
_w = glob_wildcards(
    pb_root + ("/{tissue}/pseudobulk/{counter}/{mode_dir}/deg_{dtissue}"
               "/{contrast}/{cell_type}/{stem}_DEG.csv"))
PB_UNITS = sorted({
    (_w.counter[i], mode, _w.tissue[i], _w.contrast[i], _w.cell_type[i])
    for i in range(len(_w.stem))
    for mode, mode_dir in _MODE_DIRS.items()
    if mode in pb_modes and _w.mode_dir[i] == mode_dir
    if _w.dtissue[i] == _w.tissue[i]
})


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
            pb_root, wildcards.tissue, "pseudobulk", wildcards.counter,
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
        sif("../function.yaml")
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


rule function_pseudobulk_go_result:
    """Panel rule: aggregates every discovered GO/KEGG unit."""
    input:
        [_go_out("go_up.csv").format(counter=co, mode=m, tissue=t,
                                     contrast=c, cell_type=ct)
         for (co, m, t, c, ct) in PB_UNITS]
