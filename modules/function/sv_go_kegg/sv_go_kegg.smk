include: "../../common/common.smk"


indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"

comparisons = config.get("comparisons") or {}
sg_formats = config.get("Params", {}).get("function", {}).get("sv_go_kegg", {}).get("formats") or ["png"]
if isinstance(sg_formats, str):
    sg_formats = [sg_formats]


rule sv_go_kegg:
    """GO/KEGG enrichment for one comparison's SV-affected gene list."""
    input:
        annotated_tab=indir + "/{comparison}/{comparison}_annotated.tab"
    output:
        genes=outdir + "/{comparison}/sv_go_kegg/sv_genes.txt",
        go_csv=outdir + "/{comparison}/sv_go_kegg/go_enrichment.csv",
        go_plot=[outdir + "/{comparison}/sv_go_kegg/go_enrichment." + fmt
                 for fmt in sg_formats],
        kegg_csv=outdir + "/{comparison}/sv_go_kegg/kegg_enrichment.csv",
        kegg_plot=[outdir + "/{comparison}/sv_go_kegg/kegg_enrichment." + fmt
                   for fmt in sg_formats],
        done=outdir + "/{comparison}/sv_go_kegg/sv_go_kegg.done",
    log:
        logdir + "/{comparison}/sv_go_kegg.log"
    threads: rule_res("sv_go_kegg").get("threads", 1)
    benchmark: logdir + "/{comparison}/sv_go_kegg.benchmark.tsv",
    resources:
        mem_mb=rule_res("sv_go_kegg").get("mem_mb", 8000),
    conda:
        "../function.yaml"
    container:
        sif("../function.yaml", "sv_go_kegg")
    params:
        script=ROOT_DIR + "/modules/function/sv_go_kegg/bin/sv_go-kegg.r",
        species=config.get("Params", {}).get("function", {}).get("sv_go_kegg", {}).get("species", "") or "",
        impact=config.get("Params", {}).get("function", {}).get("sv_go_kegg", {}).get("impact", "HIGH,MODERATE,LOW") or "HIGH,MODERATE,LOW",
        top=config.get("Params", {}).get("function", {}).get("sv_go_kegg", {}).get("top", 10),
        formats=sg_formats,
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_go_kegg", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start SV GO/KEGG enrichment for {wildcards.comparison} "
                f"at {current_time}")
            if params.species not in ("human", "mouse", "macaque"):
                raise RuntimeError(
                    "Params.function.sv_go_kegg.species must be one of "
                    "human/mouse/macaque; set it in the module config JSON.")
            if not os.path.isfile(input.annotated_tab):
                raise FileNotFoundError(
                    f"VEP-annotated SV tab not found: {input.annotated_tab} — "
                    "expected the modules/sv layout "
                    "<indir>/<comparison>/<comparison>_annotated.tab")
            sample_outdir = os.path.dirname(str(output.genes))
            os.makedirs(sample_outdir, exist_ok=True)

            cmd = [
                "Rscript", params.script,
                "-i", input.annotated_tab,
                "-o", sample_outdir,
                "-s", params.species,
                "--impact", params.impact,
                "--top", str(params.top),
            ]
            for fmt in params.formats:
                cmd += ["--format", str(fmt)]
            script = receipt_script(f"sv_go_kegg_{wildcards.comparison}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, "sv_go_kegg.sh")
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(f"export OMP_NUM_THREADS={threads}\n")
                f.write(" ".join(shlex.quote(str(item)) for item in cmd) + "\n")
                f.write(f'echo "SV GO/KEGG for {wildcards.comparison} '
                        f'at {current_time} completed successfully"\n')
            shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            with open(str(output.done), "w") as fh:
                fh.write(current_time + "\n")
        except Exception as e:
            rule_logger.error(f"Error during SV GO/KEGG enrichment: {e}")
            raise RuntimeError(
                f"Error during SV GO/KEGG enrichment for "
                f"{wildcards.comparison}: {e}, you can check the log file "
                f"{log_path} for more details.")
