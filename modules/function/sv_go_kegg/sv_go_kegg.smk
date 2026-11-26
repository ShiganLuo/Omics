include: "../../common/common.smk"

# GO/KEGG enrichment for SV-affected genes (sub-tool of modules/function;
# engine: bin/sv_go-kegg.r, clusterProfiler + org.{Hs,Mm,Mmu}.eg.db so
# MACAQUE is a first-class species — unlike the gseapy/Enrichr route in
# modules/sv, which only knows human/mouse).
#
# Mirrors modules/function/pseudobulk_go/pseudobulk_go.smk layout:
# one job per comparison over the VEP-annotated SV tab produced by
# sv_exp_specific.
#
# Input contract (modules/sv output):
#   <indir>/<comparison>/<comparison>_annotated.tab
# Output (per comparison):
#   <outdir>/<comparison>/{sv_genes.txt, go_enrichment.csv,
#                          go_enrichment.<fmt>, kegg_enrichment.csv,
#                          kegg_enrichment.<fmt>, sv_go_kegg.done}
#   (<fmt> per Params.function.sv_go_kegg.formats, default png)

indir = config.get("indir") or "input"
outdir = config.get("outdir") or "output"
logdir = config.get("logdir") or "log"

# 参数精确到 rule:参数挂在消费它的 rule 名下(Params.function.sv_go_kegg.<key>,
# rule 名 = 原组名;无 module 级共享回退,读点写 config 全量路径字面链)。
comparisons = config.get("comparisons") or {}
# ① parse 期读入 + 规范化(str 收成单元素 list;合法值由 schema 声明)。
# formats 为 parse 期变量(output 不能是函数,格式列表必须是模块级变量)。
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
        # ② 按格式展开:每种格式一个输出文件(base 不带扩展名)
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
        formats=sg_formats,   # ③ 转发
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("sv_go_kegg", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(
                f"Start SV GO/KEGG enrichment for {wildcards.comparison} "
                f"at {current_time}")
            # Species is wired from config, never guessed (wrong species
            # silently produces plausible-looking nonsense).
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
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
            # sentinel written by the rule (not the shell receipt) AFTER the
            # tool succeeded — the receipt's echo is not a success signal
            with open(str(output.done), "w") as fh:
                fh.write(current_time + "\n")
        except Exception as e:
            rule_logger.error(f"Error during SV GO/KEGG enrichment: {e}")
            raise RuntimeError(
                f"Error during SV GO/KEGG enrichment for "
                f"{wildcards.comparison}: {e}, you can check the log file "
                f"{log_path} for more details.")
