include: "../common/common.smk"


outdir = config.get("outdir", "output/regionplot")
logdir = config.get("logdir", "log/group")
ROOT_DIR = config.get("ROOT_DIR", ".")
samples = config.get("samples", [])
genome_cfg = config.get("genome", {})

# tracks: [{label, dir, suffix}] — one track group per entry; per-sample bw
# paths expand to {dir}/{sample_id}/{sample_id}{suffix}. Wired by each
# subworkflow from its own bw-producing rules; not a user-facing key.
_track_defs = config.get("tracks", []) or []
# regions: BED rows or a [NAME=]chrom:start-end list (auto-detected per row)
# 参数精确到 rule:regions/formats/sample_order/sample_group_dict 挂 regionplot
# 名下(Params.regionplot.regionplot.*,单 plot rule 消费;无 module 级共享回退),
# config 全量路径读取。regionplot_result 聚合 regionplot 输出文件名,
# 面板名与格式展开同一 rule 键。
# 每 rule 一个 parse 期变量
_regions = config.get("Params", {}).get("regionplot", {}).get("regionplot", {}).get("regions", None)
_sample_order = config.get("Params", {}).get("regionplot", {}).get("regionplot", {}).get("sample_order", None) or list(samples)
if isinstance(_sample_order, str):
    _sample_order = [_sample_order]
_group_dict = config.get("Params", {}).get("regionplot", {}).get("regionplot", {}).get("sample_group_dict", {}) or {}
_formats_regionplot = config.get("Params", {}).get("regionplot", {}).get("regionplot", {}).get("formats", None) or ["png", "svg"]
if isinstance(_formats_regionplot, str):
    _formats_regionplot = [_formats_regionplot]
_regions_regionplot_result = config.get("Params", {}).get("regionplot", {}).get("regionplot", {}).get("regions", None)
_formats_regionplot_result = config.get("Params", {}).get("regionplot", {}).get("regionplot", {}).get("formats", None) or ["png", "svg"]
if isinstance(_formats_regionplot_result, str):
    _formats_regionplot_result = [_formats_regionplot_result]


def _regionplot_inputs(wildcards):
    """Collect every input file; required keys raise at DAG-build time.

    Returns a dict with the explicit bw list, the regions file (BED or
    chr:start-end list), and the optional GTF (present only when
    configured).
    """
    if not _regions:
        raise ValueError(
            "Params.regionplot.regionplot.regions is empty: set it in the "
            "workflow config JSON (a BED of chrom/start/end[/name] rows or a "
            "list of [NAME=]chrom:start-end lines), then rerun.")
    files = {
        "bws": [f"{t['dir']}/{s}/{s}{t['suffix']}"
                for t in _track_defs for s in _sample_order],
        "regions": _regions,
    }
    gtf = genome_cfg.get("gtf")
    if gtf and config.get("Params", {}).get("regionplot", {}).get("gtf_track", True):
        files["gtf"] = gtf
    return files


def _panel_names(regions_file):
    """Panel names from a regions file; mirrors bin/regionplot.py parsing.

    BED rows use col4 (deduped) else chrom_start_end; list rows use the
    ``NAME=`` prefix else the sanitized chrom:start-end string.
    """
    import re as _re
    _region_re = _re.compile(r"([^:]+):([\d,]+)-([\d,]+)")

    def _sanitize(name):
        for ch in ("/", " ", ":", "="):
            name = name.replace(ch, "_")
        return name

    names = []
    seen = {}
    with open(regions_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            line = line.rstrip("\n")
            if "\t" in line:
                f = line.split("\t")
                name = (f[3].strip() if len(f) >= 4 and f[3].strip()
                        else f"{f[0]}_{f[1]}_{f[2]}")
            else:
                s = line.strip()
                if "=" in s:
                    name = s.split("=", 1)[0]
                else:
                    m = _region_re.fullmatch(s)
                    name = (f"{m.group(1)}_"
                            f"{int(m.group(2).replace(',', ''))}_"
                            f"{int(m.group(3).replace(',', ''))}")
            name = _sanitize(name)
            if name in seen:
                seen[name] += 1
                name = f"{name}.{seen[name]}"
            else:
                seen[name] = 0
            names.append(name)
    return names


rule regionplot:
    """Multi-sample region track stacking plots (pyGenomeTracks)."""
    input:
        unpack(_regionplot_inputs)
    output:
        plots = [outdir + "/{panel}/{panel}." + fmt for fmt in _formats_regionplot],
        inis = outdir + "/{panel}/tracks.ini",
        tsvs = outdir + "/{panel}/track_sources.tsv",
        bed = outdir + "/{panel}/regions.bed",
    wildcard_constraints:
        panel = r"[^/]+"
    log:
        logdir + "/regionplot/regionplot_{panel}.log"
    threads: rule_res("regionplot").get("threads", 1)
    benchmark: logdir + "/regionplot/regionplot_{panel}.benchmark.tsv",
    resources:
        mem_mb = lambda wildcards, input, attempt: rule_res(
            "regionplot", input_size=input_size_mb(input),
            attempt=attempt).get("mem_mb", 8000),
    conda:
        "regionplot.yaml"
    container:
        sif("regionplot.yaml", "regionplot")
    params:
        script = ROOT_DIR + "/modules/regionplot/bin/regionplot.py",
        track_defs = _track_defs,
        sample_order = _sample_order,
        group_dict = _group_dict,
        formats = _formats_regionplot,
        gtf_track = config.get("Params", {}).get("regionplot", {}).get("gtf_track", True),
        track_height = config.get("Params", {}).get("regionplot", {}).get("plot_height_per_track", 1.5),
        width = config.get("Params", {}).get("regionplot", {}).get("plot_width", 12),
    run:
        log_path = str(log)
        open(log_path, "a").close()
        rule_logger = setup_logger("regionplot", log_file=log_path)
        try:
            import json as _json
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start regionplot for panel {wildcards.panel} at {current_time}")
            panel_outdir = os.path.dirname(str(output.plots[0]))
            os.makedirs(panel_outdir, exist_ok=True)
            script = receipt_script(f"regionplot_{wildcards.panel}_{current_time}.sh")
            os.makedirs(os.path.dirname(script), exist_ok=True)  # sample_id 可能含路径分隔符
            plain_script = os.path.join(panel_outdir, "regionplot.sh")

            # regions file copy delivered beside the plots (traceability)
            import shutil as _shutil
            _shutil.copyfile(str(input.regions), str(output.bed))

            # track table: label<TAB>sample<TAB>bw, group-major then sample order
            track_tsv = os.path.join(panel_outdir, "track_table.tsv")
            with open(track_tsv, "w") as fh:
                fh.write("label\tsample_id\tbw_path\n")
                for t in params.track_defs:
                    for s in params.sample_order:
                        bw = f"{t['dir']}/{s}/{s}{t['suffix']}"
                        fh.write(f"{t['label']}\t{s}\t{bw}\n")

            # group map (only when any sample carries a group)
            group_map_arg = ""
            if params.group_dict:
                gm = os.path.join(panel_outdir, "group_map.tsv")
                with open(gm, "w") as fh:
                    fh.write("sample_id\tgroup\n")
                    for s in params.sample_order:
                        if s in params.group_dict:
                            fh.write(f"{s}\t{params.group_dict[s]}\n")
                group_map_arg = gm

            cmd = [
                "python", params.script,
                "-r", str(output.bed),
                "--panel", str(wildcards.panel),
                "-t", track_tsv,
                "-o", panel_outdir,
                "--formats", *params.formats,
                "--track-height", str(params.track_height),
                "--width", str(params.width),
                "-l", log_path,
            ]
            if group_map_arg:
                cmd += ["-g", group_map_arg]
            if params.gtf_track and getattr(input, "gtf", ""):
                cmd += ["--gtf", str(input.gtf), "--gtf-track"]
            with open(script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("set -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "regionplot for {wildcards.panel} at {current_time} completed successfully"\n')
            _shutil.copyfile(script, plain_script)
            if SCRIPT_ONLY:  # 测试模式只落盘回执脚本,不执行;touch 占位保 DAG 走完(事后清理)
                script_only_touch(output)
            else:
                shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during regionplot for {wildcards.panel}: {e}\n")
            raise RuntimeError(f"Error occurred during regionplot for {wildcards.panel}: {e}\n")


_panels = _panel_names(_regions_regionplot_result) if _regions_regionplot_result else []

rule regionplot_result:
    input:
        plots = expand(outdir + "/{panel}/{panel}.{fmt}",
                       panel=_panels, fmt=_formats_regionplot_result),
        tsvs = expand(outdir + "/{panel}/track_sources.tsv", panel=_panels),
