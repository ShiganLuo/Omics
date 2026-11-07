include: "../common/common.smk"

import os
import time

outdir = config.get("outdir", "output/regionplot")
logdir = config.get("logdir", "log/group")
ROOT_DIR = config.get("ROOT_DIR", ".")
samples = config.get("samples", [])
p = config.get("Params", {}).get("regionplot", {})
genome_cfg = config.get("genome", {})

# tracks: [{label, dir, suffix}] — one track group per entry; per-sample bw
# paths expand to {dir}/{sample_id}/{sample_id}{suffix}. Wired by each
# subworkflow from its own bw-producing rules; not a user-facing key.
_track_defs = config.get("tracks", []) or []
_regions_bed = p.get("regions_bed")
_formats = p.get("formats") or ["png", "svg"]
if isinstance(_formats, str):
    _formats = [_formats]
if not isinstance(_formats, list) or not all(isinstance(f, str) for f in _formats):
    raise ValueError("Params.regionplot.formats must be a list of strings")
_sample_order = p.get("sample_order") or list(samples)
if isinstance(_sample_order, str):
    _sample_order = [_sample_order]
_group_dict = p.get("sample_group_dict", {}) or {}


def _regionplot_inputs(wildcards):
    """Collect every input file; required keys raise at DAG-build time.

    Returns a dict with the explicit bw list, the regions BED, and the
    optional GTF (present only when configured).
    """
    if not _regions_bed:
        raise ValueError(
            "Params.regionplot.regions_bed is empty: set it in the workflow "
            "config JSON (a BED of chrom/start/end[/name] rows), then rerun.")
    files = {
        "bws": [f"{t['dir']}/{s}/{s}{t['suffix']}"
                for t in _track_defs for s in _sample_order],
        "regions": _regions_bed,
    }
    gtf = genome_cfg.get("gtf")
    if gtf and p.get("gtf_track", True):
        files["gtf"] = gtf
    return files


def _panel_names(regions_bed):
    """Panel names = regions BED col4 (deduped) or chrom_start_end."""
    from pathlib import Path as _P
    names = []
    seen = {}
    with open(regions_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            name = (f[3].strip() if len(f) >= 4 and f[3].strip()
                    else f"{f[0]}_{f[1]}_{f[2]}")
            name = name.replace("/", "_").replace(" ", "_")
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
        plots = [outdir + "/{panel}/{panel}." + fmt for fmt in _formats],
        inis = outdir + "/{panel}/tracks.ini",
        tsvs = outdir + "/{panel}/track_sources.tsv",
        bed = outdir + "/{panel}/regions.bed",
    wildcard_constraints:
        panel = r"[^/]+"
    log:
        logdir + "/regionplot/regionplot_{panel}.log"
    threads: rule_res("regionplot").get("threads", 1)
    resources:
        mem_mb = rule_res("regionplot").get("mem_mb", 8000),
    conda:
        "regionplot.yaml"
    container:
        sif("regionplot.yaml", "regionplot")
    params:
        script = ROOT_DIR + "/modules/regionplot/bin/regionplot.py",
        track_defs = _track_defs,
        sample_order = _sample_order,
        group_dict = _group_dict,
        formats = _formats,
        gtf_track = p.get("gtf_track", True),
        track_height = p.get("plot_height_per_track", 1.5),
        width = p.get("plot_width", 12),
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("regionplot", log_file=log_path)
        try:
            import json as _json
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start regionplot for panel {wildcards.panel} at {current_time}")
            panel_outdir = os.path.dirname(str(output.plots[0]))
            os.makedirs(panel_outdir, exist_ok=True)
            script = os.path.join(RECEIPT_DIR, f"regionplot_{wildcards.panel}_{current_time}.sh")
            plain_script = os.path.join(panel_outdir, "regionplot.sh")

            # regions BED copy delivered beside the plots (traceability)
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
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during regionplot for {wildcards.panel}: {e}\n")
            raise RuntimeError(f"Error occurred during regionplot for {wildcards.panel}: {e}\n")


_panels = _panel_names(_regions_bed) if _regions_bed else []

rule regionplot_result:
    input:
        plots = expand(outdir + "/{panel}/{panel}.{fmt}",
                       panel=_panels, fmt=_formats),
        tsvs = expand(outdir + "/{panel}/track_sources.tsv", panel=_panels),
