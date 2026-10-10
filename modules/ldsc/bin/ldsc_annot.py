#!/usr/bin/env python3
"""Prepare LDSC partitioned-heritability annotations from FIRE peak clusters.

Converts a peak BED (optionally with a cluster/class column in the name field,
e.g. "cluster1") into per-cluster annotation BEDs and a summary manifest that
the ldsc module rules consume with make_annot.py / ldsc.py --h2-cts [1].

Usage:
    python ldsc_annot.py -p fire_peaks.bed --chrom-sizes genome.chrom.sizes \
        -o results/ldsc/annot --prefix sample
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence

_SRC_DIR = Path(__file__).resolve().parent.parent.parent.parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
from src.common.util.LogUtil import setup_logger  # noqa: E402


def _get_logger(log_file: Optional[Path] = None) -> logging.Logger:
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        lg = setup_logger("LdscAnnot", log_file=str(log_file))
    else:
        lg = setup_logger("LdscAnnot")
    for handler in lg.handlers:
        if type(handler) is logging.StreamHandler:
            handler.setStream(sys.stdout)
    return lg


logger = _get_logger()


def run(args: argparse.Namespace) -> int:
    clusters: Dict[str, List[str]] = defaultdict(list)
    all_rows: List[str] = []
    with open(args.peaks) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[0], f[1], f[2]
            name = f[3] if len(f) > 3 else "all"
            cluster = name.split("|")[-1] if args.cluster_from_name else "all"
            clusters[cluster].append(f"{chrom}\t{start}\t{end}\t{name}\n")
            all_rows.append(f"{chrom}\t{start}\t{end}\t{name}\n")

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for cluster, rows in sorted(clusters.items()):
        out = args.outdir / f"{args.prefix}.{cluster}.bed"
        tmp = out.with_suffix(".bed.tmp")
        with open(tmp, "w") as fh:
            fh.writelines(rows)
        tmp.replace(out)
        manifest_rows.append(f"{cluster}\t{out}\t{len(rows)}")
        logger.info(f"Cluster {cluster}: {len(rows)} peaks -> {out}")

    all_out = args.outdir / f"{args.prefix}.all.bed"
    tmp = all_out.with_suffix(".bed.tmp")
    with open(tmp, "w") as fh:
        fh.writelines(all_rows)
    tmp.replace(all_out)
    manifest_rows.append(f"all\t{all_out}\t{len(all_rows)}")

    manifest = args.outdir / f"{args.prefix}.annot_manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("cluster\tbed\tpeaks\n")
        fh.write("\n".join(manifest_rows) + "\n")
    logger.info(f"Wrote manifest {manifest}")
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="ldsc_annot.py",
        description="Prepare LDSC annotations from FIRE peak clusters.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-p", "--peaks", type=Path, required=True,
                   help="peaks BED (name may encode cluster)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--prefix", type=str, default="sample", help="output prefix")
    p.add_argument("-c", "--cluster-from-name", action="store_true",
                   help="derive cluster from the BED name field")
    p.add_argument("-l", "--log", type=Path, default=None, help="log file")
    p.add_argument("-n", "--dry-run", action="store_true", help="only log the plan")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    global logger
    args = parse_args(argv)
    logger = _get_logger(args.log)
    logger.info(f"Command: {' '.join(sys.argv)}")
    if not args.peaks.exists():
        logger.error(f"input not found: {args.peaks}")
        return 1
    if args.dry_run:
        logger.info(f"[DRY-RUN] would split {args.peaks} into LDSC annotations -> {args.outdir}")
        return 0
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
