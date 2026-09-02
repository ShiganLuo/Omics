#!/usr/bin/env python3
"""Run computeMatrix + plotHeatmap in one step.

Supports three region modes:
  1. Pre-made BED: --regions regions.bed
  2. TSS from GTF:  --gtf genes.gtf --region-mode tss --tss-flank 1000
  3. Gene list:     --gtf genes.gtf --region-mode genes --gene-names TP53 BRCA1 --merge

Usage:
    # Pre-made BED
    python run_heatmap.py \\
        --ip-bigwig IP.bigwig --regions regions.bed --output heatmap.png

    # TSS mode
    python run_heatmap.py \\
        --ip-bigwig IP.bigwig --gtf genes.gtf --region-mode tss \\
        --tss-flank 1000 --output heatmap.png

    # Gene list mode
    python run_heatmap.py \\
        --ip-bigwig IP.bigwig --gtf genes.gtf --region-mode genes \\
        --gene-names TP53 BRCA1 --merge --output heatmap.png
"""
import argparse
import logging
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("run_heatmap")


# ============================================================
# GTF attribute parsing (shared by TSS and gene extraction)
# ============================================================

def _parse_gtf_attribute(attr_str, key):
    """Extract a value from GTF attributes string."""
    for part in attr_str.rstrip(";").split(";"):
        part = part.strip()
        if part.startswith(key + " "):
            return part.split('"')[1] if '"' in part else part.split()[-1]
    return None


# ============================================================
# BED region generation from GTF
# ============================================================

def generate_tss_bed(gtf_path, output_path, flank=1000, feature="transcript"):
    """Generate TSS BED from GTF. Returns output_path."""
    seen = set()
    count = 0
    with open(gtf_path) as gtf, open(output_path, "w") as out:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != feature:
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            gene_name = _parse_gtf_attribute(fields[8], "gene_name") or _parse_gtf_attribute(fields[8], "gene_id")
            if not gene_name:
                continue
            tss = start if strand == "+" else end
            bed_start, bed_end = max(0, tss - flank), tss + flank
            key = (chrom, bed_start, bed_end, gene_name)
            if key in seen:
                continue
            seen.add(key)
            out.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_name}\t0\t{strand}\n")
            count += 1
    logger.info(f"Generated {count} TSS regions (+-{flank}bp) from {gtf_path}")
    return output_path


def generate_genes_bed(gtf_path, output_path, gene_names, match_by="gene_name",
                       feature="exon", merge=False):
    """Extract gene/TE regions BED from GTF. Returns output_path."""
    target_names = set(gene_names)
    if not target_names:
        logger.error("No gene/TE names specified.")
        return output_path

    name_features = defaultdict(list)
    with open(gtf_path) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != feature:
                continue
            name = _parse_gtf_attribute(fields[8], match_by)
            if name and name in target_names:
                name_features[name].append((fields[0], int(fields[3]), int(fields[4]), fields[6]))

    if not name_features:
        logger.warning(f"No features found matching {match_by} in {target_names}")
        return output_path

    count = 0
    with open(output_path, "w") as out:
        for name in sorted(name_features):
            features = name_features[name]
            if merge:
                chrom = features[0][0]
                strand = features[0][3]
                min_start = min(f[1] for f in features)
                max_end = max(f[2] for f in features)
                out.write(f"{chrom}\t{min_start}\t{max_end}\t{name}\t0\t{strand}\n")
                count += 1
            else:
                for chrom, start, end, strand in features:
                    out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
                    count += 1

    found = sorted(name_features)
    missing = target_names - set(found)
    logger.info(f"Extracted {count} BED regions for {len(found)} names ({match_by})")
    if missing:
        logger.warning(f"Not found: {', '.join(sorted(missing))}")
    return output_path


# ============================================================
# Subprocess helpers
# ============================================================

def run_cmd(cmd):
    logger.info(f"Running: {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        logger.error(r.stderr)
        sys.exit(1)


# ============================================================
# Main
# ============================================================

def main():
    p = argparse.ArgumentParser("Run computeMatrix + plotHeatmap")
    # Input
    p.add_argument("--ip-bigwig", required=True)
    p.add_argument("--input-bigwig", default="")
    p.add_argument("--regions", default="", help="Pre-made BED file (mutually exclusive with --gtf)")
    p.add_argument("--output", required=True, help="Output heatmap PNG")
    # GTF-based region generation
    p.add_argument("--gtf", default="", help="GTF annotation file (generates regions on the fly)")
    p.add_argument("--region-mode", default="", choices=["tss", "genes", ""],
                   help="Region generation mode: tss or genes")
    p.add_argument("--tss-flank", type=int, default=1000, help="Flanking region around TSS (bp)")
    p.add_argument("--gene-names", action="append", default=[], help="Gene/TE names to extract (repeatable)")
    p.add_argument("--match-by", default="gene_name",
                   choices=["gene_name", "gene_id", "family_id", "class_id"],
                   help="GTF attribute to match gene names")
    p.add_argument("--merge", action="store_true", help="Merge features into single region per gene")
    # computeMatrix params
    p.add_argument("--computeMatrix", default="computeMatrix")
    p.add_argument("--mode", default="reference-point")
    p.add_argument("--reference-point", default="center")
    p.add_argument("--before", type=int, default=3000)
    p.add_argument("--after", type=int, default=3000)
    p.add_argument("--upstream", type=int, default=3000)
    p.add_argument("--downstream", type=int, default=3000)
    p.add_argument("--body-length", type=int, default=5000)
    p.add_argument("--bin-size", type=int, default=10)
    p.add_argument("--sort-using", default="mean")
    p.add_argument("--top-n", type=int, default=0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--no-missing-as-zero", action="store_true")
    # plotHeatmap params
    p.add_argument("--plotHeatmap", default="plotHeatmap")
    p.add_argument("--title", default="")
    p.add_argument("--color-map", default="YlOrRd")
    p.add_argument("--height", type=int, default=15)
    p.add_argument("--width", type=int, default=8)
    p.add_argument("--what-to-show", default="plot, heatmap and colorbar")
    p.add_argument("--dpi", type=int, default=300)
    # Keep intermediate files
    p.add_argument("--keep-matrix", default="", help="If set, save matrix.gz to this path")
    p.add_argument("--keep-regions-bed", default="", help="If set, save generated BED to this path")
    args = p.parse_args()

    # --- Resolve regions BED ---
    regions_file = args.regions
    generated_bed = None

    if args.gtf and args.region_mode:
        # Generate BED from GTF on the fly
        if args.keep_regions_bed:
            bed_path = args.keep_regions_bed
        else:
            fd, bed_path = tempfile.mkstemp(suffix=".bed", prefix="regions_")
            os.close(fd)
            generated_bed = bed_path

        if args.region_mode == "tss":
            generate_tss_bed(args.gtf, bed_path, flank=args.tss_flank)
        elif args.region_mode == "genes":
            generate_genes_bed(args.gtf, bed_path, args.gene_names,
                               match_by=args.match_by, merge=args.merge)
        regions_file = bed_path
    elif not regions_file:
        logger.error("Either --regions or --gtf + --region-mode is required.")
        sys.exit(1)

    # --- Select top N peaks ---
    top_bed = None
    if args.top_n > 0:
        logger.info(f"Selecting top {args.top_n} peaks by score")
        top_bed = tempfile.mktemp(suffix=".bed")
        subprocess.run(f"sort -k5 -nr {regions_file} | head -n {args.top_n} > {top_bed}",
                       shell=True, check=True)
        regions_file = top_bed

    # --- computeMatrix ---
    # Include gene name or mode in matrix path to avoid overwrites when running in parallel
    matrix_label = args.gene_names[0] if args.gene_names else args.region_mode or "regions"
    matrix_path = args.keep_matrix or tempfile.mktemp(suffix=f"_{matrix_label}_matrix.gz")
    signals = [args.ip_bigwig]
    if args.input_bigwig:
        signals.append(args.input_bigwig)

    cm_cmd = [
        args.computeMatrix, args.mode,
        "--binSize", str(args.bin_size),
        "--sortUsing", args.sort_using,
        "--numberOfProcessors", str(args.threads),
        "-R", regions_file,
    ]
    for s in signals:
        cm_cmd += ["-S", s]
    cm_cmd += ["-o", matrix_path]

    if args.mode == "reference-point":
        cm_cmd += ["--referencePoint", args.reference_point,
                    "--beforeRegionStartLength", str(args.before),
                    "--afterRegionStartLength", str(args.after)]
    elif args.mode == "scale-regions":
        cm_cmd += ["--regionBodyLength", str(args.body_length),
                    "--upstream", str(args.upstream),
                    "--downstream", str(args.downstream)]

    if not args.no_missing_as_zero:
        cm_cmd.append("--missingDataAsZero")

    run_cmd(cm_cmd)
    logger.info(f"Matrix: {matrix_path}")

    # --- plotHeatmap ---
    ph_cmd = [
        args.plotHeatmap,
        "-m", matrix_path,
        "-o", args.output,
        "--colorMap", args.color_map,
        "--heatmapHeight", str(args.height),
        "--heatmapWidth", str(args.width),
        "--whatToShow", args.what_to_show,
        "--dpi", str(args.dpi),
    ]
    if args.title:
        ph_cmd += ["--plotTitle", args.title]
    # Label: use sample name from ratio bigwig filename
    sample_label = os.path.splitext(os.path.basename(args.ip_bigwig))[0].replace("_IP_over_Input", "")
    ph_cmd += ["--samplesLabel", sample_label]

    run_cmd(ph_cmd)
    logger.info(f"Heatmap: {args.output}")

    # --- Cleanup ---
    if not args.keep_matrix and os.path.exists(matrix_path):
        os.unlink(matrix_path)
    if top_bed and os.path.exists(top_bed):
        os.unlink(top_bed)
    if generated_bed and os.path.exists(generated_bed):
        os.unlink(generated_bed)


if __name__ == "__main__":
    main()
