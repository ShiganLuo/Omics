#!/usr/bin/env python3
"""Run ``tidk search`` on an assembly and publish the tidk telomere deliverables.

One-engine orchestration for the tidk_scan rule: runs tidk search over the
assembly FASTA, publishes the search table as the declared TSV (placeholder
when tidk wrote no table at all), then renders the per-contig telomere figure
via plot_tidk. Declared outputs are written on every non-fatal exit path.

Usage:
    python tidk_scan.py -f asm.bp.p_ctg.fa -s s1 -o out/tidk/s1 \\
        -b s1_tidk_telomeres -t tidk -m TTAGGG -F png
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List

try:
    from .LogUtil import setup_logger
    from .plot_tidk import parse_tidk_tsv, plot_tidk
except ImportError:
    from LogUtil import setup_logger
    from plot_tidk import parse_tidk_tsv, plot_tidk

logger = setup_logger("tidk_scan", stream=sys.stdout)


def run_tidk_search(tidk: str, motif: str, sample: str, outdir: str,
                    fasta: str) -> str:
    """Run ``tidk search`` and return the windows TSV path it writes.

    Args:
        tidk: tidk executable.
        motif: telomeric repeat string to search for.
        sample: sample id (names the search output prefix).
        outdir: directory tidk writes its tables to.
        fasta: assembly FASTA to scan.

    Returns:
        Path of the ``*_telomeric_repeat_windows.tsv`` table tidk writes.
    """
    search_prefix = f"{sample}_tidk_search"
    cmd = [tidk, "search", "--string", motif,
           "--output", search_prefix, "--dir", outdir, fasta]
    logger.info("Executing: " + " ".join(shlex.quote(str(x)) for x in cmd))
    subprocess.run(cmd, check=True)
    return os.path.join(outdir, f"{search_prefix}_telomeric_repeat_windows.tsv")


def publish_search_table(windows_tsv: str, declared_tsv: str) -> None:
    """Publish the tidk search table as the declared TSV deliverable.

    Renames the windows table into place (no intermediate left behind).
    When tidk produced no table at all, writes a placeholder instead so the
    declared output exists on every exit path.

    Args:
        windows_tsv: path tidk search wrote.
        declared_tsv: declared output TSV path.
    """
    if os.path.exists(windows_tsv):
        os.replace(windows_tsv, declared_tsv)
        logger.info(f"Wrote: {declared_tsv}")
    else:
        with open(declared_tsv, "w") as f:
            f.write("# No telomeric repeats found by tidk\n")
        logger.warning(
            f"tidk search output not found: {windows_tsv}; "
            f"wrote placeholder {declared_tsv}")


def render_figure(declared_tsv: str, sample: str, out_base: str,
                  formats: List[str]) -> None:
    """Render the per-contig telomere figure from the published TSV.

    Placeholder content renders as the "no parsable data" figure.

    Args:
        declared_tsv: published tidk table (placeholder content allowed).
        sample: sample id for the figure title.
        out_base: figure output base path (no extension).
        formats: output figure formats.
    """
    contig_col, value_col, rows = parse_tidk_tsv(Path(declared_tsv))
    for path in plot_tidk(contig_col, value_col, rows, sample, out_base, formats):
        logger.info(f"Wrote: {path}")


def main() -> None:
    """Parse arguments and publish the tidk deliverables in order."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fasta", required=True,
                        help="Assembly FASTA to scan")
    parser.add_argument("-s", "--sample", required=True,
                        help="Sample id (names the tidk search prefix)")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory for all deliverables")
    parser.add_argument("-b", "--out_base", dest="out_base", required=True,
                        help="Output base name without extension")
    parser.add_argument("-t", "--tidk", default="tidk",
                        help="tidk executable (default: tidk)")
    parser.add_argument("-m", "--motif", default="TTAGGG",
                        help="Telomeric repeat string (default: TTAGGG)")
    parser.add_argument("-F", "--format", dest="formats", action="append",
                        default=None,
                        help="Figure format, repeatable (default: png)")
    args = parser.parse_args()
    formats = args.formats or ["png"]

    os.makedirs(args.outdir, exist_ok=True)
    out_base = os.path.join(args.outdir, args.out_base)
    declared_tsv = out_base + ".tsv"
    windows_tsv = run_tidk_search(args.tidk, args.motif, args.sample,
                                  args.outdir, args.fasta)
    publish_search_table(windows_tsv, declared_tsv)
    render_figure(declared_tsv, args.sample, out_base, formats)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("%s", exc)
        sys.exit(1)
