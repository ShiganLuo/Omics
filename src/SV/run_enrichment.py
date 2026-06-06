import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.LogUtil import setup_logger
import pandas as pd
import logging
from enricher.function import enrich_go, enrich_kegg
from utils.VEP_SV import read_vep_tab
import os
import logging
import subprocess
import argparse
from typing import List, Literal, Optional
PlotFormat = Literal["png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg", "pgf", "raw", "rgba"]
logger = setup_logger("SVEnrichment", level=logging.INFO)

def sv_go(
        anno_file: str,
        IMPACT_filter: List[str] = None,
        gene_col: str = "SYMBOL",
        outdir: str = "enrichment_results",
        image_formats: Optional[List[PlotFormat]] = None,
        **kwargs
) -> None:
    """Perform GO and KEGG enrichment analysis on structural variants.

    Reads a VEP-annotated tab file, filters variants by IMPACT level,
    extracts unique gene symbols, and runs GO and KEGG enrichment.

    Parameters
    ----------
    anno_file : str
        Path to the VEP-annotated SV file (tab-delimited).
    IMPACT_filter : list of str, optional
        IMPACT levels to include. Default is
        ``["HIGH", "MODERATE", "LOW"]``.
    gene_col : str, optional
        Column name containing gene symbols. Default is ``"SYMBOL"``.
    outdir : str, optional
        Output directory for enrichment results. Default is
        ``"enrichment_results"``.
    image_formats : list of str, optional
        Output image formats (e.g. ``["png", "pdf"]``). Default is
        ``["png"]``.
    **kwargs
        Additional keyword arguments passed to ``read_vep_tab``.

    Returns
    -------
    None
    """
    if IMPACT_filter is None:
        IMPACT_filter = ["HIGH", "MODERATE", "LOW"]
    if image_formats is None:
        image_formats = ["png"]
    os.makedirs(outdir, exist_ok=True)
    df = read_vep_tab(anno_file, **kwargs)
    df_filtered = df[df["IMPACT"].isin(IMPACT_filter)]
    logger.info(f"Total SVs: {len(df)}, Filtered SVs (IMPACT in {IMPACT_filter}): {len(df_filtered)}")
    try:
        enrich_go(df_filtered[gene_col].dropna().unique().tolist(), outdir=f"{outdir}/go", image_formats=image_formats)
    except Exception as e:
        logger.error(f"Error during enrichment: {e}")
    try:
        enrich_kegg(df_filtered[gene_col].dropna().unique().tolist(), outdir=f"{outdir}/kegg", image_formats=image_formats)
    except Exception as e:
        logger.error(f"Error during enrichment: {e}")

def hot_spot(
    annotated_tab: str,
    outfile: str,
    fai_file: str,
    window_size: int = 1000000,
    threshold: int = 10,
    **kwargs
) -> None:
    """Identify genomic hotspots of SVs using a sliding-window approach.

    Counts the number of SVs falling in fixed-size windows across the
    genome and outputs windows exceeding a significance threshold.

    Parameters
    ----------
    annotated_tab : str
        Path to the annotated SV file in tab-delimited format.
    outfile : str
        Path to the output file for significant hotspots (TSV format).
    fai_file : str
        Path to the reference genome index file (``.fai``) containing
        chromosome sizes.
    window_size : int, optional
        Size of the genomic window in base pairs. Default is 1000000
        (1 Mb).
    threshold : int, optional
        Minimum number of variants in a window to be considered a
        hotspot. Default is 10.
    **kwargs
        Additional keyword arguments passed to ``read_vep_tab``.

    Returns
    -------
    None
    """
    dir = os.path.dirname(outfile)
    os.makedirs(dir, exist_ok=True)
    df = read_vep_tab(annotated_tab, **kwargs)
    df["chrom"] = df["Location"].str.split(":").str[0]
    df["start"] = df["Location"].str.split(":").str[1].str.split("-").str[0].astype(int)
    df["end"] = df["Location"].str.split(":").str[1].str.split("-").str[1].astype(int)
    df = df[["chrom", "start", "end"]]
    df_fai = pd.read_csv(fai_file, sep="\t", header=None, names=["chrom", "size", "offset", "line_bases", "line_width"])
    hotspots = []
    for chrom, size in zip(df_fai["chrom"], df_fai["size"]):
        for start in range(0, size, window_size):
            end = min(start + window_size, size)
            count = df[(df["chrom"] == chrom) & (df["start"] >= start) & (df["start"] < end)].shape[0]
            hotspots.append([chrom, start, end, count])

    hotspots_df = pd.DataFrame(hotspots, columns=["chrom", "start", "end", "count"])
    significant_hotspots = hotspots_df[hotspots_df["count"] > threshold]
    significant_hotspots.to_csv(outfile, sep="\t", index=False)

def annotate_hotspots_with_bedtools(hotspot_file: str, gtf_file: str, output_file: str) -> None:
    """Annotate genomic hotspots by intersecting with a GTF file via BEDTools.

    Runs ``bedtools intersect -wa -wb`` to overlap hotspot regions with
    gene annotations from a GTF file.

    Parameters
    ----------
    hotspot_file : str
        Path to the hotspot file in BED or BED-like format.
    gtf_file : str
        Path to the GTF annotation file.
    output_file : str
        Path to save the annotated intersection output.

    Returns
    -------
    None
    """
    try:
        # Construct the BEDTools intersect command
        command = [
            "bedtools", "intersect",
            "-a", hotspot_file,
            "-b", gtf_file,
            "-wa", "-wb"
        ]

        # Run the command and capture the output
        with open(output_file, "w") as out:
            subprocess.run(command, stdout=out, check=True)

        print(f"Annotation completed. Results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running BEDTools: {e}")
    except FileNotFoundError:
        print("BEDTools is not installed or not found in PATH.")

def tab_parser(
        table_file: str,
        **kwargs
) -> None:
    """Parse a VEP tab file and print miRNA feature value counts.

    Parameters
    ----------
    table_file : str
        Path to the VEP-annotated tab-delimited file.
    **kwargs
        Additional keyword arguments passed to ``read_vep_tab``.

    Returns
    -------
    None
    """
    df = read_vep_tab(table_file, **kwargs)
    print(df["miRNA"].value_counts())
    

def parse_args():
    """Parse command-line arguments for SV enrichment analysis.

    Returns
    -------
    argparse.Namespace
        Parsed arguments including ``anno_file``, ``outdir``,
        ``gene_col``, ``IMPACT_filter``, and ``formats``.
    """
    parser = argparse.ArgumentParser(description="Run SV enrichment analysis")
    parser.add_argument("-a", "--anno_file", required=True, help="Path to vep annotated SV file (tab-delimited)")
    parser.add_argument("-g","--gtf", help="Path to GTF file for hotspot annotation with bedtools")
    parser.add_argument("-fi", "--fai", help="Path to reference genome index file (.fai) for hotspot analysis")
    parser.add_argument("-o", "--outdir", default="enrichment_results", help="Output directory for enrichment results")
    parser.add_argument("--gene_col", default="SYMBOL", help="Column name for gene symbols in the annotation file")
    parser.add_argument("--IMPACT_filter", nargs="+", default=["HIGH", "MODERATE", "LOW"], help="IMPACT levels to include in enrichment")
    parser.add_argument("-f", "--format", action="append", dest="formats", metavar="FMT", help="Image output format (png, pdf, svg, ...). Can be specified multiple times. Default: png.")
    return parser.parse_args()
def main():
    """Entry point for the SV enrichment analysis pipeline.

    Parses command-line arguments and runs GO/KEGG enrichment and
    tab file parsing.

    Returns
    -------
    None
    """
    args = parse_args()
    anno_file = args.anno_file
    outdir = args.outdir
    fai_file = args.fai
    gtf_file = args.gtf
    sv_go(anno_file, outdir=outdir, gene_col="SYMBOL", image_formats=args.formats)

    out_hot = os.path.join(outdir, "hotspots_1000.tsv")
    hot_spot(
        annotated_tab=anno_file,
        outfile=out_hot,
        fai_file=fai_file
    )
    out_annotated_hot = os.path.join(outdir, "annotated_hotspots.tsv")
    annotate_hotspots_with_bedtools(
        hotspot_file=out_hot,
        gtf_file=gtf_file,
        output_file=out_annotated_hot
    )

if __name__ == "__main__":
    main()