import pandas as pd
import logging
from enricher.function import enrich_go, enrich_kegg
from utils.VEP_SV import read_vep_tab
import os
import logging
import subprocess
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
StreamHandler = logging.StreamHandler()
StreamHandler.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(StreamHandler)

def sv_go(
        anno_file: str,
        IMPACT_filter: list = ["HIGH", "MODERATE", "LOW"],
        gene_col: str = "SYMBOL",
        outdir: str = "enrichment_results",
        **kwargs
):
    """
    Funtion: 
    """
    os.makedirs(outdir, exist_ok=True)
    df = read_vep_tab(anno_file, **kwargs)
    df_filtered = df[df["IMPACT"].isin(IMPACT_filter)]
    logger.info(f"Total SVs: {len(df)}, Filtered SVs (IMPACT in {IMPACT_filter}): {len(df_filtered)}")
    try:
        enrich_go(df_filtered[gene_col].dropna().unique().tolist(), outdir=f"{outdir}/go")
    except Exception as e:
        logger.error(f"Error during enrichment: {e}")
    try:
        enrich_kegg(df_filtered[gene_col].dropna().unique().tolist(), outdir=f"{outdir}/kegg")
    except Exception as e:
        logger.error(f"Error during enrichment: {e}")

def hot_spot(
    annotated_tab:str,
    outfile:str,
    fai_file:str,
    window_size:int = 1000000,
    threshold:int = 10,
    **kwargs
):
    """
    Function: Identify genomic hotspots of SVs by counting the number of variants in fixed-size windows across the genome, and output significant hotspots.
    Parameters:
        - annotated_tab (str): Path to the annotated SV file in tab-delimited format.
        - fai_file (str): Path to the reference genome index file (.fai) containing chromosome sizes.
        - window_size (int): Size of the genomic window in base pairs (default: 1Mb).
        - threshold (int): Minimum number of variants in a window to be considered a hotspot (default: 10).
        - out_file (str): Path to the output file for significant hotspots in BED format (default: "hotspots.bed").
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

def annotate_hotspots_with_bedtools(hotspot_file: str, gtf_file: str, output_file: str):
    """
    Annotate genomic hotspots using BEDTools intersect with a GTF file.

    Parameters:
        hotspot_file (str): Path to the hotspot file (in BED format).
        gtf_file (str): Path to the GTF file for annotation.
        output_file (str): Path to save the annotated output.
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

<<<<<<< HEAD

=======
def tab_parser(
        table_file: str,
        **kwargs
):
    df = read_vep_tab(table_file, **kwargs)
    print(df["miRNA"].value_counts())
>>>>>>> 0ea0c995979175199379f89ae752af4c876178c6
    

def main():
    # anno_file = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_annotated.tab"
    # outdir = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/enrichment"

    # sv_go(anno_file, outdir=outdir,gene_col="SYMBOL")

    # hot_spot(
    #     annotated_tab="/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_annotated.tab",
    #     outfile="/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/hotspots/hotspots_1000.tsv",
    #     fai_file="/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai"
    # )
    # annotate_hotspots_with_bedtools(
    #     hotspot_file="/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/hotspots/hotspots_1000.tsv",
    #     gtf_file="/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf",
    #     output_file="/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/hotspots/annotated_hotspots.tsv"
    # )
<<<<<<< HEAD

=======
    tab_parser("/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_annotated.tab")
>>>>>>> 0ea0c995979175199379f89ae752af4c876178c6
    pass

if __name__ == "__main__":
    main()