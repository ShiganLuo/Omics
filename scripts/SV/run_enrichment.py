import pandas as pd
import logging
from enricher.function import enrich_go
import os


logger = logging.getLogger(__name__)

def sv_go(
        anno_file: str,
        col_line_prefiex: str = "#Uploaded_variation",
        IMPACT_filter: list = ["HIGH", "MODERATE", "LOW"],
        gene_col: str = "SYMBOL",
        outdir: str = "enrichment_results"
):
    """
    Funtion: 
    """
    os.makedirs(outdir, exist_ok=True)
    header_line = None
    with open(anno_file) as f:
        for i, line in enumerate(f):
            if line.startswith(col_line_prefiex):
                header_line = i
                break

    df = pd.read_csv(anno_file, sep='\t', skiprows=header_line)
    df_filtered = df[df["IMPACT"].isin(IMPACT_filter)]
    logger.info(f"Total SVs: {len(df)}, Filtered SVs (IMPACT in {IMPACT_filter}): {len(df_filtered)}")
    enrich_go(df_filtered[gene_col].dropna().unique().tolist(), outdir=outdir)

if __name__ == "__main__":
    anno_file = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB_annotated.tab"
    outdir = "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/enrichment/go"
    sv_go(anno_file, outdir=outdir,gene_col="SYMBOL")