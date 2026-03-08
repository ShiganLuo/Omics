from utils.VEP_SV import read_vep_tab
from utils.VEP_SV import VEP_SV
import pandas as pd
import os
def tab_parser(
        table_file: str,
        cosmic_file:str,
        tab_gene_col: str = "SYMBOL",
        cosmic_gene_col:str = "GENE_SYMBOL",
        sv_type_pattern:str = r"(DEL|DUP|INV|INS|TRA)",
        **kwargs
):
    """
    Function: Parse a VEP annotation table and filter for cancer genes.
    Parameters:
        - table_file: path to the VEP annotation table (tab-delimited format)
        - cosmic_file: path to the COSMIC cancer gene list (tab-delimited format)
        - tab_gene_col: column name in the VEP table that contains gene symbols (default "SYMBOL")
        - cosmic_gene_col: column name in the COSMIC file that contains gene symbols (default "GENE_SYMBOL")
    Returns:    
        - DataFrame containing only the rows from the VEP table where the gene is in the COSMIC cancer gene list.
    """
    df_tab = read_vep_tab(table_file, **kwargs)
    
    df_coismic = pd.read_csv(cosmic_file, sep="\t")

    cancer_genes = set(df_coismic[cosmic_gene_col].unique())
    cancer_genes = set(g.upper() for g in cancer_genes if isinstance(g, str))

    df_tab["is_cancer_gene"] = df_tab[tab_gene_col].apply(lambda x: x.upper() in cancer_genes)
    df_tab = df_tab[df_tab["is_cancer_gene"] == True]

    col_line_prefix = kwargs.get("col_line_prefix", "#Uploaded_variation")
    df_tab["svtype"] = df_tab[col_line_prefix].str.extract(sv_type_pattern).fillna("OTHER")

    return df_tab[[tab_gene_col, "svtype"]]

def prepare_oncoprint_data(
    tab_files:dict[str],
    cosmic_file:str,
    gene_col:str = "SYMBOL",
    sample_order:list[str] = None,
    **kwargs
):
    """
    Function: Prepare data for OncoPrint visualization by parsing multiple VEP annotation tables and creating a binary matrix of gene alterations across samples.
    """
    df_list = []
    for sample, tab_file in tab_files.items():
        df_sample = tab_parser(tab_file, cosmic_file=cosmic_file, tab_gene_col=gene_col, **kwargs)
        df_sample["sample"] = sample
        df_list.append(df_sample)
        df = pd.concat(df_list, ignore_index=True)

    df = (
        df.groupby([gene_col, "sample"])["svtype"]
        .apply(lambda x: ",".join(sorted(set(x))))
        .reset_index()
    )

    # pivot 成 OncoPrint 矩阵
    matrix = df.pivot(index=gene_col, columns="sample", values="svtype")
    if sample_order is not None:
        matrix = matrix.reindex(columns=sample_order)
    return matrix
        

def main():
    # vcfs = {
    #     "PlaB_P6": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/PlaB_P6.sv.vcf.gz",
    #     "DMSO_P6": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB06/unphased/DMSO_P6.sv.vcf.gz",
    #     "PlaB_P20": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/PlaB_P20.sv.vcf",
    #     "DMSO_P20": "/data/pub/zhousha/Totipotent20251031/data/Pacbio/PlaB20/unphased/DMSO_P20.sv.vcf"
    # }
    # for sample, vcf in vcfs.items():
    #     out_tab = f"/data/pub/zhousha/Totipotent20251031/PacBio/SV/sta/{sample}_annotated.tab"
    #     sample_annotation(vcf, out_tab)
    tab_files = {
        "PlaB P6": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB06_vs_DMSO06/PlaB_annotated.tab",
        "PlaB P20": "/data/pub/zhousha/Totipotent20251031/PacBio/SV/PlaB20_vs_DMSO20/PlaB_annotated.tab"
    }
    cosmic_file = "/data/pub/zhousha/Totipotent20251031/data/Cancer/Cosmic_CancerGeneCensusHallmarksOfCancer_v103_GRCh38.tsv"
    oncoprint_matrix = prepare_oncoprint_data(tab_files, cosmic_file=cosmic_file,sample_order=["PlaB P6", "PlaB P20"])
    oncoprint_matrix.to_csv("/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/oncoprint_matrix.csv")
    pass

if __name__ == "__main__":
    main()