import pandas as pd
import logging
import sys
from pathlib import Path
from typing import Union
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	stream=sys.stdout,
	datefmt='%Y-%m-%d %H:%M:%S'
)


def geneIDAnnotation(gtf_path: Union[str, Path]) -> pd.DataFrame:
    gtf = pd.read_csv(
        gtf_path,
        sep="\t", 
        comment="#", 
        header=None, 
        names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    )

    gtf['gene_id'] = gtf['attribute'].str.extract(r'gene_id\s*"(.*?)"')
    gtf['gene_name'] = gtf['attribute'].str.extract(r'gene_name\s*"(.*?)"')
    gtf['gene_type'] = gtf['attribute'].str.extract(r'gene_type\s*"(.*?)"')
    gtf_filtered = gtf.dropna(subset=["gene_id"])
    gtf_filtered = gtf_filtered.drop_duplicates(
        subset=["gene_id", "gene_name", "gene_type"], 
        keep="first"
    )

    return gtf_filtered[["gene_id", "gene_name", "gene_type"]]

def renameIndex(
        df_annotation:pd.DataFrame,
        df:pd.DataFrame,
        col:str,
) -> pd.DataFrame:
    """
    rename df index from gene_id to gene_name
    """
    df_map = df_annotation.drop_duplicates(subset=["gene_id"], keep="first")
    gene_id_to_name_map = df_map.set_index('gene_id')['gene_name']
    new_index = df[col].map(gene_id_to_name_map)
    df_new = df.copy()
    df_new.index = new_index
    df_new.index.name = 'gene_name'
    df_new.drop(columns=[col],inplace=True)
    return df_new

if __name__ == '__main__':
    # gtf =  "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/gencode.vM36.primary_assembly.annotation.gtf"
    # outfile = "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/geneIDAnnotation.csv"
    gtf = "/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/gencode.v47.primary_assembly.annotation.gtf"
    outfile = "/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/geneIDAnnotation.csv"
    df = geneIDAnnotation(gtf)
    df.to_csv(outfile, sep="\t", index=False, header=True)
