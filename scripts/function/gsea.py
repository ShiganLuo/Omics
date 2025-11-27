import gseapy as gp
import pandas as pd

def run_gsva(
        count_matrix:str,
        gmt_pathway:str,
        outdir:str
):
    es = gp.gsva(
    data = count_matrix,
    gene_sets = gmt_pathway,
    outdir = outdir
    )
    return es

if __name__ == "__main__":
    run_gsva(
        count_matrix= "/home/luosg/Data/genomeStability/output/result/matrix/tpm230.tsv",
        gmt_pathway= "/home/luosg/Data/genomeStability/data/final.gmt",
        outdir= "/home/luosg/Data/genomeStability/output/result/gsva230"
    )
