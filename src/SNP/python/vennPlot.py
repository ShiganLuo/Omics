from typing import List
import pandas as pd
import sys
import matplotlib.pyplot as plt
sys.path.append(r'pyvenn')
from functools import reduce
import venn
from matplotlib_venn import venn3

def plotVenn(fileList:List,col:str,outfile:str,outIntersect:str):
    vennList = []
    for file in fileList:
        df = pd.read_csv(file,sep="\t")
        purList = df[col].to_list()
        vennList.append(purList)


    venn3([set(g) for g in vennList], set_labels=("2C_GFPp_SA1","TPS","ciTotiSC"),set_colors=("red", "green", "blue"),alpha=0.4)

    plt.savefig(outfile)
    plt.close()
    intersection = reduce(lambda x, y: set(x) & set(y), vennList)
    dfIntersect = df[df[col].isin(intersection)]
    dfIntersect.to_csv(outIntersect,sep="\t",index=False)
    return intersection
if __name__ == '__main__':
    files1 = ["/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/SNP/go/up_GO.csv",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE183522/SNP/go/up_GO.csv",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE185005/SNP/go/up_GO.csv"]
    plotVenn(files1,"ID","/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upGo_venn.png",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upGo_venn.csv")

    files2 = ["/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/SNP/go/down_GO.csv",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE183522/SNP/go/down_GO.csv",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE185005/SNP/go/down_GO.csv"]
    plotVenn(files2,"ID","/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/downGo_venn.png",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/downGo_venn.csv")

    files3 = ["/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/SNP/kegg/up_kegg.csv",
            "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE183522/SNP/kegg/up_kegg.csv",
            "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE185005/SNP/kegg/up_kegg.csv"]
    plotVenn(files3,"ID","/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upKegg_venn.png",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/upKegg_venn.csv")
    
    files4 = ["/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE166216/SNP/kegg/down_kegg.csv",
        "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE183522/SNP/kegg/down_kegg.csv",
        "/ChIP_seq_2/StemCells/RNASNP_PT/output/GSE185005/SNP/kegg/down_kegg.csv"]
    
    plotVenn(files4,"ID","/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/downKegg_venn.png",
             "/ChIP_seq_2/StemCells/RNASNP_PT/output/plot/mouse/SNP/downKegg.csv")

