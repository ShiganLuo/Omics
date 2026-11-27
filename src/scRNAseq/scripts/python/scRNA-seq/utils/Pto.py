from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
import datetime
import argparse
h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"
figurepath = "/home/lsg/Data/glioblastoma/output/new/figure/Pto"
def Pto(adata,out,cell_type):
    sc.tl.diffmap(adata)
    root_ixs = adata.obsm["X_diffmap"][:, 3].argmin()
    sc.pl.scatter(
        adata,
        basis="diffmap",
        color=[cell_type],
        components=[2, 3],
    )
    plt.savefig(f"{figurepath}/{out}.png")
    adata.uns["iroot"] = root_ixs
    sc.tl.dpt(adata)
    sc.pl.dpt_groups_pseudotime()

if __name__ == '__main__':
    start = datetime.datetime.now()
    SC = sc.read_h5ad(f"{h5ad}/SC-cnv.h5ad")
    Pto(SC,"SC","major_celltype")
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")