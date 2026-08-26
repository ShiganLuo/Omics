import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import datetime
import argparse
gtf="/home/lsg/Data/glioblastoma/data/ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
figurepath = "/home/lsg/Data/glioblastoma/output/new/figure/cnv"
h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"
def infercnv(adata,out,cell_type):
    cnv.io.genomic_position_from_gtf(gtf,adata)
    print(adata.var.loc[:, ["chromosome", "start", "end"]].head())
    sc.pl.umap(adata, color=cell_type)
    plt.savefig(f"{figurepath}/{out}.png")
    # We provide all immune cell types as "normal cells".
    cnv.tl.infercnv(
        adata,
        reference_key=cell_type,
        reference_cat=[
            "Macrophages",
            "Mast Cell",
        ],
        window_size=250,
    )
    cnv.pl.chromosome_heatmap(adata, groupby=cell_type)
    plt.savefig(f"{figurepath}/{out}-heatmap.png")
    ## leiden
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata)
    cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True)
    plt.savefig(f"{figurepath}/{out}-heatmap-leiden.png")
    ## umap based on cnv cluster
    cnv.tl.umap(adata)
    cnv.tl.cnv_score(adata)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
    ax4.axis("off")
    cnv.pl.umap(
        adata,
        color="cnv_leiden",
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,
    )
    cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(adata, color=cell_type, ax=ax3)
    fig.savefig(f"{figurepath}/{out}-umap-cnv.png")
    ## umap based on transcriptomics data
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
    # ax4.axis("off")
    # sc.pl.umap(adata, color="cnv_leiden", ax=ax1, show=False)
    # sc.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
    # sc.pl.umap(adata, color=cell_type, ax=ax3)
    # fig.savefig(f"{figurepath}/{out}-umap-transcriptomics.png")
    ## clasify tumor cells
    adata.obs["cnv_status"] = "normal"
    adata.obs.loc[adata.obs["cnv_leiden"].isin(["0", "1", "2", "3", "4", "5", "6", "7", "10"]), "cnv_status"] = (
        "tumor"
    )
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
    cnv.pl.umap(adata, color="cnv_status", ax=ax1, show=False)
    sc.pl.umap(adata, color="cnv_status", ax=ax2, show=False)
    sc.pl.umap(adata, color=cell_type, ax=ax3, show=False)
    sc.pl.umap(adata, color="cnv_score", ax=ax4)
    fig.savefig(f"{figurepath}/{out}-cell_status.png")
    cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "tumor", :])
    plt.savefig(f"{figurepath}/{out}-heatmap-tumor.png")
    cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "normal", :])
    plt.savefig(f"{figurepath}/{out}-heatmap-normal.png")  
    adata.write_h5ad(f"{h5ad}/{out}-cnv.h5ad")
if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    SC = sc.read_h5ad(f"{h5ad}/SC-ann.h5ad")
    # TE = sc.read_h5ad(f"{h5ad}/TE-bbknn.h5ad")
    # GBM = sc.read_h5ad(f"{h5ad}/GBM-0.3-bbknn.h5ad")
    # ad_dict = {"SC":SC,"TE":TE,"GBM":GBM}
    ad_dict = {"SC":SC}
    for out,adata in ad_dict.items():
        infercnv(adata,out,"major_celltype")
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")