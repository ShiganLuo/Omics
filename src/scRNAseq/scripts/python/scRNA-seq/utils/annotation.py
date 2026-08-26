import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import omicverse as ov
import celltypist
import argparse
from celltypist import models
figure_path = "/home/lsg/Data/glioblastoma/output/new/figure/annotation"
h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"
tableout = "/home/lsg/Data/glioblastoma/output/new/table"
marker_genes = {
    "Oligodendrocytes": ["MBP","MOG","MAG","CLDN14","KLK6","EML1"],
    "Astrocytes": ["GFAP","SLC1A2","ACSL6","AGT","AQP4","APOE"],
    "Endothelial Cell":["CD93","VWF","EMCN","EGFL7","FLT1","ID3"],
    "T Cells": ["TRBC2","CD3D","CD3G","CD3E","IL7R","LTB"],
    "T Regs": ["IKZF2","FOXP3","CCR4","IL2RA","CTLA4","ENTPD1"],
    "NK Cells": ["NKG7","KLRF1","KLRD1","GNLY","NCR1","GZMA"],
    "Mast Cells": ["KIT","HSD11B1","TPSAB1","IL1RL1","HDC","SLC29A1"],
    "B Cells": ["PXK","MS4A1","CD19","CD74","CD79A","IGHD"],
    "Microglia": ["ITGAM","CX3CR1","TMEM119","P2RY12","AIF1","CSF1R"],
    "Macrophages": ["CD68","NAAA","MARCH1","JAML","TYROBP","CD163"],
    "Perivasc. Macs.": ["CD163","MRC1"],
    "Neutrophils-MDSCs": ["CSF3R","S100A8","IL1R2","TREM1","CEACAM1","HP"],
    "Dentric Cells": ["ITGAX","ZBTB46","CD86","LAMP3","CD83","CD1A"],
    "Neuron": ["PRPH","DISP2","CNTNAP2","TUBB3","RBFOX3","CHAT"],
    "Neuron progenitors": ["SOX1","	SP9","HES5","ARX","EOMES","NXPH1"]
}
marker_genes = {
    "Oligodendrocytes": ["MBP","MOG","MAG","KLK6"],
    "Astrocytes": ["AGT","AQP4","APOE"],
    "Endothelial Cell":["CD93","VWF","EGFL7","FLT1"],
    "Mast Cells": ["KIT","TPSAB1"],
    "Macrophages": ["CD68","TYROBP","CD163"],
    "Neuron": ["CNTNAP2","NXPH1","HES5"],
}
def Show_Markers(adata,out):
    sc.pl.correlation_matrix(adata, 'leiden', figsize=(5,3.5))
    plt.savefig(f"{figure_path}/{out}-correlation.png",dpi = 300)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.settings.verbosity=2    
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
    markers = sc.get.rank_genes_groups_df(adata, group=None, pval_cutoff=.05, log2fc_min=.25)
    markers.to_csv(f"{tableout}/{out}-all_markers.csv", index=False,header=True)
    
    top5 = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
    top5.to_csv(f"{tableout}/{out}-top5_markers.csv",index=False,header =True)
    fig=plt.figure(figsize=(20,54),dpi=100)
    for i in  top5.columns: 
        plt.subplot(4, 7, int(i)+1) #做一个3*3的图 range（9）从0开始，因需要从1开始，所以i+1
        sc.pl.rank_genes_groups_violin(adata, groups=str(i), n_genes=5,show=False)
        plt.tight_layout()
        plt.axis = 'off' #关闭坐标 让图更美观
    fig.savefig(f"{figure_path}/top5-markers.png")      
    
    adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
    
    sc.pl.rank_genes_groups_matrixplot(adata, n_genes=3, use_raw=False, vmin=-3, vmax=3, cmap='bwr', layer='scaled')
    plt.savefig(f"{figure_path}/{out}-genes-matrix.png",dpi = 300)
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=3, cmap='bwr')
    plt.savefig(f"{figure_path}/{out}-genes-violin.png",dpi = 300)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, values_to_plot='logfoldchanges', min_logfoldchange=3, vmax=7, vmin=-7, cmap='bwr')
    plt.savefig(f"{figure_path}/{out}-genes-dot.png",dpi = 300)
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, use_raw=False, swap_axes=True, show_gene_labels=False,
                                vmin=-3, vmax=3, cmap='bwr')
    plt.savefig(f"{figure_path}/{out}-genes-heatmap.png",dpi = 300)
    return adata
def auto_omi(ad,out):
#    sc.pl.dotplot(ad, marker_genes,groupby="leiden", standard_scale="var")
#    plt.savefig("dot.png")
# marker(SC)
    scsa=ov.single.pySCSA(adata=ad,
                        foldchange=1.5,
                        pvalue=0.01,
                        celltype='normal',
                        target='cellmarker',
                        tissue='All',
                        model_path='./pySCSA_2023_v2_plus.db'                    
        )
    sc.pp.log1p(ad) 
    ad.uns['log1p']['base']=10
    res=scsa.cell_anno(clustertype='leiden',
                cluster='all',rank_rep=True)
    res.head()
    # scsa.cell_anno_print()
    scsa.cell_auto_anno(ad,clustertype='leiden',
                        key='scsa_celltype_cellmarker')
    ov.utils.embedding(ad,
                    basis='X_umap',
                    color=[ "leiden","scsa_celltype_cellmarker"],
                    title=['Cell type'],
                    palette=ov.palette()[1:],
                    show=False,frameon='small',wspace=0.35)
    plt.tight_layout()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)

    plt.savefig(f"{figure_path}/{out}.png",dpi = 300)

def handful_annotate(adata,out):
    print(adata.var.head())
    print(adata.obs.head)
    marker_genes_in_data = dict()
    for ct, markers in marker_genes.items():
        markers_found = list()
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found
    print(marker_genes_in_data)
    sc.tl.dendrogram(adata, groupby="leiden")
    sc.pl.dotplot(
        adata,
        groupby="leiden",
        var_names=marker_genes_in_data,
        dendrogram=True,
        standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
        figsize=(12, 10),
    )
    plt.tight_layout()
    plt.savefig(f"{figure_path}/{out}-dot.png",dpi = 300)
    for ct,gene in marker_genes_in_data.items():
        print(f"{ct}:{gene}")  # print cell subtype name
        if len(gene) > 0:
            sc.pl.umap(
                adata,
                color=marker_genes_in_data[ct],
                vmin=0,
                vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
                sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
                frameon=False,
                cmap="Reds",  
            )
            plt.tight_layout()
            plt.savefig(f"{figure_path}/{out}-{ct}.png")
    cluster2annotation = {
     '13': 'Oligodendrocytes',
     '15': 'Endothelial Cell',
     '14': 'Macrophages',#Microglia 
     '5': 'Astrocytes',
     '16': 'Mast Cell',
     '10': 'Astrocytes',#?
     '4': 'Astrocytes',
     '9': 'Astrocytes',
     '12': 'Astrocytes',#？
     '3': 'Astrocytes',
     '5': 'Astrocytes',
     '0': 'Neuron',#?
     '7': 'Neuron',#?
     '6': 'Neuron',
     '8': 'Neuron',
     '2': 'Astrocytes',
     '1': 'Neuron',#?
     '11': 'Neuron',#?
    }
    adata.obs['major_celltype'] = adata.obs['leiden'].map(cluster2annotation).astype('category')
    ov.utils.embedding(adata,
                basis='X_umap',
                color=[ "leiden","major_celltype"],
                title=['Clusters','Cell types'],
                palette=ov.palette()[:],wspace=0.55,
                show=False,frameon='small',)
    plt.savefig(f"{figure_path}/{out}-handfulann.png",dpi = 300)
    adata.write_h5ad(f"{h5ad}/{out}-ann.h5ad")
def auto_cet(adata,out):
    adata_celltypist = adata.copy()  # make a copy of our adata
    adata_celltypist.X = adata.layers["count"]  # set adata.X to raw counts
    sc.pp.normalize_per_cell(
        adata_celltypist, counts_per_cell_after=10**4
    )  # normalize to 10,000 counts per cell
    sc.pp.log1p(adata_celltypist)  # log-transform
    # make .X dense instead of sparse, for compatibility with celltypist:
    adata_celltypist.X = adata_celltypist.X.toarray()
    models.download_models(
    force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"]
    )
    model_low = models.Model.load(model="Immune_All_Low.pkl")
    model_high = models.Model.load(model="Immune_All_High.pkl")
    predictions_high = celltypist.annotate(
    adata_celltypist, model=model_high, majority_voting=True
    )
    predictions_high_adata = predictions_high.to_adata()
    adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
    adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]
    predictions_low = celltypist.annotate(
    adata_celltypist, model=model_low, majority_voting=True
    )
    predictions_low_adata = predictions_low.to_adata()
    adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "majority_voting"
    ]
    adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
        adata.obs.index, "conf_score"
    ]
    sc.pl.umap(
    adata,
    color=["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"],
    frameon=False,
    sort_order=False,
    wspace=1,
    )
    plt.savefig(f"{figure_path}/{out}-1.jpg")
    sc.pl.umap(
    adata,
    color=["celltypist_cell_label_fine", "celltypist_conf_score_fine"],
    frameon=False,
    sort_order=False,
    wspace=1,
    )
    plt.savefig(f"{figure_path}/{out}-2.jpg")
    sc.pl.dendrogram(adata, groupby="celltypist_cell_label_fine")
    plt.savefig(f"{figure_path}/{out}-3.jpg")
if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    SC = sc.read_h5ad(f"{h5ad}/SC-bbknn.h5ad")
    # TE = sc.read_h5ad(f"{h5ad}/TE-bbknn.h5ad")
    # GBM = sc.read_h5ad(f"{h5ad}/GBM-0.3-bbknn.h5ad")
    # ad_dict = {"SC":SC,"TE":TE,"GBM":GBM}
    ad_dict = {"SC":SC}
    for out,adata in ad_dict.items():
        # auto_omi(adata,f"{out}-auto")
        handful_annotate(adata,out)
        # Show_Markers(adata,out)
        # auto_cet(adata,out)
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")