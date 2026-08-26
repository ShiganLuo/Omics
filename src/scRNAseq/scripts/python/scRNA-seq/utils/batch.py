import scanpy as sc
import matplotlib.pyplot as plt
import datetime
import argparse
figure_path = "/home/lsg/Data/glioblastoma/output/new/figure/batch"
h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"

def Run_Normalization(adata,n_neighbors,n_pcs,resolution,out,fig=0):
    # count can be assign to X;if you don't filter adata with high_variable genes
    adata.layers['count'] = adata.X.copy()
    # annotation need complete var and X
    # adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, 
                            n_top_genes=3000,
                            flavor='seurat',
                            subset=False, 
                            batch_key='sample')
    
    if fig:
        sc.pl.highly_variable_genes(adata)
        plt.savefig(f"{figure_path}/{out}-high_variable_genes.png")
    # avoid wrong and speed;didn't affect the obs.cpy():avoid warning
    # adata = adata[:, adata.var.highly_variable].copy()
    
    #scale;if there exist some extrem value in X,below should be added;otherwise not
    # sc.pp.scale(adata, max_value=10)
    
    sc.tl.pca(adata,use_highly_variable=True)
    #futer versions
    # sc.tl.pca(adata,mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden',flavor="igraph", 
             n_iterations=2, 
             directed=False)                                                                                                                                                 
    
    if fig:
        figure, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        plt.subplots_adjust(wspace=.5)
        sc.pl.umap(adata, color=['sample'], frameon=False, ax=ax1, show=False)
        sc.pl.umap(adata, color=['leiden'], frameon=False, legend_loc='on data', ax=ax2, show=False)
        figure.savefig(f"{figure_path}/{out}-umap.png")
    return adata
def Run_batchRemove(adata,n_neighbors,n_pcs,resolution,out,fig=0,methods=0):
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    if methods == 0:
        sc.external.pp.bbknn(adata, batch_key='sample')
    elif methods == 1:
        sc.external.pp.harmony_integrate(adata, batch_key='sample')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden',flavor="igraph", 
             n_iterations=2, 
             directed=False)
    if fig:
        figure, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        plt.subplots_adjust(wspace=.5)
        sc.pl.umap(adata, color=['sample'], frameon=False, ax=ax1, show=False)
        sc.pl.umap(adata, color=['leiden'], frameon=False, legend_loc='on data', ax=ax2, show=False)
        figure.savefig(f"{figure_path}/{out}-umap-BBKNN.png")    
        sc.pl.umap(adata, color='leiden', add_outline=True, legend_loc='on data',
                    legend_fontsize=12, legend_fontoutline=2,frameon=False,
                    title='clustering of cells', palette='Set1')
        plt.savefig(f"{figure_path}/{out}-BBKNN-re.png",dpi = 300)
    adata.write_h5ad(f"{h5ad}/{out}-bbknn.h5ad")
    return adata

##############################################
def combine(ad_list,number=0):
    #just align data;the number decide how to combine
    adata = sc.concat(ad_list,axis=number,join="outer",merge='same')
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata
def flow(adata,n_neighbors,n_pcs,resolution,out):
    adata = Run_Normalization(adata,n_neighbors,n_pcs,resolution,out)
    Run_batchRemove(adata,n_neighbors,n_pcs,resolution,out)


if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    samples = ["GBM27","GBM28","GBM29"]
    GBM_SC = []
    GBM_TE = []
    for i in samples:
        adata_SC = sc.read_h5ad(f"{h5ad}/{i}-SC-QC.h5ad")
        adata_TE = sc.read_h5ad(f"{h5ad}/{i}-TE-QC.h5ad")
        #TE may don't need QC,because it already filter after generating by scTE
        # adata_TE = sc.read_h5ad(f"{h5ad}/{i}-TE-raw.h5ad")
        GBM_SC.append(adata_SC)
        GBM_TE.append(adata_TE)
    SC = combine(GBM_SC,0)    
    TE = combine(GBM_TE,0)
    # from the same sample,how I align the TE and SC ?
    GBM = combine([SC,TE])
    del samples,GBM_SC,GBM_TE
    flow(SC,50,50,0.6,"SC-0.6")
    flow(TE,50,50,0.6,"TE-0.6")
    flow(GBM,50,50,0.6,"GBM-0.6")
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")