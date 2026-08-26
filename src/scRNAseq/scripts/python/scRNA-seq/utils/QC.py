import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
import anndata as ad
import matplotlib.pyplot as plt
import scrublet as scr
from math import sqrt
import pandas as pd
import datetime
figure_path = "/home/lsg/Data/glioblastoma/output/new/figure/QC"
h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"
tableout = "/home/lsg/Data/glioblastoma/output/new/table"
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
def lowquality(adata,sample,fig=0):
    print(f"------------------------{sample}---------------------")
    ### metrix caculate
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
    # print(f"统计信息：\n{adata}")
    ####plot above
    
    if fig:
        sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
        plt.savefig(f"{figure_path}/{sample}-total_counts.png")
        sc.pl.violin(adata, "pct_counts_mt",show=False)
        plt.savefig(f"{figure_path}/{sample}-pct_counts_mt.png")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",show=False)
        plt.savefig(f"{figure_path}/{sample}-scatter-before.png")
        plt.close()

    ### qc function
    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    # print(f"outlier:\n{adata.obs.outlier.value_counts()}")
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
    )
    # print(f"outlier-mt:\n{adata.obs.mt_outlier.value_counts()}")
    # print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    # print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    ### plot after filter
    if fig:
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
        plt.savefig(f"{figure_path}/{sample}-scatter-after.png")
        plt.close()
    return adata
# def ambientRNA():
sim_doublet_ratio = 2
def Compute_Doublet(adata,rate,out,fig=0,tab=0):
    counts_matrix = adata.to_df()
    n_cells = adata.shape[0]
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=rate,
                        n_neighbors = round(0.5 * sqrt(n_cells)),
                        sim_doublet_ratio = sim_doublet_ratio)
    ### annoy-1.15.1
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
    if fig:
        scrub.plot_histogram()
        plt.savefig(f"{figure_path}/{out}-histogram.png",dpi=300, bbox_inches='tight')
        plt.close()
    scrub.set_embedding('TSNE', scr.get_tsne(scrub.manifold_obs_, 0.5,10))
    if fig:
        scrub.plot_embedding('TSNE', order_points=True)
        plt.savefig(f"{figure_path}/{out}-predicted_doublets.png")
        plt.close()
    
    out_df = pd.DataFrame()
    out_df['barcodes'] = counts_matrix.index
    out_df['doublet_scores'] = doublet_scores
    out_df['predicted_doublets'] = predicted_doublets
    if tab:
        out_df.to_csv(f"{tableout}/{out}-doublet.txt", index=False,header=True)
        out_df.head()
    return out_df,doublet_scores,predicted_doublets
def Filter_cells(ad,doublet_scores,predicted_doublets):
    ad.obs["doublet_scores"] = doublet_scores
    ad.obs["predicted_doublets"] = predicted_doublets
    #~  可以作为取反的功能
    ad = ad[~ad.obs.predicted_doublets, :]
    return ad


if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    samples = ["GBM27","GBM28","GBM29"]
    for i in samples:
        adata_SC = sc.read_h5ad(f"{h5ad}/{i}-SC-raw.h5ad")
        adata_TE = sc.read_h5ad(f"{h5ad}/{i}-TE-raw.h5ad")
        adata_SC = lowquality(adata_SC,f"{i}-SC")
        adata_TE = lowquality(adata_TE,f"{i}-TE")
        SC_df,SCdoublet_scores,SCpredicted_doublets = Compute_Doublet(adata_SC,0.06,f"{i}-SC")
        TE_df,TEdoublet_scores,TEpredicted_doublets = Compute_Doublet(adata_TE,0.06,f"{i}-TE")
        adata_SC = Filter_cells(adata_SC,SCdoublet_scores,SCpredicted_doublets)
        adata_TE = Filter_cells(adata_TE,TEdoublet_scores,TEpredicted_doublets)
        adata_SC.write(f'{h5ad}/{f"{i}-SC"}-QC.h5ad')
        adata_TE.write(f'{h5ad}/{f"{i}-TE"}-QC.h5ad')
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")