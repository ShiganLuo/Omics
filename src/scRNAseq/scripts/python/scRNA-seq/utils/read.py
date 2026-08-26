import numpy as np
import scanpy as sc
import anndata
import pandas as pd
import scipy as sp
import argparse
import datetime
h5ad = "/home/lsg/Data/glioblastoma/output/new/h5ad"
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

def sparsify(filename,sample):
    data = pd.read_csv(filename, index_col=0, header=0)
    data.index = data.index.astype(str)
    genes = data.columns
    cells = data.index
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')

    '''
    oh = open('gene_names.{0}.tsv'.format(os.path.split(filename)[1]), 'w')
    for g in genes:
        oh.write('%s\n' % g)
    oh.close()
    '''

    print('Loaded {0}'.format(filename))
    ad = anndata.AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data
    ad.obs['n_genes'] = (ad.X > 0).sum(axis=1).A1  # 计算每个细胞的基因数量
    ad.obs['n_counts'] = ad.X.sum(axis=1).A1 #axis=1沿行方向
    ad.obs['sample'] = np.full(ad.n_obs, sample)
    # ad.var_names_make_unique()
    # ad.obs_names_make_unique()
    return ad
# the output of cellranger
def cellranger(filename,sample):
    ad = sc.read_10x_mtx(
    filename,  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
    ad.obs['sample'] = np.full(ad.n_obs, sample)
    # ad.var_names_make_unique()
    # ad.obs_names_make_unique()
    return ad
# the output of scTE
def scTE(filename,sample):
    data = pd.read_csv(filename, index_col=0, header=0)
    data.index = data.index.astype(str)
    genes = data.columns
    cells = data.index
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')
    ad = anndata.AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data
    ad.obs['sample'] = np.full(ad.n_obs, sample)
    # ad.obs['n_TE'] = (ad.X > 0).sum(axis=1).A1  # 计算每个细胞的转座子数量
    # ad.obs['n_counts'] = ad.X.sum(axis=1).A1 #axis=1沿行方向
    # ad.var_names_make_unique()
    # ad.obs_names_make_unique()
    return ad

if __name__ == '__main__':
    start = datetime.datetime.now()
    # parser = argparse.ArgumentParser(description="A script to process cell cycle data.")
    # parser.add_argument('--options', type=str, required=True, help='options to execute procedure')
    # parser.add_argument('--input', type=str, required=True, help='Path to input file')
    # parser.add_argument('--out', type=str, required=True, help='Path to out file')
    samples = ["GBM27","GBM28","GBM29"]
    for i in samples:
        SC_file = f"/home/lsg/Data/glioblastoma/output/cellranger/{i}/outs/filtered_feature_bc_matrix"
        TE_file = f"/home/lsg/Data/glioblastoma/output/scTE/cellranger_yes/{i}.csv"
        print(f"read: {SC_file}")
        print(f"read: {TE_file}")
        ad_SC = cellranger(SC_file,i)
        ad_TE = scTE(TE_file,i)
        ad_SC.write_h5ad(f"{h5ad}/{i}-SC-raw.h5ad")
        ad_TE.write_h5ad(f"{h5ad}/{i}-TE-raw.h5ad")
    end = datetime.datetime.now()
    print("程序运行时间："+str((end-start).seconds/3600)+"h")