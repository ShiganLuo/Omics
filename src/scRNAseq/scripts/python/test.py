from anndata import AnnData
import scanpy as sc
TE1 = sc.read_10x_mtx(
    "/home/lsg/Data/glioblastoma/output/star/GBM27Solo.out/Gene/raw",
    var_names="gene_symbols",
    cache=True,
)
print(TE1)
print(TE1.obs)
print(TE1.var)
