from anndata import AnnData
import scanpy as sc
TE1 = sc.read_10x_mtx(
    "/home/lsg/Data/glioblastoma/output/star/GBM27Solo.out/Gene/raw",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)
print(TE1)
print(TE1.obs)
print(TE1.var)
# TE2 = sc.read_10x_mtx(
#     "/home/lsg/Data/glioblastoma/output/star/GBM28Solo.out/Gene/raw",  # the directory with the `.mtx` file
#     var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
#     cache=True,  # write a cache file for faster subsequent reading
# )
# TE3 = sc.read_10x_mtx(
#     "/home/lsg/Data/glioblastoma/output/star/GBM29Solo.out/Gene/raw",  # the directory with the `.mtx` file
#     var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
#     cache=True,  # write a cache file for faster subsequent reading
# )
