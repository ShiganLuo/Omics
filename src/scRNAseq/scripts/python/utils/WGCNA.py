# %% [markdown]
# ## 1.extract gene expression

import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("/home/lsg/Data/glioblastoma/output/new/h5ad/SC-bbknn.h5ad")
normalized_data = pd.DataFrame(adata.X.toarray(), 
                               index=adata.obs_names, 
                               columns=adata.var_names)

print(normalized_data.head())
normalized_data.to_csv("gene_expression.csv")


import PyWGCNA
geneExp = '/home/lsg/Data/glioblastoma/gene_expression.csv'
pyWGCNA_5xFAD = PyWGCNA.WGCNA(name='5xFAD', 
                              species='human', 
                              geneExpPath=geneExp, 
                              outputPath='',
                              save=True)
pyWGCNA_5xFAD.geneExpr.to_df().head(5)

pyWGCNA_5xFAD.preprocess()


