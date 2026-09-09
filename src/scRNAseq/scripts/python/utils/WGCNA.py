# %% [markdown]
# ## 1.extract gene expression

# %%
import scanpy as sc
import pandas as pd

# Assuming 'adata' is your AnnData object
# Extract the normalized expression matrix
adata = sc.read_h5ad("/home/lsg/Data/glioblastoma/output/new/h5ad/SC-bbknn.h5ad")
normalized_data = pd.DataFrame(adata.X.toarray(), 
                               index=adata.obs_names, 
                               columns=adata.var_names)

# View or save the data
print(normalized_data.head())  # Display the first few rows
normalized_data.to_csv("gene_expression.csv")  # Save as CSV

# %% [markdown]
# ## 2.pyWGCNA

# %%
import PyWGCNA
geneExp = '/home/lsg/Data/glioblastoma/gene_expression.csv'
pyWGCNA_5xFAD = PyWGCNA.WGCNA(name='5xFAD', 
                              species='human', 
                              geneExpPath=geneExp, 
                              outputPath='',
                              save=True)
pyWGCNA_5xFAD.geneExpr.to_df().head(5)

# %%
pyWGCNA_5xFAD.preprocess()


