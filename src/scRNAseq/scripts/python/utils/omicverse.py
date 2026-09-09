# %%
import omicverse as ov
import scanpy as sc
import scvelo as scv

ov.utils.ov_plot_set()

# %%
adata = ov.read('/home/lsg/Data/glioblastoma/kang.h5ad')
adata

# %%
#quantity control
adata=ov.pp.qc(adata,
              tresh={'mito_perc': 0.2, 'nUMIs': 500, 'detected_genes': 250},doublets_method='sccomposite',)
#normalize and high variable genes (HVGs) calculated
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=2000,)

#save the whole genes and filter the non-HVGs
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]

#scale the adata.X
ov.pp.scale(adata)

#Dimensionality Reduction
ov.pp.pca(adata,layer='scaled',n_pcs=50)

adata.obsm['X_mde']=ov.utils.mde(adata.obsm['scaled|original|X_pca'])

# %%
ov.utils.embedding(adata,
                   basis='X_mde',
                    frameon='small',
                   color=['replicate','label','cell_type'])

# %%
adata_harmony=ov.single.batch_correction(adata,batch_key='replicate',methods='harmony',n_pcs=50)
adata


# %%
adata.obsm["X_mde_harmony"] = ov.utils.mde(adata.obsm["X_harmony"])

# %%
ov.utils.embedding(adata,
                basis='X_mde_harmony',frameon='small',
                color=['replicate','label','cell_type'],show=False)

# %%
test_adata=adata[adata.obs['cell_type'].isin(['CD4 T cells'])]
dds=ov.bulk.pyDEG(test_adata.to_df(layer='lognorm').T)
dds.drop_duplicates_index()
print('... drop_duplicates_index success')

# %%
treatment_groups=test_adata.obs[test_adata.obs['label']=='stim'].index.tolist()
control_groups=test_adata.obs[test_adata.obs['label']=='ctrl'].index.tolist()
dds.deg_analysis(treatment_groups,control_groups,method='ttest')

# %%
# -1 means automatically calculates
dds.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=10)

# %%
dds.plot_volcano(title='DEG Analysis',figsize=(4,4),
                 plot_genes_num=8,plot_genes_fontsize=12,)

# %%
dds.plot_boxplot(genes=['NELL2','TNFSF13B'],treatment_groups=treatment_groups,
                control_groups=control_groups,figsize=(2,3),fontsize=12,
                 legend_bbox=(2,0.55))

# %%
ov.utils.embedding(test_adata,
                   basis='X_mde_harmony',
                    frameon='small',
                   color=['label','NELL2','TNFSF13B'])

# %%
test_adata=adata[adata.obs['cell_type'].isin(['CD4 T cells'])]
dds=ov.bulk.pyDEG(test_adata.to_df(layer='counts').T)
dds.drop_duplicates_index()
print('... drop_duplicates_index success')

treatment_groups=test_adata.obs[test_adata.obs['label']=='stim'].index.tolist()
control_groups=test_adata.obs[test_adata.obs['label']=='ctrl'].index.tolist()
dds.deg_analysis(treatment_groups,control_groups,method='DEseq2')

# %%
# -1 means automatically calculates
dds.foldchange_set(fc_threshold=1.2,
                   pval_threshold=0.05,
                   logp_max=10)

# %%
dds.plot_volcano(title='DEG Analysis',figsize=(4,4),
                 plot_genes_num=8,plot_genes_fontsize=12,)

# %%
dds.plot_boxplot(genes=['IFI6','RGCC'],treatment_groups=treatment_groups,
                control_groups=control_groups,figsize=(2,3),fontsize=12,
                 legend_bbox=(2,0.55))

# %%
ov.utils.embedding(test_adata,
                   basis='X_mde_harmony',
                    frameon='small',
                   color=['label','IFI6','RGCC'])

# %%
meta_obj=ov.single.MetaCell(adata,use_rep='X_harmony',n_metacells=250,
                           use_gpu=True)

# %%
meta_obj.initialize_archetypes()


