---
name: "scanpy-post-cluster-qc-annotate"
description: "Post-clustering QC and annotation workflow for scRNA-seq (ovaries/hystera, macaque). Iterative filtering, TE-dominated cluster handling, resolution tuning."
tags: ["scRNA-seq", "QC", "annotation", "scanpy", "macaque", "ovaries", "hystera"]
---

# Post-Clustering QC and Annotation Workflow

## Overview

After initial clustering (`scRNAseq.py --mode cluster`), clusters often contain low-quality cells or misidentified cell types. This skill documents the iterative QC and annotation workflow validated on macaque ovaries and hystera data.

**Core principle:** Cluster → Analyze → Flag → Filter cells (not clusters) → Re-cluster → Re-annotate.

## Workflow Summary

```
v1: Initial cluster (resolution=0.8)
  → Analyze cluster quality
  → Flag problematic clusters
  → Filter individual cells within flagged clusters
  → Save RAW counts (not processed)
v2: Re-cluster filtered data (resolution=0.8)
  → Check if problems resolved
  → If TE clusters still dominate (>20%), try higher resolution
v3: Higher resolution (resolution=1.5)
  → More granular cell type identification
  → Manual annotation if needed
```

## Step 1: Cluster Quality Analysis

For each cluster, check:

| Metric | Threshold | Flag |
|--------|-----------|------|
| Mean genes per cell | < 800 | low_genes |
| Mean UMI per cell | < 3000 | low_counts |
| Top markers are ENSMMUG | > 50% | high_unannotated |
| Top markers are MT genes | > 30% | high_mt_genes |
| Top markers are RPS/RPL | — | ribosomal_contamination |
| Top markers are TE elements | — | te_dominated (scTE only) |
| Cell type = Unknown | — | unknown |

**Code:**
```python
for cluster in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
    mask = adata.obs["leiden"] == cluster
    ct = adata.obs[mask]["cell_type"].iloc[0]
    mean_genes = adata.obs[mask]["n_genes_by_counts"].mean()
    mean_counts = adata.obs[mask]["total_counts"].mean()
    # Get top markers
    marker_df = sc.get.rank_genes_groups_df(adata, group=cluster)
    top20 = marker_df.head(20)["names"].tolist()
    pct_ensmmug = sum(1 for g in top20 if g.startswith("ENSMMUG")) / len(top20) * 100
```

## Step 2: Cell-Level Filtering (NOT Cluster Removal)

**Critical:** Never remove entire clusters. Filter individual cells within flagged clusters.

```python
CELL_THRESHOLDS = {"min_genes": 800, "min_counts": 3000, "max_pct_mt": 20.0}

keep_mask = pd.Series(True, index=adata.obs.index)
for cluster in flagged_clusters:
    cluster_mask = adata.obs["leiden"] == cluster
    cluster_cells = adata.obs[cluster_mask]
    cell_qc = (
        (cluster_cells["n_genes_by_counts"] >= CELL_THRESHOLDS["min_genes"])
        & (cluster_cells["total_counts"] >= CELL_THRESHOLDS["min_counts"])
        & (cluster_cells["pct_counts_mt"] <= CELL_THRESHOLDS["max_pct_mt"])
    )
    keep_mask[cluster_cells[~cell_qc].index] = False
```

## Step 3: Save RAW Counts

**Critical:** Save raw counts from the original merged file, NOT processed data.

```python
# CORRECT
adata_raw = ad.read_h5ad("merged.h5ad")
adata_filtered = adata_raw[filtered_cell_indices].copy()
adata_filtered.write_h5ad("filtered.h5ad")

# WRONG - processed data causes "Bin edges must be unique" error
adata_filtered.write_h5ad("filtered.h5ad")  # Don't do this
```

## Step 4: Re-cluster and Re-annotate

```bash
# Re-cluster
python scRNAseq.py --mode cluster \
    --input filtered.h5ad --output clustered.h5ad \
    --n-pcs 50 --n-neighbors 50 --resolution 0.8 \
    --batch-method harmony --batch-key sample_id --auto-n-pcs

# Re-annotate (filename must contain _clustered.h5ad)
python annotate_all.py --input clustered.h5ad --tissue ovaries --counter scTE
```

## Resolution Tuning

| Resolution | Clusters | Use Case |
|------------|----------|----------|
| 0.8 | ~17-19 | Default, good for initial clustering |
| 1.5 | ~26 | When TE clusters >20% and may contain multiple cell types |

**When to increase resolution:**
- ERVK_high or Alu_high clusters are >20% of total cells
- You suspect they contain multiple cell types masked by TE expression
- After filtering, some cell types are still not identified

**Ovaries example:**
- v1 (res=0.8): 19 clusters, ERVK_high=20.4%, Alu_high=11.7%
- v2 (res=0.8): 17 clusters after filtering, ERVK_high=27%
- v3 (res=1.5): 26 clusters, identified Cumulus, Pericyte, Mesothelial, T_cell

## Manual Annotation Override

When automatic annotation fails (marker overlap between cell types):

```bash
python annotate_all.py --input clustered.h5ad --tissue ovaries --counter scTE \
    --manual-annotation '{"1":"Myofibroblast"}'
```

**When to use:**
- Cell types share markers (Myofibroblast vs Smooth_muscle share ACTA2, MYH11, TAGLN)
- Automatic scoring can't distinguish cell types
- Strong biological evidence for specific annotation

**Code handling:**
```python
# Add new categories if needed
new_cats = set(manual_annotation.values()) - set(adata.obs["cell_type"].cat.categories)
if new_cats:
    adata.obs["cell_type"] = adata.obs["cell_type"].cat.add_categories(new_cats)
```

## TE-Dominated Clusters (scTE-specific)

### ERVK_high
- Top markers: MacERVK2_LTR1c, MacNERVK2-int, etc.
- Quality normal (genes 2500-2600)
- Does NOT match any known cell type
- Not Granulosa (FOXL2, CYP19A1, FSHR very low)
- Not Luteal (STAR, CYP11A1 very low)
- Keep as ERVK_high, do NOT assign cell type

### Alu_high
- Top markers: AluSz, AluSx, AluSx1, AluY, etc.
- May indicate low quality or TE dysregulation
- Keep as Alu_high

### Sample Distribution Check
```python
for ct in ['ERVK_high', 'Alu_high']:
    mask = adata.obs['cell_type'] == ct
    print(adata.obs[mask]['sample_id'].value_counts())
# If >70% from one sample → technical issue
```

## Tissue-Specific Considerations

### Ovaries (Macaca mulatta)
- Key cell types: Oocyte, Granulosa, Cumulus, Theca, Luteal, Stromal, Endothelial, Smooth_muscle, Pericyte, Macrophage, T_cell, NK_cell
- Myofibroblast: SM+Stromal hybrid, needs manual annotation
- Lymphatic_endo: MMRN1, CCL21, PROX1, LYVE1 — annotate separately from blood endothelial
- Cumulus: FST, NR5A2, PPARG, CRHBP — specialized granulosa
- Luteal: STAR, CYP11A1, HSD3B1 — steroidogenesis
- Very low MT% is normal for macaque data

### Hystera (Macaca mulatta)
- Key cell types: Epithelial_luminal, Epithelial_glandular, Stromal, Smooth_muscle, Pericyte, Endothelial, Macrophage, Uterine_NK, T_cell, Schwann_cell
- rRNA contamination: LSU_rRNA_eukarya, SSU_rRNA_eukarya in top markers → filter
- Epithelial splits into luminal and glandular
- Schwann_cell: PMP22, S100B, CDH19

## Annotation Algorithm (annotate_all.py)

**Specificity-weighted scoring:**
- Each marker's weight = 1 / (number of cell types it appears in)
- A gene unique to one cell type scores 1.0
- A gene shared by N cell types scores 1/N
- Rank-weighted: genes ranked higher get more weight (rank0=1.0, rank10=0.5, rank49=0.17)

**Formula:**
```
score = Σ (specificity_weight × rank_weight) for matched markers
rank_weight = 1.0 / (1.0 + rank × 0.1)
```

## Validated Cell Type Markers

### Ovary Markers (OVARY_MARKERS)
| Cell Type | Key Markers | Reference |
|-----------|-------------|-----------|
| Oocyte | ZP1, ZP2, ZP3, FIGLA, NOBOX, GDF9, BMP15 | Wang 2024 PMID:38191526 |
| Granulosa | FOXL2, CYP19A1, FSHR, AMH | Zhao 2026 PMID:42552301 |
| Cumulus | FST, NR5A2, PPARG, CRHBP | DEG analysis |
| Theca | CYP17A1, STAR, DLK1, INSL3 | Zhao 2026 PMID:42552301 |
| Luteal | STAR, CYP11A1, HSD3B1, PTCH2, GPC5 | Corpus luteum |
| Stromal | COL1A1, COL1A2, DCN, LUM, PDGFRA | Zhao 2026 PMID:42552301 |
| Endothelial | PECAM1, VWF, CDH5, ERG, KDR | Zhao 2026 PMID:42552301 |
| Lymphatic_endo | MMRN1, CCL21, PROX1, LYVE1 | Wigle 1999 PMID:10499794 |
| Smooth_muscle | ACTA2, MYH11, TAGLN, CNN1 | Zhao 2026 PMID:42552301 |
| Pericyte | RGS5, PDGFRB, NOTCH3, MCAM | Zhao 2026 PMID:42552301 |
| Myofibroblast | ACTA2, MYH11, TAGLN, POSTN, IGFBP5, SFRP1 | PMC 2024 |
| Mesothelial | MSLN, ITLN1, WT1, UPK3B | Nature 2022 |
| Macrophage | CD68, CD163, CSF1R, MRC1 | Zhao 2026 PMID:42552301 |
| T_cell | CD3E, CD3D, CD3G, CD4, CD8A | Zhao 2026 PMID:42552301 |

### Hystera Markers (HYSTERA_MARKERS)
| Cell Type | Key Markers | Reference |
|-----------|-------------|-----------|
| Epithelial_luminal | EPCAM, KRT18, CDH1, PAEP, SLC2A1 | Zhao 2026 PMID:42552301 |
| Epithelial_glandular | EPCAM, PAEP, SPP1, LIF, PAX8, ESR1 | Zhao 2026 PMID:42552301 |
| Stromal | COL1A1, DCN, LUM, WNT4, IGFBP1 | Zhao 2026 PMID:42552301 |
| Uterine_NK | NKG7, GNLY, KLRB1, CSF1, XCL1 | Zhao 2026 PMID:42552301 |
| Schwann_cell | PMP22, S100B, CDH19, PLP1 | Chen 2021 PMID:34237024 |

## Key Scripts

| Script | Location | Purpose |
|--------|----------|---------|
| scRNAseq.py | workflow/Omics/modules/scanpy/bin/ | Main pipeline (qc/merge/cluster/annotate) |
| annotate_all.py | output/.../5_combine_h5ad/ | Annotation with specificity-weighted scoring |
| test_scTE_qc.py | output/.../5_combine_h5ad/hystera/ | scTE post-clustering QC |
| test_cluster_qc_v2.py | output/.../5_combine_h5ad/hystera/ | Cell Ranger post-clustering QC |
| plot_cluster_qc.py | output/.../5_combine_h5ad/hystera/ | QC visualization |

## Common Pitfalls

1. **Removing entire clusters** → Filter individual cells instead
2. **Saving processed data** → Always save raw counts for re-clustering
3. **Not re-clustering after filtering** → Cluster IDs change, must re-cluster
4. **Not re-annotating after re-clustering** → Cluster IDs change, must re-annotate
5. **annotate_all.py filename** → Input must contain `_clustered.h5ad`
6. **Assigning cell types to TE clusters** → Keep as ERVK_high/Alu_high
7. **Trusting UMAP distances** → Check PCA distances for cell type relationships
8. **Not checking co-expression** → Mixed markers might be doublets or Myofibroblast
9. **Removing cells based on co-expression alone** → Check doublet score first
10. **Low MT% filtering** → Macaque data has very low MT%, don't use as primary filter

## Real-World Results

### Ovaries scTE (v3, resolution=1.5)
- 14943 cells, 26 clusters
- Removed: 4873 cells (24.6%) in v1 filtering
- Identified: ERVK_high (22.6%), Stromal (21.6%), Endothelial (19.3%), Smooth_muscle (12.9%), Pericyte (5.5%), Cumulus (4.5%), Alu_high (3.2%), Macrophage (2.4%), Proliferating (2.0%), Luteal (2.0%), Mesothelial (1.7%), Myofibroblast (1.3%), Epithelial (1.3%), T_cell (1.0%)

### Hystera scTE (v2)
- 41009 cells, 19 clusters
- Removed: 6455 cells (13.6%) — rRNA contamination and low quality
- Identified: Stromal (~27%), Smooth_muscle (~17%), Endothelial (~16%), Epithelial_luminal (~13%), Pericyte (~10%), Uterine_NK (~7%), Epithelial_glandular (~4%), Macrophage (~2%)

---

*Generated from macaque ovaries/hystera scRNA-seq analysis (2026-09-06)*
*Key references: Zhao et al. (2026) Cell Discov. PMID:42552301, Wang et al. (2024) Stem Cell Res Ther. PMID:38191526*
