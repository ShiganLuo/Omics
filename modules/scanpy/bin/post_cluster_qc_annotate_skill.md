---
name: "scanpy-post-cluster-qc-annotate"
description: "Post-clustering QC and annotation workflow for scRNA-seq (ovaries/hystera, macaque). Iterative filtering, TE-dominated cluster handling, resolution tuning."
tags: ["scRNA-seq", "QC", "annotation", "scanpy", "macaque", "ovaries", "hystera"]
---

# Post-Clustering QC and Annotation Workflow

## Overview

After initial clustering (`scRNAseq.py --mode cluster`), clusters often contain low-quality cells or misidentified cell types. This skill documents the iterative QC and annotation workflow validated on macaque ovaries and hystera data.

**Core principle:** Cluster → Analyze → Flag → Filter cells (not clusters) → Re-cluster → Re-annotate.

## Auto Mode (Integrated)

The full workflow is now integrated into `scRNAseq.py --mode auto`:

```bash
python scRNAseq.py --mode auto \
    --input merged.h5ad \
    --output annotated.h5ad \
    --tissue ovaries \
    --llm-method openai --llm-model gpt-4o \
    --resolution 0.8 --max-iterations 3 \
    --min-genes 800 --min-counts 3000 --max-pct-mt 20 \
    --batch-method harmony --batch-key sample_id --auto-n-pcs
```

**Per iteration:**
1. Cluster (same as mode_cluster)
2. AI-annotate each cluster (LLM + PubMed literature search for references)
3. Analyze cluster quality (programmatic + AI quality_flag)
4. Flag problematic clusters
5. Filter individual cells within flagged clusters
6. Repeat until clean or max_iterations reached

**Outputs:**
- `{output}.h5ad` — final annotated h5ad
- `{output}_reports/annotation_report.tsv` — per-cluster annotation
- `{output}_reports/references.tsv` — per-marker PubMed references

## Workflow Summary

```
v1: Initial cluster (resolution=0.8)
  → AI analyzes cluster quality (LLM + PubMed)
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

## Key Scripts

| Script | Location | Purpose |
|--------|----------|---------|
| scRNAseq.py | workflow/Omics/modules/scanpy/bin/ | Main pipeline (qc/merge/cluster/annotate/auto/advanced/de) |
| scRNAseq.py --mode auto | workflow/Omics/modules/scanpy/bin/ | Autonomous cluster→AI annotate→QC→filter→re-cluster |
| annotate_all.py | output/.../5_combine_h5ad/ | Annotation with specificity-weighted scoring |

## Common Pitfalls

1. **Removing entire clusters** → Filter individual cells instead
2. **Saving processed data** → Always save raw counts for re-clustering
3. **Not re-clustering after filtering** → Cluster IDs change, must re-cluster
4. **Not re-annotating after re-clustering** → Cluster IDs change, must re-annotate
5. **Assigning cell types to TE clusters** → Keep as ERVK_high/Alu_high
6. **Trusting UMAP distances** → Check PCA distances for cell type relationships
7. **Not checking co-expression** → Mixed markers might be doublets or Myofibroblast
8. **Removing cells based on co-expression alone** → Check doublet score first
9. **Low MT% filtering** → Macaque data has very low MT%, don't use as primary filter

---

*Generated from macaque ovaries/hystera scRNA-seq analysis (2026-09-09)*
*Key references: Zhao et al. (2026) Cell Discov. PMID:42552301, Wang et al. (2024) Stem Cell Res Ther. PMID:38191526*
