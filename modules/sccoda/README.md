# sccoda — cell-type differential abundance (scCODA)

Bayesian Dirichlet-Multinomial compositional analysis of cell-type
abundance across two conditions (Büttner et al. 2021, Nat Comm;
https://www.sc-best-practices.org/conditions/compositional/). One module
per atomic analysis step: this module wraps `bin/sccoda.py` (pertpy
scCODA) and is consumed by the scRNAseq subworkflow.

## Inputs

Per tissue x counter, one annotated h5ad with `obs['sample_id']` and
`obs['cell_type']` — the `scanpy_auto` output:

```
<outdir>/common/5_combine_h5ad/{tissue}/{tissue}_{counter}_auto.h5ad
```

When cluster-annotation overrides are configured for the tissue
(`Params.scanpy.<counter>.auto.cluster_annotations`), the rule also
depends on the `*.auto.h5ad.reannotated` marker so it never reads a
pre-override label set.

## Configuration (workflow config, `Params.sccoda.<counter>`)

| key | default | meaning |
|---|---|---|
| `enabled` | `false` | gate; node.py registers the report target only when enabled |
| `compare` | `[]` | `tissue:ref_sample:treat_sample[=label]` entries (ref FIRST). **Empty = shared from `Params.scanpy.<counter>.pseudobulk.compare`** (identical entry shape — one compare design serves both analyses); a hand-filled list wins |
| `ref_celltype` | `automatic` | scCODA reference cell type |
| `formula` | `condition` | design formula (must include `condition`) |
| `credible_threshold` | `0.9` | posterior cutoff for `is_credible` |
| `num_samples` | `2000` | NUTS posterior samples |
| `num_warmup` | `500` | NUTS warmup samples |

A report target is registered per tissue only when the counter block is
enabled AND at least one compare entry carries that tissue prefix.

## Outputs

```
<outdir>/common/5_combine_h5ad/{tissue}/sccoda/{counter}/
  sccoda_report.html                     summary report (declared output)
  da/{tissue}_{label}_da.csv             per-cell-type log2FC + final_prob
  da/{tissue}_{label}_effects_barplot.png
  da/{tissue}_{label}_stacked_barplot.png
```

`label` is the `=label` part of the compare entry (default
`{ref}_vs_{treat}`). Effects are log2(treat / ref); the CSV covariate
column tags the treatment with `[T.]`. The report is written on every
code path, including zero-credible-effect runs.

## CLI (standalone)

```bash
python modules/sccoda/bin/sccoda.py \
  --h5ad ovaries=/data/project/output/common/5_combine_h5ad/ovaries/ovaries_scTE_auto.h5ad \
  --compare ovaries:luanchao-21310:luanchao-11238=Aging_vs_Youth \
  --out-dir ./da_sccoda/ovaries/scTE \
  --ref-celltype automatic --formula condition \
  --credible-threshold 0.9 --num-samples 2000 --num-warmup 500
```

Run inside the module env: `apptainer exec --bind /data/16T:/data/16T
/data/16T/Database/env/sccoda/sccoda.sif python .../sccoda.py ...`.

## Environment

`sccoda.yaml` (python 3.11, scanpy, pertpy, jax, numpyro) — built once as
`/data/16T/Database/env/sccoda/sccoda.sif` and registered in the workflow
config `env` map as `"sccoda"`. The scanpy SIF does NOT carry pertpy, so
this module needs its own image.
