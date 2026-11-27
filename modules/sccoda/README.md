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
| `compare` | `{}` | JSON design `{tissue: {contrast: entry}}`. `entry` = DESeq2 `group_pairs` keys (`control_group_name`, `experimental_group_name`, `control_samples`, `experimental_samples`) or an ordered `groups: [{name, samples}]` list for a multi-group contrast (first group = scCODA base level; each later group vs base; stacked barplot draws one consecutive bar per group and sums a group's samples into its bar). The contrast key names the outputs. **Empty = shared from `Params.scanpy.<counter>.scanpy_pseudobulk.compare`** (bare-sample string entries, converted to single-sample group_pairs); a hand-filled design wins |
| `ref_celltype` | `automatic` | scCODA reference cell type |
| `formula` | `condition` | design formula (must include `condition`) |
| `credible_threshold` | `0.9` | posterior cutoff for `is_credible` |
| `num_samples` | `2000` | NUTS posterior samples |
| `num_warmup` | `500` | NUTS warmup samples |
| `formats` | `["png"]` | plot output formats written simultaneously (matplotlib savefig set: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba); each plot is written per format as `da/{prefix}_<name>.<fmt>`. CLI `--format` (repeatable; `-f` is taken by `--formula`) |

A report target is registered per tissue only when the counter block is
enabled AND at least one compare entry carries that tissue prefix.

## Outputs

```
<outdir>/diff/abundance/{tissue}/{counter}/
  compare_design.json                        design slice this run consumed
  sccoda_report.html                          summary report (declared output)
  da/da_sccoda_{contrast}_da.csv              per-cell-type log2FC + final_prob
  da/da_sccoda_{contrast}_effects_barplot.png
  da/da_sccoda_{contrast}_stacked_barplot.png
```

`{contrast}` is the design's contrast key. Effects are
log2FC = log2(group / control); the CSV covariate column tags non-base
groups with `[T.]`. Figure conventions:

- **stacked barplot** — no title; one bar per group (samples of a group
  summed), x ticks labeled with the group names. DA status: the `↑`/`↓`
  arrows sit DIRECTLY AFTER each cell-type legend name (centered on the
  name's row, so every arrow unambiguously belongs to its cell type),
  one glyph per comparison, in that comparison's color (e.g. a blue `↓`
  = youth vs aged — color, not position, identifies the comparison, so
  multi-group contrasts stay unambiguous). Below, a DA key legend maps
  each color to its meaning: one row per comparison, the `↑↓` glyph pair
  in that color as the handle (leftmost, left-aligned with the cell-type
  swatches) and the meaning as the label (`ctrl vs exp`; ↑ = higher in
  the experimental group, ↓ = lower in it).
- **effects barplot** — no title; one panel per compared group in group
  order, sorted by log2FC per panel. The x label spells out the
  direction (`log2-fold change (youth vs aged)`, control FIRST); both
  axis ends are annotated relative to the EXPERIMENTAL group
  (`← lower in aged` / `higher in aged →`). No per-cell-type arrow
  marks — the bars' side already shows the direction. A panel whose
  effects are all zero carries an explicit grey note ("all effects are
  zero (nothing selected)") instead of rendering blank.

The report is written on every code path, including zero-credible-effect
runs.

## CLI (standalone)

```bash
python modules/sccoda/bin/sccoda.py \
  --h5ad ovaries=/data/16T/luosg/output/luanchao_mulatta/scRNAseq/common/5_combine_h5ad/ovaries/ovaries_scTE_auto.h5ad \
  --compare-json /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/abundance/ovaries/scTE/compare_design.json \
  --out-dir /data/16T/luosg/output/luanchao_mulatta/scRNAseq/diff/abundance/ovaries/scTE \
  --ref-celltype automatic --formula condition \
  --credible-threshold 0.9 --num-samples 2000 --num-warmup 500
```

`compare_design.json` (pairwise + multi-group contrast in one design):

```json
{
  "ovaries": {
    "youth_vs_aged": {
      "control_group_name": "youth",
      "experimental_group_name": "aged",
      "control_samples": ["luanchao-21310-10XSC3"],
      "experimental_samples": ["luanchao-11238-10XSC3"]
    },
    "aging_series": {
      "groups": [
        {"name": "youth", "samples": ["s1", "s2"]},
        {"name": "aged", "samples": ["s3", "s4"]}
      ]
    }
  }
}
```

Run inside the module env: `apptainer exec --bind /data/16T:/data/16T
/data/16T/Database/env/sccoda/sccoda.sif python .../sccoda.py ...`.

## Environment

`sccoda.yaml` (python 3.11, scanpy, pertpy, jax, numpyro) — built once as
`/data/16T/Database/env/sccoda/sccoda.sif` and registered in the workflow
config `env` map as `"sccoda"`. The scanpy SIF does NOT carry pertpy, so
this module needs its own image.
