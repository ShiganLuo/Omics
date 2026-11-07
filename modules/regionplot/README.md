# regionplot

Region-focused multi-sample track stacking plots. Given a regions BED and a
set of per-sample bigWigs already produced by the workflow, renders one
pyGenomeTracks figure per region: one wiggle line per sample, optional group
coloring + legend (from the meta `group` column or hand-filled
`Params.regionplot.sample_group_dict`), and an optional gene-model track from
the genome GTF. One image per region panel (PNG + SVG by default).

## Inputs

- Per-sample bigWigs via the `tracks` contract (wired by the subworkflow,
  not user-facing): `[{label, dir, suffix}]` expands to
  `{dir}/{sample_id}/{sample_id}{suffix}`.
- `Params.regionplot.regions_bed`: BED of `chrom/start/end[/name]` rows;
  col4 (sanitized, deduplicated) names the panel, else `chrom_start_end`.

## Outputs (`{outdir}/regionplot/{panel}/`)

- `<panel>.png` / `<panel>.svg` — the stacked-track figure
- `tracks.ini` — pyGenomeTracks config (rerun to restyle)
- `track_sources.tsv` — audit table (sample / track group / bw path)
- `regions.bed` — copy of the panel's region row
- `group_map.tsv` + `group_legend.{png,svg}` — only when groups exist

## Config keys (module `Params.regionplot`)

| key | default | meaning |
|---|---|---|
| `regions_bed` | null | regions BED; required when enabled |
| `sample_group_dict` | {} | sample -> group; auto-filled from meta `group` by node.py (`_register_regionplot`), hand-filled config wins |
| `sample_order` | [] | plotting order; empty = meta order |
| `formats` | ["png","svg"] | output image formats |
| `gtf_track` | true | add gene-model track when genome GTF exists |
| `plot_height_per_track` | 1.5 | inches per sample track |
| `plot_width` | 12 | figure width in inches |

Group colors are assigned per group (tab10 family, lightness-ramped shades
within a group); track labels read `{sample_id} ({group})` when grouped.

## Wiring (which workflows consume it)

Wired by every workflow that already produces bigWigs — Fiberseq
(accessibility / m6a_density / coverage), PacVar (coverage), CLIP (plus /
minus), MERIP / ncRNAseq / PeakCalling (igv bigwig). The subworkflow block is
unconditional; run/skip is gated by `Params.regionplot.enabled` (default
false) in node.py's `_register_regionplot`, which also auto-fills the group
map and panel outfiles.

## Standalone run

```
python modules/regionplot/bin/regionplot.py \
    -r regions.bed --panel MYC_locus \
    -t track_table.tsv -g group_map.tsv \
    --gtf genes.gtf --gtf-track \
    -o out/MYC_locus --formats png svg
```

`track_table.tsv`: `label<TAB>sample_id<TAB>bw_path` (group-major order =
plot order). Engine runs in the `regionplot.sif` container
(`pygenometracks=3.9`, built from `regionplot.def`; build dir
`/data/16T/Database/env/regionplot/`).
