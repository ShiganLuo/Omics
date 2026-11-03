# function/sv_go — SV-affected genes GO/KEGG enrichment (macaque-aware)

GO/KEGG enrichment over the gene list carried by a VEP-annotated SV tab
(`<comparison>_annotated.tab`, produced by modules/sv `sv_exp_specific`).
Engine is `bin/sv_go-kegg.r` (clusterProfiler + `org.Hs/Mm/Mmu.eg.db`),
so **macaque is a first-class species** (`-s macaque`, KEGG organism `mcc`) —
unlike the gseapy/Enrichr route in modules/sv `sv_exp_enrichment`, which only
knows human/mouse.

Mirrors `modules/function/pseudobulk_go` (same module family, same wiring
shape): the parent subworkflow registers one `sv_go_kegg` job per SV
comparison and gates it behind `Params.function.sv_go`.

## Inputs

| key | meaning |
| --- | --- |
| `indir` | root containing `<comparison>/<comparison>_annotated.tab` (modules/sv outdir) |
| `comparisons` | comparison name map (from Params.sv.comparisons; drives target registration) |

## Parameters (`Params.function.sv_go`)

| key | default | meaning |
| --- | --- | --- |
| `enabled` | `false` | gate; parent wiring only registers outputs when true |
| `species` | — (required) | `human` / `mouse` / `macaque` |
| `impact` | `HIGH,MODERATE,LOW` | VEP IMPACT levels whose genes enter the list |
| `top` | `10` | top N terms plotted |
| `format` | `png` | plot format (`png`/`pdf`) |

## Outputs (per comparison)

`<outdir>/<comparison>/`

- `sv_genes.txt` — the extracted gene list (always written, even when empty)
- `go_enrichment.csv` / `kegg_enrichment.csv` — full clusterProfiler result tables
- `go_enrichment.png` / `kegg_enrichment.png` — top-N bar plots (empty list → "No enriched terms" panel)
- `sv_go_kegg.done` — sentinel
- `sv_go_kegg.sh` — replayable receipt

## Standalone example config

```json
{
  "ROOT_DIR": "/data/16T/luosg/workflow/Omics",
  "OUTDIR_ROOT": "/data/16T/luosg/output/luanchao_mulatta/Fiberseq",
  "env": {"env_dir": "/data/16T/Database/env", "function": "/data/16T/Database/env/function/function.sif"},
  "indir": "/data/16T/luosg/output/luanchao_mulatta/Fiberseq/variation/germline_sv_analysis",
  "outdir": "/data/16T/luosg/output/luanchao_mulatta/Fiberseq/variation/germline_sv_analysis/sv_go",
  "logdir": "/data/16T/luosg/output/luanchao_mulatta/Fiberseq/log/group",
  "comparisons": {"youth_vs_aged": {}},
  "Params": {
    "function": {
      "sv_go": {"species": "macaque"}
    }
  }
}
```

Run:

```bash
snakemake -s modules/function/sv_go/sv_go.smk --configfile sv_go.json \
  --cores 4 --sdm apptainer \
  --singularity-args "--bind /data/16T/Database,/data/16T/luosg/output"
```
