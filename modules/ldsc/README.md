# ldsc

Partitioned heritability (LDSC) interface for FIRE peak sets [1]. Implements
the standard partitioned-LDSC command sequence (munge_sumstats -> ldsc.py
--h2-cts) with all external resources supplied via the config JSON.

## Rules

- `ldsc_annot`   : split a peak BED into per-cluster annotation BEDs (bin/ldsc_annot.py)
- `ldsc_munge`   : `munge_sumstats.py` per GWAS trait
- `ldsc_h2_cts`  : `ldsc.py --h2-cts` per trait (cluster-specific heritability)
- `ldsc_result`  : panel of `*.cts.results`

## Files to fill in config JSON (Params.ldsc)

| Key | Content |
|---|---|
| `annot_bed` | peak set to test (e.g. merged FIRE peaks; name col can encode clusters) |
| `annot_prefix` | output prefix (default `peaks`) |
| `annot_cluster_from_name` | true to split by BED name column |
| `gwas` | `{trait_name: /path/to/sumstats.tsv}` (add one entry per trait) |
| `merge_alleles` | HapMap3 SNP list for `munge_sumstats.py --merge-alleles` |
| `ref_ld_chr` | reference LD score prefix dir (e.g. eur_w_ld_chr/) |
| `w_ld_chr` | regression weights dir (e.g. eur_w_ld_chr/) |

`Procedure.ldsc` / `Procedure.munge_sumstats` override executable paths.

## Sources

- [1] https://doi.org/10.1016/j.crmeth.2024.100911 (LD score enrichment analysis)
