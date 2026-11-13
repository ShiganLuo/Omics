# sigtest

Per-feature two-group significance testing with star labels, multiple-testing
correction, and a bracket-annotated grouped bar plot. Unifies the legacy
`src/common/statistics_plot/Python/statistical_test.py` (assumption checks +
automatic sample-based test selection) and `bar_significance.py` (count-based
contingency tests, p-to-star labels, display-point brackets) into one module.

## Input

One long-format TSV/CSV via config key `sigtest_table`:

- sample mode (`method=auto|t-test|welch|mannwhitney|wilcoxon`): one row per
  replicate `feature, group, value`. `auto` checks Shapiro-Wilk normality and
  Levene equal-variance per feature, then picks paired t / Wilcoxon
  (`paired=true`), Student t, Welch t, or Mann-Whitney U. Reports Cohen's d
  or rank-biserial r.
- count mode (`method=chi2|fisher`): exactly one summarized count per
  feature x group; builds 2x2 tables against group totals. chi2 falls back
  to Fisher when expected frequencies hit zero.

## Outputs

Under `outdir/{contrast}/`:

- `significance_test.tsv` — per feature: n/mean/std per group, test used,
  statistic, p_value, p_adjusted (bh/bonferroni/none), stars, significant,
  effect size, note (skipped features keep their reason here)
- `significance_bars.png` — grouped bars (SEM error bars) with significance
  brackets drawn in display points; y headroom auto-fitted

Stars: p<1e-4 `****`, <1e-3 `***`, <1e-2 `**`, <0.05 `*`, else `ns`.

## Config keys

`sigtest_table` (required file path), `contrasts` (list; the panel rule
`sigtest_result` expands over it), and `Params.sigtest.*`:
`feature_col` (default `feature`), `group_col` (`group`), `value_col`
(`value`), `groups` ([g1, g2]; default: first two distinct labels), `method`
(`auto`), `paired` (false), `alpha` (0.05), `correction` (`none`), `formats`
(repeatable image formats, default `png`; an explicit empty value = no
figure), `ylabel` (`value`), `legend` ([l1, l2]), `annotate_ns` (false),
`no_errorbars` (false).

## Example

Standalone module run (config needs `ROOT_DIR`, `OUTDIR_ROOT`, `env`
mapping, as any `*_config` dict would provide):

```bash
snakemake -s modules/sigtest/sigtest.smk --configfile=cfg.json \
    --target-jobs sigtest_result:contrast=P6_vs_P20 --cores 1
```

with `cfg.json` (the exact config used in verification):

```json
{
    "ROOT_DIR": "/data/16T/luosg/workflow/Omics",
    "OUTDIR_ROOT": "<workdir>/module_run",
    "indir": "<workdir>",
    "outdir": "<workdir>/module_run/sigtest",
    "logdir": "<workdir>/module_run/log",
    "sigtest_table": "<workdir>/expr_long.tsv",
    "contrasts": ["P6_vs_P20"],
    "env": {"sigtest": "/data/16T/Database/env/sv/sv.sif"},
    "resources": {"sigtest_run": {"threads": 1, "mem_mb": 4000}},
    "Params": {"sigtest": {"feature_col": "gene", "group_col": "group",
                           "value_col": "expr", "method": "auto",
                           "correction": "bh", "ylabel": "Expression (TPM)"}}
}
```

Direct bin-script use:

```bash
python modules/sigtest/bin/sigtest.py -i expr_long.tsv -o sig_out/ \
    -f gene -g Control Mutant -m auto -c bh -y "Expression (TPM)"
python modules/sigtest/bin/sigtest.py -i sv_counts.tsv -o sig_out/ \
    -f svtype -m chi2 --no-errorbars -y "SV count"
```

The module is group-level (no `{sample_id}` wildcard): wire it from a
subworkflow with `logdir` pointing at `{logdir}/group`, and derive
`contrasts` from the workflow's group pairs.
