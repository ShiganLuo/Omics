# burden module

Mutation-burden calculations for somatic variant calls: tumor mutational
burden (TMB) and indel burden from filtered somatic SNV/InDel VCFs,
structural-variant (SV) burden from somatic SV VCFs, and copy-number (CNV)
burden from CNVkit segment tables, plus a cross-sample summary merger.

Logic ported and hardened from
`/data/16T/luosg/workflow/genome-stability/index/burden` (hard-coded paths
removed, CLI parameterized, mixed-language comments translated, atomic
output writes, provenance columns added).

## Rules

| Rule | Scope | Input | Output |
|---|---|---|---|
| `burden_snv_indel` | per sample | filtered somatic VCF (`sample_snv_vcf_dict`) | `{outdir}/{sample}/{sample}.snv_indel_burden.tsv` |
| `burden_sv` | per sample | somatic SV VCF (`sample_sv_vcf_dict`) | `{outdir}/{sample}/{sample}.sv_burden.tsv` |
| `burden_cnv` | per sample | CNVkit `.cns` (`sample_cnv_dict`) | `{outdir}/{sample}/{sample}.cnv_burden.tsv` |
| `burden_aggregate` | group | all per-sample TSVs | `{outdir}/burden_summary.tsv` |
| `burden_result` | panel | all of the above | — |

Samples are **tumor** sample ids; for paired calling the config maps each
tumor to its pair VCF (e.g. `KSR-1 -> .../E14-1_vs_KSR-1/E14-1_vs_KSR-1.filtered.vcf.gz`).

## Config keys

| Key | Meaning |
|---|---|
| `sample_snv_vcf_dict` | `{sample: filtered somatic VCF}` |
| `sample_sv_vcf_dict` | `{sample: somatic SV VCF}` |
| `sample_cnv_dict` | `{sample: CNVkit .cns}` |
| `Params.burden.min_vaf` | minimum tumor allele fraction (default 0.05) |
| `Params.burden.size_mb` | fixed coding territory in Mb (default 34.0) |
| `Params.burden.gtf_feature` | `exon` or `CDS` for the GTF denominator (default `exon`) |
| `Params.burden.include_sex_chroms` | count X/Y too (default false: autosomes only) |
| `Params.burden.min_abs_log2` | CNV altered threshold when no `call` column (default 0.3) |
| `Params.burden.annotate` | VEP annotation + meaningful-variant filtering for TMB (default true) |
| `Params.burden.vep_cache` | VEP cache root directory (default `~/.vep`) |
| `Params.burden.vep_species` | VEP species name (auto-filled from the genome by node.py) |
| `Params.burden.vep_assembly` | VEP assembly name (auto-filled from the genome by node.py) |
| `genome.targets` | capture targets BED; becomes the per-Mb denominator |
| `genome.gtf` | gene annotation GTF; merged feature length becomes the denominator |

Denominator precedence for the per-megabase rates: `genome.targets` >
`genome.gtf` > `Params.burden.size_mb`. The chosen source and territory are
recorded per row (`territory_mb`, `denominator`). Variant counting and
territory computation share one chromosome universe (canonical human/mouse
autosomes; X/Y only with `include_sex_chroms`).

SNV/InDel counting keeps PASS (or `.`) records only (override with
`--include-filtered`) and applies a VAF floor on the tumor sample column
(FORMAT `AF`, falling back to `AD`). With annotation enabled (default) the
VCF is first annotated with Ensembl VEP (`bin/vep_ann.py`, cache-based
offline; the annotated copy is kept as `<sample>.vep.vcf`) and only
meaningful variants are counted — VEP `CSQ` consequences in the
protein-impacting / core-splice whitelist (missense, frameshift, stop
gained/lost, splice acceptor/donor, inframe indels, ...). The annotated
VCF must carry `CSQ`; otherwise the run raises. With `--no-annotate` VEP
is skipped and every passing allele counts without consequence filtering.
CNV segments count as altered via the `call` column when present, else
`|log2| >= min_abs_log2`.

## Outputs

- `{sample}.snv_indel_burden.tsv`: `sample snv_count indel_count tmb_count
  frameshift_count tmb_per_mb snv_per_mb indel_per_mb territory_mb denominator`
- `{sample}.sv_burden.tsv`: `sample sv_total sv_DEL sv_DUP sv_INV sv_INS
  sv_BND sv_OTHER` (INFO/SVTYPE buckets)
- `{sample}.cnv_burden.tsv`: `sample cnv_total_segments cnv_gain_segments
  cnv_loss_segments cnv_altered_length_mb cnv_altered_gene_count`
- `burden_summary.tsv`: outer join of all per-sample tables on `sample`

Every table is written even when all counts are zero (header + zero row).

## Standalone CLI

```bash
# TMB with VEP annotation + meaningful-variant filtering (default)
python modules/burden/bin/burden.py snv-indel \
    -i somatic_snv_indel/E14-1_vs_KSR-1/E14-1_vs_KSR-1.filtered.vcf.gz \
    -s KSR-1 -o KSR-1.snv_indel_burden.tsv -b targets.bed -v 0.05 \
    -c ~/.vep -p mus_musculus -y GRCm39

# Same, without annotation (count every passing allele)
python modules/burden/bin/burden.py snv-indel \
    -i somatic_snv_indel/E14-1_vs_KSR-1/E14-1_vs_KSR-1.filtered.vcf.gz \
    -s KSR-1 -o KSR-1.raw_burden.tsv --no-annotate

python modules/burden/bin/burden.py sv \
    -i somatic_sv/E14-1_vs_KSR-1/results/variants/somaticSV.vcf.gz \
    -s KSR-1 -o KSR-1.sv_burden.tsv

python modules/burden/bin/burden.py cnv \
    -i cnv/KSR-1/KSR-1.cns -s KSR-1 -o KSR-1.cnv_burden.tsv -m 0.3

python modules/burden/bin/burden.py aggregate \
    -n *.snv_indel_burden.tsv -v *.sv_burden.tsv -c *.cnv_burden.tsv \
    -o burden_summary.tsv
```

All subcommands accept `-l/--log` (log file) and `-w/--overwrite` (rewrite
existing output; without it an existing output is kept).

## Environment

Conda env `burden.yaml` (python + pandas). Build the SIF once:

```bash
apptainer build /data/16T/Database/env/burden/burden.sif \
    modules/burden/burden.def
```

and register it in the workflow config `env` map as
`"burden": "/data/16T/Database/env/burden/burden.sif"` (already present in
`config/Mutation.json`).
