# Omics Workflow

English | [中文](README.md)

![Snakemake](https://img.shields.io/badge/Snakemake-9.16.3%2B-1f7a4d)
![Python](https://img.shields.io/badge/Python-3.12-3776ab)
![Apptainer](https://img.shields.io/badge/Apptainer-1.5.4-e2761b)
![Workflows](https://img.shields.io/badge/workflows-16-6a5acd)
![Modules](https://img.shields.io/badge/modules-86-2f7d32)

16 omics analysis workflows behind one entry point: `python run.py -w <workflow>`. Transcriptomics, single-cell, ChIP-seq, mutation, long-read, proteomics, and population genomics all share the same reusable modules — run whole pipelines in one shot, or any single module via `run_module.py`. Every run's parameters are persisted in `raw.json` for full provenance, conda / Apptainer environments assemble themselves, and AI can pick up troubleshooting and extension directly.

> Documentation convention: `README.md` is for humans; each `<name>.md` in the repo is a skill loaded by AI agents.

## Contents

- [Highlights](#highlights)
- [Quick Start](#quick-start)
- [Supported Workflows](#supported-workflows)
- [Data & Support](#data--support)
- [Environment Dependencies](#environment-dependencies)
- [Directory Structure](#directory-structure)
- [Running Guide](#running-guide)
- [Modification Guide](#modification-guide)
- [Roadmap](#roadmap)

## Highlights

### Full parameter provenance, reproducible results

All run parameters are versioned in `output/<workflow>/raw.json`, input → output dependencies are described by Snakemake, and intermediate results are automatically recoverable. Every run's parameters and versions can be traced completely afterwards.

### Schema validation plus runtime overrides

All tool parameters are progressively disclosed in `config/<workflow>.json`, human-readable and directly editable, with legal values constrained by `schema.json` (invalid values are flagged on the spot, list-type parameters are recognized automatically). Prefer not to edit files? Override nested fields at runtime from the CLI:

```bash
python workflow/Omics/run.py -m data/meta/fastq -w CLIP -o output \
  --Params.trim_galore.quality 10
```

Image output format is configurable too: png by default, pdf for publication figures, or several formats at once.

### One entry point, module assembly

A dozen-plus analysis types converge on one entry point and share `modules/` and `config/`. New workflows build on existing modules instead of reinventing the wheel.

### Run whole pipelines or single modules

Run the entire pipeline in one shot with `run.py`, or just one module: `run_module.py` drives any module standalone from its own json input, with `--until`, `--target`, `--rules`, `--samples` to scope exactly what runs, and the same `--dry-run` / `--sdm` runtime controls:

```bash
# full pipeline
python workflow/Omics/run.py -m data/meta.tsv -w RNAseq -o output

# a single module
python workflow/Omics/run_module.py -m macs3 --config my_macs3.json -o out/
```

### Add branches without touching the mainline

Adding a branch only declares new outputs, never modifies existing inputs. The DAG describes dependencies rather than imperative execution, so the mainline never drifts. Workflow changes are traceable via git; rerun decisions are cached by Snakemake.

### Precise reruns

`--until <rule>` stops at a chosen stage, `--forcerun` force-reruns any rule (with `rule:wildcard=value` down to exact wildcard combinations), and `--target-jobs` can target any intermediate product as an endpoint:

```bash
# Stop once deduped BAMs are done
python workflow/Omics/run.py -m data/meta.tsv -w Mutation -o output \
  --until mutation_markduplicates

# Precisely rerun a rule with wildcards
python workflow/Omics/run.py -m data/meta.tsv -w scRNAseq -o output \
  --counters scTE cellranger --aligner cellranger \
  --forcerun scanpy_auto:tissue=Uterus,counter=scTE
```

### No manual software installation

Software required by each workflow needs no manual installation: in conda mode Snakemake creates environments automatically from `modules/*/*.yaml`; in container mode pre-built SIF images are used directly, with bind paths derived automatically by `run.py`.

### Controllable resource allocation

The `resources` field controls threads and memory per rule (`resources.default` as fallback, overridden per rule name), preventing host memory exhaustion or CPU saturation.

### Built-in test mode

`--test <workflow>` (or `--test all`) runs config parsing, `raw.json` generation, and a dry-run on bundled test data without real execution. A few seconds to verify a workflow works.

### AI-friendly, with an LLM auto mode

All information is persisted in JSON (configs, `raw.json`, schemas) and log files, so the full context of every run is available without digging through code. Combined with dedicated skills (`<dirname>.md`), AI can diagnose failures, improve performance, and add analysis modules. `scRNAseq` goes further with an `auto` mode: built-in LLM decisions pick the number of PCs, tune clustering resolution, filter low-quality clusters, run LLM cell-type annotation, and produce an audit report. Just provide the data and tissue name:

```bash
python workflow/Omics/run.py -m data/meta/fastq -w scRNAseq -o output \
  --auto --tissue ovaries
```

## Quick Start

1. Install the base environment: snakemake >= 9.16.3 (plus apptainer for container mode, see [Environment Dependencies](#environment-dependencies)).

   ```bash
   conda create -n smk -c conda-forge -c bioconda "snakemake>=9.16.3" -y
   conda activate smk
   ```

2. Prepare input: pass a metadata file (sample, species, library layout, etc.) or a fastq directory, from which `run.py` parses sample information itself.
3. Pick an environment mode: conda mode takes `--conda-prefix` (`run.py` adds `--use-conda` for you); container mode takes `--sdm` (runs via SIF images, bind paths derived automatically).
4. Launch:

   ```bash
   # conda mode
   python workflow/Omics/run.py \
     -m data/meta/fastq \
     -w CLIP \
     -o output \
     -t 48 \
     --log log/CLIP.log \
     --conda-prefix /path/to/conda/envs \
     --Params.trim_galore.quality 10

   # container mode
   python workflow/Omics/run.py \
     -m data/meta/fastq \
     -w CLIP \
     -o output \
     -t 48 \
     --sdm \
     --rerun-triggers mtime
   ```

Add `--dry-run` first to review the jobs before letting it run. To sanity-check the workflow itself, `--test CLIP` does a side-effect-free pass on bundled test data (see [Running Guide](#running-guide)).

## Supported Workflows

`run.py` selects a workflow via `-w/--workflow_name`:

| Category | Workflow | Description | Typical Outputs |
| --- | --- | --- | --- |
| Transcriptomics | `RNAseq` | Bulk transcriptomics (quantification, differential analysis, fusion genes, variant calling) | count matrix, TE expression, DESeq2 differential analysis, GO/KEGG/GSEA enrichment, fusion genes (arriba), germline SNV/INDEL, transcript assembly (StringTie), HTML report |
| Transcriptomics | `ncRNAseq` | Non-coding RNA analysis (incl. tailer and small RNA) | demultiplexer, trim_galore, hisat2/star alignment, gatk dedup, small RNA extraction, star 3pass, tailer analysis, ncRNA expression |
| Transcriptomics | `LRtranscriptome` | Long-read transcriptome analysis | minimap2 alignment, quantification (skippable) |
| Transcriptomics | `CoCulture` | Co-culture sample analysis (species separation + TE expression) | SOAPnuke QC, hisat2 alignment, ngs_disambiguate species separation, TEtranscripts quantification |
| Epigenetics / modification | `MERIP` | MeRIP-seq / m6A-seq analysis | deduped BAM, peak results, IGV visualization |
| Epigenetics / modification | `CLIP` | iCLIP / CLIP-seq analysis | QC, alignment, PureCLIP, bedGraph / bigWig, IGV pages |
| Epigenetics / modification | `tRNAseq` | tRNA modification-induced mismatch sequencing (mim-tRNAseq) | coverage, modification quantification, CCA analysis, DESeq2 differential expression |
| Epigenetics / modification | `Fiberseq` | Fiber-seq epigenetic analysis | methylation detection, nucleosome positioning, FIRE active regions, data extraction, SV comparison analysis (exp-specific SV, enrichment, circos, OncoPrint) |
| Epigenetics / modification | `PeakCalling` | ChIP-seq / ChIRP-seq / DIP-seq peak calling | fastqc, trim_galore, bowtie2 alignment, gatk dedup, IGV/wig visualization, MACS3 peak calling, FRiP score, deeptools heatmaps, Homer annotation, peak-TE overlap, summary report |
| Epigenetics / modification | `RIPseq` | RIP-seq RNA-protein binding site analysis (IP vs Input) | fastqc, trim_galore, hisat2 alignment, gatk dedup, IGV/wig visualization, MACS3 peak calling, FRiP score, deeptools heatmaps, Homer annotation, peak-TE overlap, summary report (PPT + Excel) |
| Variants / genomes | `Mutation` | Somatic/germline mutation analysis (tumor vs normal) and cfDNA analysis | fastqc, cutadapt, bwa-mem2, GATK (Mutect2 somatic + germline), mutation spectra, cfDNA fragment length, manta SV, cnvkit CNV |
| Variants / genomes | `PacVar` | PacBio long-read variant calling | SV VCF, SNP VCF, phasing results, telomere length (4 methods), centromere stats, SV comparison analysis (exp-specific SV, enrichment, circos, OncoPrint), optional eccDNA reconstruction (DECOIL) and visualization (stat plots + cross-sample comparison + interactive IGV HTML) |
| Variants / genomes | `Population_genomics` † | Population genomics analysis | GATK HaplotypeCaller + joint genotyping, bcftools filtering, plink2/vcftools statistics, admixture population structure, easySFS history inference, PopLDdecay LD decay |
| Single-cell / spatial | `scRNAseq` | Single-cell RNA-seq analysis | Cell Ranger/STARsolo alignment, scTE TE quantification, Scanpy QC/clustering/annotation/advanced analysis |
| Single-cell / spatial | `spatial_transcriptomics` † | Spatial transcriptomics analysis | Scanpy QC, clustering, advanced analysis (Squidpy spatial neighbors, Moran's I) |
| Interaction | `KARRseq` | Kethoxal-Assisted RNA-RNA interaction sequencing | RNA-RNA interaction pairs file |
| Proteomics | `QuantMS` | Quantitative proteomics (TMT/LFQ/DIA) | mzTab quantification, MSstats statistics |

† `Population_genomics` and `spatial_transcriptomics` are not yet wired into `run.py -w`; run them directly via `snakemake -s subworkflow/<name>.smk --configfile <raw.json>`.

Detailed workflow descriptions: [subworkflow](./subworkflow/README.md)

## Data & Support

The reference databases the workflows rely on (genomes, indices, germline resources, etc.)
live under the `Database/` directory and are **not distributed with this repository**.
To obtain the database (species reference genomes, MGP/gnomAD germline variation
resources, prebuilt SIF container images, and more), contact:
**[goodzyx36@gmail.com](mailto:goodzyx36@gmail.com)**

For more bioinformatics content, visit:

- **Blog**: <https://blog.shiganluo.top/home> — analysis notes and pipeline walkthroughs
- **Bioinformatics cloud platform**: <https://bio.shiganluo.top/> — run analysis pipelines online

## Environment Dependencies

The base environment needs **snakemake** (>= 9.16.3) and **apptainer**. The specific software used by each workflow does not need manual installation —
in conda mode Snakemake creates environments automatically from `modules/*/*.yaml`;
in container mode pre-built SIF images are used.

### snakemake (conda install)

```bash
conda create -n smk -c conda-forge -c bioconda "snakemake>=9.16.3" -y
conda activate smk

snakemake --version   # expect >= 9.16.3
```

### apptainer

apptainer is unnecessary if you only use conda mode (`--conda-prefix`); `--sdm` container mode
requires the `apptainer` (or `singularity`) command. For other installation methods see the official docs:
<https://github.com/apptainer/apptainer/blob/main/INSTALL.md>.

The official **deb package** (setuid build) is recommended; it depends on `fakeroot` and `uidmap`:

```bash
# Dependencies
sudo apt install fakeroot uidmap

# Download debs from https://github.com/apptainer/apptainer/releases (not the trixie+ suffixed ones)
# Both apptainer and apptainer-suid must be installed: the suid package provides the
# setuid-root starter-suid, the key file for setuid mode; without it you fall back to
# unprivileged mode, which cannot start containers on Ubuntu 24.04
sudo dpkg -i ./apptainer_1.5.4_amd64.deb ./apptainer-suid_1.5.4_amd64.deb
```

- Ubuntu 24.04 requires the setuid build; the unprivileged build (including conda installs) fails
  to start containers due to AppArmor restrictions (see [docs/apptainer.md](docs/apptainer.md) Troubleshooting).
- Watch PATH precedence after installation: if apptainer was also installed in a conda env,
  remove it to make sure the system installation is used.

Verified working combination: Python 3.12 + Snakemake 9.22.0 + Apptainer 1.5.4 (setuid build).

Full container-mode documentation — SIF builds, path resolution, bind derivation, manual Snakemake invocation, troubleshooting, and the `apptainer exec` toolbox — lives in **[docs/apptainer.md](docs/apptainer.md)**.

## Directory Structure

- `run.py`: unified entry script; reads metadata, merges config, and invokes Snakemake.
- `run.sh`: archive of previously run commands (all commented out by default); uncomment to reuse.
- `run_module.py`: single-module entry point, used with each module's own json input.
- `node.py`: result registration and config dispatch per workflow.
- `config/`: template config files per workflow plus schema constraints.
- `subworkflow/`: Snakemake main entry for each analysis workflow.
- `modules/`: reusable module definitions (86 of them); conventions in [modules/modules.md](modules/modules.md).
- `src/`: shared utilities and standalone analysis scripts (e.g. `common/util/EnvUtil.py`).
- `docs/`: design documents and topic guides.
- `assests/`: bundled test data and test meta templates.

## Running Guide

### Run examples

`run.sh` in the repo records commands that have actually been run; reference it or uncomment what you need:

```bash
bash workflow/Omics/run.sh
```

To call Snakemake manually, mirror the command form `run.py` assembles.

### Running specific steps

`run.py` provides `--forcerun` and `--target-jobs` for precisely rerunning specific rules.

#### `--forcerun`: force-rerun rules

Accepts rule names, file paths, or `rule:wildcard=value` format. Rule names are auto-prefixed.

```bash
# Rerun a rule (all wildcards)
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun function_gsea

# Rerun the rule producing a file path
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun /path/to/output/some_file.txt

# Rerun multiple rules at once
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun trimming_Paired function_gsea
```

Rule name auto-completion: rules are renamed in subworkflows via `use rule ... as <WorkflowName>_...`. `run.py` auto-adds the prefix, so users only write the original rule name:

| User input | Converted to |
| --- | --- |
| `function_gsea` | `RNAseq_function_gsea` |
| `RNAseq_function_gsea` | `RNAseq_function_gsea` (already prefixed, unchanged) |
| `all` | `all` (special rule name, unchanged) |
| `/path/to/file` | `/path/to/file` (file path, no prefix) |

#### `--forcerun` precise reruns (with wildcards)

For rules with wildcards (e.g. `scanpy_auto` has `{tissue}` and `{counter}`), snakemake's `--forcerun` does not support bare rule names. Use the `rule:wildcard=value` format; `run.py` converts it to `--target-jobs RULE:wildcard=value --force`:

```bash
# Rerun the scTE scanpy_auto for Uterus tissue
python workflow/Omics/run.py \
  -m data/meta.tsv -w scRNAseq -o output --sdm \
  --counters scTE cellranger --aligner cellranger \
  --forcerun scanpy_auto:tissue=Uterus,counter=scTE

# Rerun multiple wildcard combinations at once
python workflow/Omics/run.py \
  -m data/meta.tsv -w scRNAseq -o output --sdm \
  --counters scTE cellranger --aligner cellranger \
  --forcerun \
    scanpy_auto:tissue=Uterus,counter=scTE \
    scanpy_auto:tissue=Uterus,counter=cellranger
```

#### `--touch`: update output timestamps

After rerunning a rule, `--touch` refreshes downstream output timestamps so the next execution does not rerun downstream rules again:

```bash
# Mark all outputs as up to date (no actual execution)
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --touch
```

`--touch` and `--dry-run` cannot be used together; when both are given, a warning is printed and `--touch` is ignored.

#### `--list-target-rules`: list runnable top-level rules

Lists rules that can be used as targets (top-level rules without wildcards). Useful for quickly seeing a workflow's entry points.

```bash
python workflow/Omics/run.py \
  -m data/meta.tsv -w scRNAseq -o output --sdm \
  --list-target-rules
```

Example output:

```
all
scRNAseq_cellranger_ref
scRNAseq_scTE_build_index
```

#### `--list-rules`: list all rules

Lists all rules in the workflow (including wildcarded intermediate rules), fully showing the rules available in the DAG.

```bash
python workflow/Omics/run.py \
  -m data/meta.tsv -w scRNAseq -o output --sdm \
  --list-rules
```

Example output (more complete than `--list-target-rules`, includes wildcarded rules):

```
all
scRNAseq_cellranger_ref
scRNAseq_scTE_build_index
scRNAseq_scanpy_auto
scRNAseq_scanpy_cluster
...
```

**`--list-target-rules` vs `--list-rules`**:

- `--list-target-rules`: only shows rules directly usable as `snakemake` targets (no wildcards); good for daily inspection.
- `--list-rules`: shows all rules (including wildcarded ones); good for debugging and understanding the full DAG.

#### Viewing the DAG and debugging

```bash
# View the jobs that would run (combine with --dry-run)
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun function_gsea --dry-run

# Pass through extra snakemake arguments
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun function_gsea \
  --snakemake-args --printshellcmds

# Generate a DAG plot
snakemake -s workflow/Omics/subworkflow/RNAseq.smk \
  --configfile output/RNAseq/raw.json \
  --cores 1 --sdm apptainer \
  --until RNAseq_function_gsea --forcerun RNAseq_function_gsea \
  --dag | dot -Tpng > dag.png
```

#### Skipping specific rules

To skip a rule and its downstream, pass snakemake's `--omit-from` through `--snakemake-args`:

```bash
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --snakemake-args --omit-from RNAseq_generate_report
```

### Test mode (`--test`)

`run.py` has a built-in test mode to quickly check whether a workflow can correctly parse config, generate `raw.json`, and pass a Snakemake dry-run.

```bash
# Test a single workflow
python workflow/Omics/run.py \
  --test ncRNAseq \
  -o output

# Test all workflows
python workflow/Omics/run.py \
  --test all \
  -o output
```

Test mode automatically:

- uses `workflow/Omics/assests/test/meta_<workflow>.tsv`
- generates test paths from `workflow/Omics/assests/test/data/`
- writes output to `output/test/` (or `test/` under your chosen output directory)
- enables `--dry-run`
- uses a local `.conda` prefix to avoid polluting shared environments

If a workflow has no corresponding test meta file, `run.py` skips it with a warning.

### `run.py` Parameters

- `-m, --meta`: metadata file or fastq directory.
- `-w, --workflow_name`: workflow name; one of `CoCulture`, `MERIP`, `RNAseq`, `ncRNAseq`, `scRNAseq`, `tRNAseq`, `CLIP`, `Mutation`, `PacVar`, `KARRseq`, `PeakCalling`, `QuantMS`, `Fiberseq`, `LRtranscriptome` (`Population_genomics` and `spatial_transcriptomics`: see the Supported Workflows footnote).
- `-o, --output_dir`: output directory.
- `-t, --threads`: Snakemake thread count.
- `--dry-run`: generate the plan only, no execution.
- `--log`: log file path.
- `--conda-prefix`: conda package cache directory (required in conda mode, omitted in container mode).
- `--sdm`: use the Apptainer container backend (SIF images). When set, `run.py` does not add `--use-conda` and instead runs via SIF images resolved from `container:` directives. Bind paths are auto-generated from the config JSON as `--singularity-args --bind ...`.
- `--singularity-args`: manually specify singularity/apptainer arguments, e.g. `--singularity-args '--bind /path1,/path2'` (used with `--sdm`, overrides the auto-generated bind list).
- `--rerun-trigger` (alias `--rerun-triggers`): snakemake rerun triggers, defaulting to all five (same as snakemake's default). When specified, only the listed triggers are checked:

  | Option | Meaning |
  | --- | --- |
  | `code` | rerun when rule definition code (.smk files) changed |
  | `input` | rerun when input file contents (hashes) changed. **Requires historical hashes in `.snakemake/metadata/`; if metadata is empty nothing can be compared and all rules rerun** |
  | `mtime` | rerun when input mtime is newer than output |
  | `params` | rerun when rule params changed |
  | `software-env` | rerun when the conda environment changed |

  For example, `--rerun-triggers mtime` reruns only when inputs are newer than outputs, checking neither code/input/params/software-env; `--rerun-trigger input` looks only at content hashes and depends on metadata.

- `--conda-frontend`: `conda` or `mamba` (default mamba).
- `--forcerun`: force-rerun given rules or file paths; format: `RULE`, `/path/to/file`, or `RULE:wildcard=value`. Rule names get the workflow prefix automatically.
- `--target-jobs`: locate jobs by wildcard constraints; format: `RULE:WILDCARD1=VALUE,...`. Used with `--forcerun`.
- `--touch`: update output file timestamps without running. Conflicts with `--dry-run`.
- `--list-target-rules`: list rules usable as targets (no wildcards), no execution.
- `--list-rules`: list all rules (including wildcarded ones), no execution.
- `--no-schema-validate`: disable schema-aware extra_args type coercion (enabled by default).
- `--snakemake-args`: extra arguments passed through to Snakemake, placed after this flag, e.g. `--snakemake-args --keep-going --rerun-incomplete`.

`run.py` also supports extra arguments passed through to the config file:

- `--key value`
- `--key=value`
- Nested fields: `--Params.trim_galore.quality 10`

### Running a single module (`run_module.py`)

A module (`modules/<name>/`) is an independent unit: it has its own Snakefile and json input, unrelated to any workflow's meta.tsv. `run_module.py` parses no meta and does no workflow dispatch — it runs one module directly, with the module's own Snakefile as the entrypoint and its json as the sole configuration. Environment handling (conda/SIF resolution, bind-path collection, cross-server Database path migration, memory detection) reuses the same logic as `run.py`:

```bash
# Run cutadapt with the module's own json (smk/json auto-discovered,
# polygenomes/ variant preferred)
python workflow/Omics/run_module.py -m cutadapt \
  -i /path/to/fastq_dir -o /path/to/cutadapt_run \
  --samples SampleA SampleB --sdm

# Use your own config and override parameters at runtime
python workflow/Omics/run_module.py -m cutadapt \
  --config my_cutadapt.json \
  -o /path/to/out --samples SampleA SampleB \
  --Params.cutadapt.quality 30 --sdm

# Stop after a given rule / build a single target
python workflow/Omics/run_module.py -m star --samples SampleA --sdm --until star_align
python workflow/Omics/run_module.py -m macs3 --samples SampleA --sdm \
  --target /path/to/out/SampleA/SampleA_peaks.narrowPeak

# List module rules; test one module / all modules (dry-run validation only)
python workflow/Omics/run_module.py -m cutadapt --list-rules
python workflow/Omics/run_module.py --test cutadapt -o /path/to/module_test
python workflow/Omics/run_module.py --test all -o /path/to/module_test
```

Input convention: module JSONs describe inputs relative to `indir` (e.g. cutadapt expects `indir/<sample>/<sample>_1.fq.gz`); point `indir` at your data directory and `outdir` at the run directory (via `-i` / `-o`, editing the json, or dotted overrides).

Target selection: most module Snakefiles have no `rule all` (their first rule carries wildcards, which snakemake rejects as a default target). Without `--until` / `--target`, the script expands each rule's output paths (substituting `--samples` values for `{sample_id}`) and generates a wrapper Snakefile (`<outdir>/.run_module/all.smk`, adding only a wildcard-free `rule all` — the module itself is never modified) before invoking Snakemake. Per-sample rules therefore require `--samples`, otherwise the script aborts with "No output patterns generated".

The effective config is written to `<outdir>/<module>.config.json` (fully traceable); the module's original json in the repo stays untouched.

#### `run_module.py` Parameters

- `-m, --module`: module name; the smk/json pair is auto-discovered under `modules/<name>/` (the `polygenomes/` variant preferred, falling back to other resolvable entry smks). Optional when `--snakefile` is given (inferred from the path).
- `--config`: module config JSON; defaults to the sibling `<stem>.json` next to the snakefile.
- `--snakefile`: use this exact module Snakefile (e.g. `modules/star/star.smk` for the single-genome variant instead of the default `polygenomes/` one); its sibling json becomes the config.
- `-o, --outdir`: run/output directory (overrides config outdir; default `./<module>_run`).
- `-i, --indir`: input directory (overrides config indir).
- `-t, --threads`: thread count (default 8).
- `--samples`: `sample_id` wildcard values (space-separated); expands per-sample rules into concrete targets so a module without `rule all` runs directly.
- `--rules`: restrict target expansion to these rules (module-local names); default: every wildcarded per-sample rule.
- `--genome`: restrict to this `genome` wildcard value (multi-genome modules).
- `--until`: stop after this rule (upstream is built automatically).
- `--target`: build this output file/rule only (no wrapper Snakefile is generated in this mode; the module smk runs directly).
- `--list-rules`: list the module's rules and exit.
- `--dry-run`, `--sdm`, `--conda-prefix`, `--conda-frontend`, `--singularity-args`, `--rerun-trigger`, `--snakemake-args`: same meanings as in `run.py`. Note `--conda-prefix` is mandatory without `--sdm` (the script errors out explicitly).
- `--test [MODULE]`: module test mode. Dry-run validates every standalone module (or one named module): config loading + schema validation, target expansion, and DAG construction. Modules carry no fixed test inputs, so DAG stops on missing inputs count as passes; parse errors, wiring errors, and wrong config key names are failures. Exits non-zero on any failure.
- `--log`: log file path.
- Dotted overrides: same as `run.py` — `--key value`, `--key=value`, `--Params.cutadapt.quality 30`; comma-separated values become lists (e.g. `--genome_paired_samples.GRCm39=S1,S2`).

## Modification Guide

### Conventions

#### Input conventions

- To ensure wildcards are generated correctly, `sample_id` must not contain `.` to avoid misidentification.
- Single-end files: usually recognized as a single `fq.gz` / `fastq.gz` file.
- Paired-end files: usually recognized as pairs `*_1.fq.gz` / `*_2.fq.gz`, or `*_R1.fq.gz` / `*_R2.fq.gz`.
- `trim_galore` is only a wrapper command; `cutadapt` is still required at runtime.
- Species names containing spaces are normalized internally in the config, e.g. `Mus musculus` -> `Mus_musculus`.

#### Output conventions

Each run generates results under `output/<workflow_name>/`, along with `raw.json` and a `log/` directory.

### General modification workflow

When adding a new workflow, also update:

- `config/<workflow>.json` and `config/<workflow>.schema.json` (overall input and constraints)
- `subworkflow/<workflow>.smk` (the workflow body)
- the json, yaml, and smk files in `modules` (functionality and constraints)
- the dispatch logic in `node.py` (result registration)

Conventions for adding functionality: [modules/modules.md](modules/modules.md)

## Roadmap

- [x] Integrate all previously used workflows
- [x] Wrap actual execution in shell for HPC compatibility
- [x] Add project skill docs
- [x] Improve meta design
- [x] Add JSON value validation module
- [x] Schema-aware extra_args type coercion
- [x] Dual shell logging for all pipelines: keep the final executed shell script in the results directory, and log every executed script under log/cmd; currently supported for scRNAseq, Fiberseq
- [ ] Multi-genome parallel runs; supported: RNAseq
- [ ] Balanced resource consumption (no host memory exhaustion or CPU saturation); initially supported: Fiberseq, RNAseq; next step is tying data volume to compute resources
- [ ] Multiple selectable image output formats
