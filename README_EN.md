# Omics Workflow

English | [中文](README.md)

This is the unified entry point for all `Omics` workflows. `run.py` in this directory selects the corresponding Snakemake subworkflow based on sample metadata and `workflow_name`, and automatically generates the `raw.json` config file each workflow needs.

- `README.md` / `README_EN.md` are documentation for humans
- `<dirname>.md` files are skills loaded by agents

## Advantages

Each advantage below comes with concrete usage — CLI commands, parameters, or workflow entry points; see "Quick Start" and "Running Guide" for full examples.

**Full parameter provenance, reproducible results**
All run parameters are versioned in `output/<workflow>/raw.json`, input → output dependencies are described by Snakemake, and intermediate results are automatically recoverable — every run's parameters and versions can be fully traced afterwards.

**Highly customizable configuration**
All tool parameters are progressively disclosed in `config/<workflow>.json`, human-readable and directly editable; or skip editing files entirely and override nested fields at runtime from the CLI (e.g. `--Params.trim_galore.quality 10`):

```bash
python workflow/Omics/run.py -m data/meta/fastq -w CLIP -o output \
  --Params.trim_galore.quality 10
```

**One entry point, highly reusable**
A dozen-plus analysis types (RNA-seq, scRNA-seq, ChIP-seq, WGS, population genomics, etc.) converge on a single entry `python run.py -w <workflow>`, sharing `modules/` and `config/`. New workflows build on existing modules instead of reinventing the wheel.

**Add branches without touching the mainline**
Adding a branch only declares new outputs, never modifies existing inputs — the DAG describes dependencies rather than imperative execution, so the mainline never drifts. Workflow changes are traceable via git; rerun decisions are cached by Snakemake.

**Exploratory analysis supported (LLM auto mode)**
The `scRNAseq` workflow offers an `auto` mode — built-in LLM decisions that automatically pick the number of PCs, tune clustering resolution, filter low-quality clusters, call an LLM for cell-type annotation, and produce an audit report. Just provide the data and tissue name:

```bash
python workflow/Omics/run.py -m data/meta/fastq -w scRNAseq -o output \
  --auto --tissue ovaries
```

**Precise reruns**
Output-driven DAG — `--until <rule>` stops at a chosen stage, `--forcerun` force-reruns any rule (supporting `rule:wildcard=value` down to exact wildcard combinations), and `--target-jobs` can target any intermediate product as an endpoint:

```bash
# Stop once deduped BAMs are done
python workflow/Omics/run.py -m data/meta.tsv -w Mutation -o output \
  --until mutation_markduplicates

# Rerun a single rule (no wildcards)
python workflow/Omics/run.py -m data/meta.tsv -w RNAseq -o output \
  --forcerun function_gsea

# Precisely rerun a rule with wildcards
python workflow/Omics/run.py -m data/meta.tsv -w scRNAseq -o output \
  --counters scTE cellranger --aligner cellranger \
  --forcerun scanpy_auto:tissue=Uterus,counter=scTE
```

**Controllable resource allocation (currently Fiberseq only)**
The `resources` field controls threads and memory per rule (`resources.default` as fallback, overridden per rule name), preventing host memory exhaustion or CPU saturation.

**No manual software installation**
Software required by each workflow needs no manual installation: in conda mode Snakemake creates environments automatically from `modules/*/*.yaml`; in container mode pre-built SIF images are used directly, with bind paths derived automatically by `run.py`.

**Built-in test mode**
`--test <workflow>` (or `--test all`) runs config parsing, `raw.json` generation, and a dry-run on bundled test data without real execution — a fast way to verify a workflow works.

**AI-friendly**
Thanks to the unified entry point and structured design, all information is persisted in JSON (configs, `raw.json`, schemas) and log files — the full context of every run is available without digging through code. Combined with dedicated skills (`<dirname>.md`), AI can directly analyze the workflows: diagnosing failures, improving performance, adding analysis modules, and more.

## Directory Structure

- `run.py`: unified entry script; reads metadata, merges config, and invokes Snakemake.
- `run.sh`: one-click run examples commonly used in this project.
- `config/`: template config files per workflow.
- `subworkflow/`: Snakemake main entry for each analysis workflow.
- `modules/`: reusable module definitions.
- `example/`: example configs and example workflow files.

## Supported Workflows

`run.py` selects a workflow via `-w/--workflow_name`:

| Workflow | Description | Typical Outputs |
| --- | --- | --- |
| `CoCulture` | Co-culture sample analysis (species separation + TE expression) | SOAPnuke QC, hisat2 alignment, ngs_disambiguate species separation, TEtranscripts quantification |
| `MERIP` | MeRIP-seq / m6A-seq analysis | deduped BAM, peak results, IGV visualization |
| `RNAseq` | Bulk transcriptomics (quantification, differential analysis, fusion genes, variant calling) | count matrix, TE expression, DESeq2 differential analysis, GO/KEGG/GSEA enrichment, fusion genes (arriba), germline SNV/INDEL, transcript assembly (StringTie), HTML report |
| `ncRNAseq` | Non-coding RNA analysis (incl. tailer and small RNA) | demultiplexer, trim_galore, hisat2/star alignment, gatk dedup, small RNA extraction, star 3pass, tailer analysis, ncRNA expression |
| `CLIP` | iCLIP / CLIP-seq analysis | QC, alignment, PureCLIP, bedGraph / bigWig, IGV pages |
| `Mutation` | Somatic/germline mutation analysis (tumor vs normal) and cfDNA analysis | fastqc, cutadapt, bwa-mem2, GATK (Mutect2 somatic + germline), mutation spectra, cfDNA fragment length, manta SV, cnvkit CNV |
| `PacVar` | PacBio long-read variant calling | SV VCF, SNP VCF, phasing results, telomere length (4 methods), centromere stats, optional eccDNA reconstruction (DECOIL) and visualization (stat plots + cross-sample comparison + interactive IGV HTML) |
| `KARRseq` | Kethoxal-Assisted RNA-RNA interaction sequencing | RNA-RNA interaction pairs file |
| `PeakCalling` | ChIP-seq / ChIRP-seq / DIP-seq peak calling | fastqc, trim_galore, bowtie2 alignment, gatk dedup, IGV/wig visualization, MACS3 peak calling, FRiP score, deeptools heatmaps, Homer annotation, peak-TE overlap, summary report |
| `QuantMS` | Quantitative proteomics (TMT/LFQ/DIA) | mzTab quantification, MSstats statistics |
| `tRNAseq` | tRNA modification-induced mismatch sequencing (mim-tRNAseq) | coverage, modification quantification, CCA analysis, DESeq2 differential expression |
| `scRNAseq` | Single-cell RNA-seq analysis | Cell Ranger/STARsolo alignment, scTE TE quantification, Scanpy QC/clustering/annotation/advanced analysis |
| `spatial_transcriptomics` | Spatial transcriptomics analysis | Scanpy QC, clustering, advanced analysis |
| `Population_genomics` | Population genomics analysis | GATK HaplotypeCaller + joint genotyping, bcftools filtering, plink2/vcftools statistics, admixture population structure, easySFS history inference, PopLDdecay LD decay |
| `Fiberseq` | Fiber-seq epigenetic analysis | methylation detection, nucleosome positioning, FIRE active regions, data extraction |

Detailed workflow descriptions: [subworkflow](./subworkflow/README.md)


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
# the unprivileged mode, which cannot start containers on Ubuntu 24.04
sudo apt dpkg -i ./apptainer_1.5.4_amd64.deb ./apptainer-suid_1.5.4_amd64.deb
```

- Ubuntu 24.04 requires the setuid build; the unprivileged build (including conda installs) fails
  to start containers due to AppArmor restrictions (see "Troubleshooting").
- Watch PATH precedence after installation: if apptainer was also installed in a conda env,
  remove it to make sure the system installation is used.

Verified working combination: Python 3.12 + Snakemake 9.22.0 + Apptainer 1.5.4 (setuid build).

## Quick Start

1. Prepare input data.
   - If passing a metadata file, it should usually contain sample, species, and library layout information.
   - If passing a fastq directory, `run.py` parses sample information directly from the directory.
2. Prepare the runtime environment: choose conda mode or container (Apptainer/SIF) mode.
   - **conda mode**: specify `--conda-prefix`; `run.py` automatically adds `--use-conda`.
   - **container mode**: specify `--sdm`; `run.py` omits `--use-conda` and lets Snakemake run jobs via SIF images resolved from `container:` directives. SIF images are generated from conda yamls by `EnvUtil.py` (see `src/common/util/EnvUtil.py`); the `sif()` function in `common.smk` looks up paths by yaml file stem.
3. Launch the workflow with `run.py`.

Recommended example (conda mode):

```bash
python workflow/Omics/run.py \
  -m data/meta/fastq \
  -w CLIP \
  -o output \
  -t 48 \
  --log log/CLIP.log \
  --conda-prefix /data/pub/zhousha/env/mutation_0.1 \
  --Params.trim_galore.quality 10
```

Container mode example:

```bash
python workflow/Omics/run.py \
  -m data/meta/fastq \
  -w CLIP \
  -o output \
  -t 48 \
  --sdm \
  --rerun-triggers mtime
```

In container mode `run.py` automatically extracts paths from the config JSON and generates
`--singularity-args --bind ...`; no manual specification needed. To override, use `--singularity-args`.

To inspect the workflow without executing it, add `--dry-run`.

### Apptainer Container Mode in Detail

#### How it works

Each module's conda environment YAML (`modules/<module>/<name>.yaml`) can be built into an Apptainer SIF image. Snakemake calls the corresponding SIF automatically at runtime via `container:` directives.

#### Building SIF images

Use `EnvUtil.py` to generate `.def` files from conda YAMLs and build `.sif` images:

```bash
# Generate .def files (no build)
python workflow/Omics/src/common/util/EnvUtil.py gen \
  --modules-dir workflow/Omics/modules \
  --outdir /path/to/env

# Build existing .def files into .sif images
python workflow/Omics/src/common/util/EnvUtil.py build \
  --modules-dir workflow/Omics/modules \
  --outdir /path/to/env

# Full pipeline: YAML → .def → .sif (recommended)
python workflow/Omics/src/common/util/EnvUtil.py all \
  --modules-dir workflow/Omics/modules \
  --outdir /path/to/env
```

Built SIF layout:

```
/path/to/env/
  bowtie2/
    bowtie2.sif
  macs3/
    macs3.sif
  fastqc/
    fastqc.sif
  ...
```

#### SIF path resolution

The `sif()` function in `common.smk` handles path resolution:

```python
# modules/common/common.smk
def sif(yaml_filename: str) -> str:
    """Look up SIF path from config["env"] mapping"""
```

The `env` field in the config JSON maps YAML stems → SIF paths:

```json
{
  "env": {
    "env_dir": "/path/to/env",
    "bowtie2": "/path/to/env/bowtie2/bowtie2.sif",
    "macs3": "/path/to/env/macs3/macs3.sif"
  }
}
```

Modules without an explicit mapping fall back to `<env_dir>/<module>/<yaml_stem>.sif`.

#### Automatic bind path derivation

In `run.py --sdm` mode, `_collect_bind_paths()` recursively scans all path strings in the config JSON, derives the list of directories to bind-mount, and merges parent directories to minimize the number of binds.

Automatic binds include:
- Input data directories (raw FASTQ paths)
- Output directories
- Reference genome directories
- Software/env directories
- `/tmp` (temporary files)

To add extra bind paths:

```bash
python workflow/Omics/run.py \
  -m data/meta.tsv -w PeakCalling -o output \
  --sdm \
  --singularity-args '--bind /extra/data,/extra/ref'
```

#### Running Snakemake + Apptainer manually

Calling Snakemake directly without `run.py`:

```bash
snakemake -s workflow/Omics/subworkflow/PeakCalling.smk \
  --configfile output/Rnp/PeakCalling/raw.json \
  --cores 48 \
  --sdm apptainer \
  --singularity-args '--bind /data,/home/user/Database,/home/user/Data,/tmp' \
  --rerun-triggers code input mtime params software-env
```

#### Troubleshooting

- **`command not found`**: the tool inside the SIF is not on PATH; check that the `%environment` or `%post` section of the `.def` file sets PATH correctly.
- **Permission errors**: make sure the `--bind` target directories are readable/writable by the current user.
- **Stale SIFs**: after a conda YAML update, rebuild the SIF (`EnvUtil.py all`).
- **`Could not write info to setgroups: Permission denied` / `Error while waiting event for user namespace mappings: no event received` on container start**:
  the unprivileged apptainer installed via conda relies on user namespaces, which Ubuntu 24.04 blocks by default
  (`kernel.apparmor_restrict_unprivileged_userns=1`), causing startup failure. Fix: switch to the setuid
  `apptainer-suid` build (see "Environment Dependencies").
  After the fix, if the unprivileged conda apptainer still precedes it on PATH, the same errors reappear.
  A simultaneous `fuse2fs not found` message is an unrelated hint and can be ignored.

#### Direct apptainer execution

Tools inside SIF images can be invoked directly via `apptainer exec` without going through Snakemake. Each SIF's `%runscript` is set to `exec "$@"`, so the passed command runs directly.

Basic syntax:

```bash
apptainer exec \
  --bind /data,/home/user/Database,/tmp \
  /path/to/env/<module>/<module>.sif \
  <tool> [arguments...]
```

Common tool examples:

```bash
BIND="--bind /home/luosg/Data,/home/luosg/Database,/tmp"
ENV="/home/luosg/Database/env"

# FastQC quality control
apptainer exec $BIND $ENV/fastqc/fastqc.sif \
  fastqc --threads 8 -o output/QC/ sample_1.fq.gz sample_2.fq.gz

# Trim Galore adapter trimming
apptainer exec $BIND $ENV/trim-galore/trim-galore.sif \
  trim_galore --quality 25 --paired -o output/trimmed/ sample_1.fq.gz sample_2.fq.gz

# Bowtie2 alignment
apptainer exec $BIND $ENV/bowtie2/bowtie2.sif \
  bowtie2 -x /path/to/index -1 sample_1.fq.gz -2 sample_2.fq.gz --very-sensitive -S output.sam

# SAMtools processing
apptainer exec $BIND $ENV/samtools/samtools.sif \
  samtools sort -@ 8 -o output.bam input.sam

# GATK MarkDuplicates
apptainer exec $BIND $ENV/gatk/gatk4.sif \
  gatk MarkDuplicates -I input.bam -O output.bam -M metrics.txt

# MACS3 peak calling
apptainer exec $BIND $ENV/macs3/macs3.sif \
  macs3 callpeak -t IP.bam -c Input.bam -f BAM -g mm -n sample --outdir peaks/

# HOMER annotation
apptainer exec $BIND $ENV/homer/homer.sif \
  annotatePeaks.pl peaks.narrowPeak mm39 -gtf annotation.gtf > annotation.txt

# BEDTools intersection
apptainer exec $BIND $ENV/bedtools/bedtools.sif \
  bedtools intersect -a peaks.bed -b TE.gtf -wa -wb > overlap.bed

# deeptools computeMatrix
apptainer exec $BIND $ENV/deeptools_heatmap/deeptools_heatmap.sif \
  computeMatrix reference-point -S signal.bigwig -R peaks.bed -a 3000 -b 3000 -o matrix.gz

# deeptools plotHeatmap
apptainer exec $BIND $ENV/deeptools_heatmap/deeptools_heatmap.sif \
  plotHeatmap -m matrix.gz -o heatmap.png --whatToShow 'heatmap and colorbar'

# FRiP score (bedtools + samtools combo)
apptainer exec $BIND $ENV/samtools/samtools.sif \
  samtools view -c -F 4 aligned.bam  # total mapped reads
apptainer exec $BIND $ENV/bedtools/bedtools.sif \
  bedtools intersect -a aligned.bam -b peaks.bed -u | samtools view -c  # reads in peaks
```

Debugging the container environment:

```bash
# Enter an interactive shell in the container
apptainer shell --bind /home/luosg/Data,/home/luosg/Database,/tmp \
  /home/luosg/Database/env/bowtie2/bowtie2.sif

# Check commands available inside the container
apptainer exec $ENV/bowtie2/bowtie2.sif which bowtie2
apptainer exec $ENV/bowtie2/bowtie2.sif bowtie2 --version
```

### Test Mode (`--test`)

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

## Running Guide

### Run examples

Existing example scripts in the repo can be found in `run.sh`. Its typical execution form is:

```bash
bash workflow/RNA-SNP/run.sh
```

If you need to call Snakemake manually, refer to the command form `run.py` finally assembles.

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

Prevents re-executing downstream rules when re-running after a rule has been rerun:

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

### `run.py` Parameters

- `-m, --meta`: metadata file or fastq directory.
- `-w, --workflow_name`: workflow name; one of `CoCulture`, `MERIP`, `RNAseq`, `ncRNAseq`, `CLIP`, `Mutation`, `PacVar`, `KARRseq`, `PeakCalling`, `QuantMS`, `tRNAseq`, `scRNAseq`, `Fiberseq`.
- `-o, --output_dir`: output directory.
- `-t, --threads`: Snakemake thread count.
- `--dry-run`: generate the plan only, no execution.
- `--log`: log file path.
- `--conda-prefix`: conda package cache directory (required in conda mode, omitted in container mode).
- `--sdm`: use the Apptainer container backend (SIF images). When set, `run.py` does not add `--use-conda` and instead runs via SIF images resolved from `container:` directives. Bind paths are auto-generated from the config JSON as `--singularity-args --bind ...`.
- `--singularity-args`: manually specify singularity/apptainer arguments, e.g. `--singularity-args '--bind /path1,/path2'` (used with `--sdm`, overrides the auto-generated bind list).
- `--rerun-trigger`: snakemake rerun trigger, default `input`. Options and meanings:

  | Option | Meaning |
  | --- | --- |
  | `code` | rerun when rule definition code (.smk files) changed |
  | `input` | rerun when input file contents (hashes) changed. **Requires historical hashes in `.snakemake/metadata/`; if metadata is empty nothing can be compared and all rules rerun** |
  | `mtime` | rerun when input mtime is newer than output |
  | `params` | rerun when rule params changed |
  | `software-env` | rerun when the conda environment changed |

  Without `--rerun-trigger`, Snakemake uses all five triggers by default. `--rerun-trigger input` means **only** input content changes are checked — no code/mtime/params/software-env checks; lighter but dependent on metadata.

- `--conda-frontend`: `conda` or `mamba`.
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

- When adding a new workflow, also update:
  - `config/<workflow>.json`
  - `config/<workflow>.schema.json`
  - `subworkflow/<workflow>.smk`
  - the json, yaml, and smk files in `modules`
  - the dispatch logic in `node.py`
- Default values for tool parameters follow the software defaults or values adapted for the workflow.
- `" ".join(cmd). cmd must not contain None.


## Roadmap

- [x] Integrate all previously used workflows
- [x] Wrap actual execution in shell for HPC compatibility
- [x] Add project skill docs
- [x] Improve meta design
- [x] Add JSON value validation module
- [x] Schema-aware extra_args type coercion
- [ ] Multi-genome parallel runs; supported: RNAseq
- [ ] Balanced resource consumption (no host memory exhaustion or CPU saturation); supported: Fiberseq
