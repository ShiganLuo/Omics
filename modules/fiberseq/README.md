# fiberseq（Fiber-seq 下游分析家族模块，布局 C：父目录 + 子模块）

单工具链（python + numpy/scipy/matplotlib/pysam + UCSC bedGraphToBigWig），
所有子模块共享父目录 `fiberseq.yaml`（一个 SIF）。分析方法来源：
Cell Rep Methods 2024 [1]、Genome Res 2024 [2]、Nat Genet 2025 [3]。

## 子模块（每个 .smk 一个原子分析，独立可复用）

| 子模块 | 规则 | 输出（`fiberseq_analysis/<子模块>/{sample_id|contrast}/`） |
|---|---|---|
| `annotate/` | fiberseq_annotate_run | peaks_annotated.tsv / peaks_summary.tsv |
| `motif/` | fiberseq_motif_scan, fiberseq_motif_profile | motif_sites.bed / motif_profile.tsv / 谱图 |
| `footprint/` | fiberseq_footprint_run | footprint_{categories,profiles,scores}.tsv + PDF |
| `nucpos/` | fiberseq_nucpos_run | nuc_offset_{per_anchor,summary}.tsv + PDF |
| `quant/` | fiberseq_quant_run, fiberseq_quant_hap, fiberseq_quant_merge | region_quant.tsv / matrix/ |
| `diff/` | fiberseq_diff_run | diff_*.tsv + volcano + 基因表（供 modules/function） |
| `mods/` | fiberseq_mods_longtable, fiberseq_mods_sigtest | {metric}_long.tsv + {contrast}/significance_test.tsv + significant_objects.{tsv,bed}（BED 可直接喂 regionplot regions） |
| `coactuation/` | fiberseq_coactuation_run | coactuation*.tsv |
| `haplotype/` | fiberseq_haplotype_diff | haplotype_diff.tsv |
| `censat/` | fiberseq_censat_run | censat_*.tsv + PDF |
| `track/` | fiberseq_track_run | *.bw + trackDb.txt |

共享脚本在 `bin/`，共享库 `src/common/util/FiberlibUtil.py`。

## 绘图输出（formats 参数）

各绘图脚本/规则支持 `formats` 多格式同时输出（modules.md #9）：每种格式
一个文件 `<base>.<fmt>`；合法值 `png pdf svg ps eps tif tiff jpg jpeg pgf
raw rgba`（matplotlib savefig 全集），非法值在解析期即报错。dpi 统一传递
（矢量格式由 matplotlib 忽略）。默认值保持既有输出格式不变（默认 `pdf`；
`bin/utils/region_mod_compare.py` 默认 `png, svg`）。

| 参数路径 | CLI | 默认 | 图（base，位于各自输出目录） |
|---|---|---|---|
| `Params.fiberseq_motif.formats` | `--format`（可重复） | pdf | `{sample_id}.motif_profile` |
| `Params.fiberseq_footprint.formats` | `--format`（可重复） | pdf | `{sample_id}.footprint_profile` |
| `Params.fiberseq_nucpos.formats` | `-f/--format`（可重复） | pdf | `{sample_id}.nucpos` |
| `Params.fiberseq_diff.formats` | `-f/--format`（可重复） | pdf | `volcano_{contrast}` |
| `Params.fiberseq_censat.formats` | `--format`（可重复） | pdf | `{sample_id}.censat_profile` |

注：`censat_profile.py`/`footprint_profile.py`/`motif_profile.py` 的 `-f`
已被 `--flank` 占用，故这三个脚本只有长选项 `--format`。
`bin/utils/region_mod_compare.py`（两组 region 修饰比较，独立 CLI，未接入
.smk）同样接受 `-f/--format`，图 base 为 `fig_class_metrics` /
`fig_region_scatter` / `fig_profiles_accessibility` / `fig_profiles_mcpg`。

## 跨子模块输入约定

- 元素文件（indir = `common/4_elements`，唯一权威路径）：
  `{indir}/{sample_id}/{sample_id}.{m6a,nuc,msp,cpg}.bed.gz`、
  `.fire.bed.gz`（bgzip+tabix）、`.fire_peaks.bed`
- 分析产物根 = `fiberseq_analysis/`（子目录无 `fiberseq_` 前缀）：
- motif 位点：`{analysis}/motif/{sample_id}/{sample_id}.motif_sites.bed`
- 定量矩阵：`{analysis}/quant/matrix/`；hap 定量：`{analysis}/quant/{sample_id}/`
- peak 注释：`{analysis}/annotate/{sample_id}/`（diff 用，annot_sample 可指定）
- 命令回执双写：outdir 里 `<rule>_<sample>.sh`（无时间戳，每次覆盖，
  与结果同在，可直接复制重跑）；`{logdir}/cmd/{sample_id|contrast}/` 里
  `<rule>_<sample>_<timestamp>.sh`（带时间戳历史备份）

## 环境

conda/SIF 共用 `fiberseq.yaml`（SIF：`fiberseq.sif`）。子模块 smk 以
`include: "../../common/common.smk"` 直接引入共享工具（与 modules/openms 一致），
`conda: "../fiberseq.yaml"`。

## Sources

- [1] https://doi.org/10.1016/j.crmeth.2024.100911
- [2] https://doi.org/10.1101/gr.279095.124
- [3] https://doi.org/10.1038/s41588-025-02324-w
