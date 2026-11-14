# fibertools — Fiber-seq 单分子染色质分析模块(prepare + analysis)

fibertools-rs(`ft`)驱动的 Fiber-seq 分析,目录含两个 stem:

- `fibertools_prepare`:m6A 预测 → 核小体/MSP → FIRE 调用(带 SPRQ 化学
  自动检测);
- `fibertools_analysis`:元素 BED 提取、QC、FIRE peaks(FDR 控制)、
  单倍型拆分。

共享资源:`fibertools.yaml`/`.def`(conda/容器)、`Params.fibertools`
参数段、`fibertools.json`(公共 config 示例)、`bin/` 引擎脚本
(`call_peaks.py`、`shuffle_bed.py`、`distributed_plot.py`)、
`fibertools_build.sh`/`fibertools.post.sh`(容器构建)、`FIBERSEQ_ANALYSIS.md`
(分析设计记录)。目录下 `polygenomes/` 子模块不在本文档范围。

## 1. 模块输入输出与参数

### 1.1 输入

fibertools_prepare(形态 A:`indir` + 布局):

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}.bam` | PacBio HiFi CCS BAM(含动力学信息)。自动检测是否已带 m6A:① @RG 头含 SPRQ 结合试剂编号 `103-496-900`;② 前 100 reads 含 `MM:Z:` 标签。检测在容器内 shell 进行,DAG 期另有 apptainer 探针(`_has_m6a`)做 fire 输入路由 |

fibertools_analysis:

| 输入 | 说明 |
|---|---|
| `{indir}/{sample_id}/{sample_id}{bam_suffix}` | 已调用 m6a/nuc/FIRE 的比对 BAM;`bam_suffix` 默认 `.fiberseq.fire.aligned.bam`(json 示例为 `.sorted.bam`) |
| `genome.fai` | fdr 模式必需(chrom sizes,构建打乱零分布);缺失时 `_require()` 报错并指名 `genome.fai` |

### 1.2 输出

prepare(`{outdir}/{sample_id}/`):

```
{sample_id}.fiberseq.bam        # ft predict-m6a 产物(SPRQ 输入则软链直通)
{sample_id}.fiberseq.nuc.bam    # ft add-nucleosomes(核小体 + MSP)
{sample_id}.fiberseq.fire.bam   # ft fire(FIRE 调用入 Aq/Ma 标签)
```

analysis(`{outdir}/{sample_id}/`,`{hap}` ∈ hap1/hap2):

```
{sample_id}.m6a.bed.gz / .nuc.bed.gz / .msp.bed.gz / .cpg.bed.gz   # ft extract
{sample_id}.fire.bed.gz                                           # ft fire --extract(BED9)
{sample_id}.fire_peaks.bed / {sample_id}.fire_fdr_table.tsv        # call_peaks.py
{sample_id}.qc.tsv                                                 # ft qc
{sample_id}.{hap}.bam(.bam.bai)                                    # ft_split_hap
{sample_id}.{hap}.m6a/nuc/msp/cpg.bed.gz                           # ft_extract_hap
```

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.fibertools` | 可执行文件本身 | `ft` | fibertools-rs 可执行路径 |
| `Params.fibertools.ont` | `--ont` | false | ONT 数据(ft fire 相关规则) |
| `Params.fibertools.add_nucleosomes_manual` | — | null | true 时强制 fire 消费 `.fiberseq.nuc.bam`(否则按 m6A 自动检测路由) |
| `Params.fibertools.call_peaks_mode` | `--mode` | `fdr` | `fdr` / `fire_frac` 两种 FIRE 峰判定模式 |
| `Params.fibertools.sd_cov` | `--sd-cov` | 5.0 | fdr 模式:覆盖阈值(零分布标准差倍数) |
| `Params.fibertools.max_fdr` | `--max-fdr` | 0.05 | fdr 模式:FDR 上限 |
| `Params.fibertools.min_fire_frac` | `--min-fire-frac` | 0.1 | fire_frac 模式:FIRE 覆盖占比下限 |
| `Params.fibertools.haps` | `--haps` | false | 按单倍型分别出峰 |
| `Params.fibertools.acf` | `ft qc --acf` | false | QC 追加 ACF 指标 |
| `bam_suffix` | — | `.fiberseq.fire.aligned.bam` | analysis 输入 BAM 后缀 |
| `genome.fai` | `--fai` | — | fdr 模式必填 |
| `tmp_dir` | `--tmp-dir` | `{outdir}/tmp` | scratch 目录 |

`call_peaks.py` 引擎另有 `--samtools`(`samtools`)、`--seed`(42)、
`--log`、`--sample-id`;`shuffle_bed.py` 有 `-g`(chrom sizes)、
`-k/-K`(同染色体内打乱,默认开)、`-s`(seed 42)。资源默认:
predict_m6a/add_nucleosomes threads=16,其余多为 8(均经 `rule_res()` 可
覆盖)。

直接调用示例(bin/ 引擎脚本):

```bash
# FIRE peaks(FDR 模式)
python bin/call_peaks.py \
  --bam s1.fiberseq.fire.aligned.bam \
  --out-peaks s1.fire_peaks.bed --fdr-table s1.fire_fdr_table.tsv \
  --mode fdr --fai genome.fa.fai --shuffle-script bin/shuffle_bed.py \
  --sd-cov 5.0 --max-fdr 0.05 --threads 8 --ft ft \
  --tmp-dir tmp --sample-id s1

# 打乱零分布
python bin/shuffle_bed.py -i fibers.bed -g genome.fa.fai -o null.bed -s 42
```

## 2. 算法原理与生物学意义

Fiber-seq 在单分子上同时读出核小体足迹、开放染色质(MSP)、CpG 甲基化
与调控元件(FIRE)。准备链(`fibertools_prepare`):

- `ft predict-m6a`:利用 HiFi 动力学信号预测 m6A 位点(足迹外 m6A 标记
  调控蛋白触达);SPRQ 化学已带 m6A 标签时自动跳过(软链直通);
- `ft add-nucleosomes`:按 m6A 足迹分布推断核小体定位,MSP = 不被核小体
  占据的区段;
- `ft fire`:把 MSP 分类为调控元件(FIRE,精度 q 见 fiberseq/FIRE 与
  fiberbrowser/README.md 的公式说明)。

分析链(`fibertools_analysis`):`ft extract` 提取 m6a/nuc/msp/cpg 四类
BED12;`call_peaks.py` 在 **fdr 模式**下取纤维跨度、用 `shuffle_bed.py`
按染色体打乱构建零分布(固定 seed 可重复),按覆盖 sd 倍数(`sd-cov`)与
`max-fdr` 决定峰阈值,输出峰 BED 与阈值表;**fire_frac 模式**按 FIRE
覆盖占比 ≥ `min-fire-frac` 出峰;`--haps` 按单倍型分别出峰;
`ft_split_hap`/`ft_extract_hap` 按 HP 标签拆分单倍型并重提取元素。结果
阅读:fire_peaks.bed 为高置信调控元件集合(FDR 表给出逐阈值对应 FDR);
qc.tsv 汇总 m6A/核小体/FIRE 计数等健全性指标;BED 四件套是下游
fiberseq/mods、censat、track 等模块的公共输入。

## 3. 参考文献

1. fibertools-rs(`ft`)—— https://github.com/fiberseq/fibertools-rs;
   文档 https://fiberseq.github.io/
2. FIRE 分析流程 —— https://github.com/fiberseq/FIRE
   (FIRE score/precision/q 约定)
