# peak_te_overlap — 峰与 TE 注释重叠分析模块

把 MACS3 narrowPeak 与 TE 注释 GTF 求交:TE 位点视角(亚家族重叠与
IP/Input 富集)、峰视角(每峰含多少 TE、TE 覆盖比例),并出富集图。
引擎为模块私有脚本 `bin/intersect_te.py`、`bin/peak_centric_overlap.py`、
`bin/plot_te_overlap.py`(内部用 bedtools intersect)。同目录 `polygenomes/`
为多基因组版子模块(独立 .smk),本 README 覆盖 `peak_te_overlap.smk`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `peaks_indir` + 布局参数 | 形态 A 拼接:`{peaks_indir}/{sample_id}/{sample_id}_peaks.narrowPeak`(MACS3 输出) |
| `bam_indir`(可选) | IP BAM `{bam_indir}/{sample_id}/{sample_id}.sorted_markdup.bam`;空则不做读段计数 |
| `sample_ip_input_map` | IP → Input 样本;配对 Input BAM `{bam_indir}/{input}/{input}.sorted_markdup.bam` 供富集计算 |
| `genome.te_gtf` | TE 注释 GTF(平铺键),必需,parse 期校验 |
| `samples` | 样本列表(汇总图输入) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_peak_te_overlap.bed          # peak×TE 交集全表
{outdir}/{sample_id}/{sample_id}_te_subfamily_overlap.tsv     # TE 位点级重叠 + IP/Input 读段计数
{outdir}/{sample_id}/{sample_id}_peak_centric_te.tsv          # 峰级:TE 数/类型/覆盖比例
{outdir}/te_subfamily_overlap_combined.png                    # 汇总图
{outdir}/te_subfamily_overlap_combined.tsv                    # 汇总表
{outdir}/{ip}_vs_{input}_enrichment.png                       # 每 IP:Input 富集图(非 combine 模式)
{outdir}/.peak_te_overlap_fig.done                            # 汇总完成标记(touch)
```

日志:`{logdir}/{sample_id}/peak_te_overlap.log`、`peak_centric_overlap.log`、
`{logdir}/../group/peak_te_overlap/peak_te_overlap.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `outdir` / `logdir` | — | `output` / `logs` | 输出/日志目录 |
| `peaks_indir` | `--peak` | `peaks` | 峰目录 |
| `bam_indir` | `--ip-bam`/`--input-bam` | `""` | BAM 目录;空则跳过读段计数 |
| `samples` | — | `[]` | 样本列表 |
| `sample_ip_input_map` | `--ip-input-pairs IP:Input;...` | `{}` | 富集配对 |
| `genome.te_gtf` | `--te-gtf` | 必需 | TE 注释 GTF |
| `Procedure.bedtools` | `--bedtools` | `bedtools` | bedtools 可执行路径 |
| `Params.peak_te_overlap.method` | `--method` | `count` | 富集口径(`interval`/`count`/`reads`;脚本默认 `reads`,smk 固定传 `count`) |
| `Params.peak_te_overlap.sort_by` | `--sort-by` | `te_length` | 排序(`enrichment`/`te_length`) |
| `Params.peak_te_overlap.top_n` | `--top-n` | `30` | 图中展示的亚家族数 |
| `Params.peak_te_overlap.combine` | `--combine` | `false` | true 时单张合并图(此时无 per-pair 富集图) |

环境:`conda: "peak_te_overlap.yaml"`(模块目录内)。

直接调用示例(脱离流程,bin/ 脚本 CLI 与 .smk 一致):

```bash
python3 bin/intersect_te.py \
  --peak s1_peaks.narrowPeak --te-gtf TE.gtf --bedtools bedtools \
  --overlap-bed s1_peak_te_overlap.bed --subfamily-tsv s1_te_subfamily_overlap.tsv \
  --sample-id s1 --ip-bam s1.sorted_markdup.bam --input-bam input.sorted_markdup.bam
python3 bin/peak_centric_overlap.py \
  --peak s1_peaks.narrowPeak --te-gtf TE.gtf \
  --output s1_peak_centric_te.tsv --bedtools bedtools --sample-id s1
python3 bin/plot_te_overlap.py \
  --tsvs s1_te_subfamily_overlap.tsv,s2_te_subfamily_overlap.tsv \
  --output outdir --combined-tsv outdir/te_subfamily_overlap_combined.tsv \
  --method count --sort-by te_length --top-n 30 \
  --ip-input-pairs "s1:input1;s2:input2"
```

## 2. 算法原理与生物学意义

bedtools intersect 把 narrowPeak 与 TE 注释求交。两个互补视角:

- **TE 位点视角**(`intersect_te.py`):每个 TE 位点记录重叠峰数、重叠比例、
  TE class(由 GTF 属性解析),配 IP/Input BAM 的读段计数,汇总图给
  log2(IP/Input) 亚家族富集(上图平均 TE 长度)——回答"哪些 TE 亚家族
  被富集"。TE 在衰老、肿瘤中常重获开放染色质/组蛋白修饰,峰-TE 重叠富集
  是识别被调控 TE 的标准切入。
- **峰视角**(`peak_centric_overlap.py`):每峰的 TE 数、TE class 集合、
  TE 覆盖碱基数/覆盖比例——回答"哪些峰落在 TE 上",直接指导实验验证
  (挑含 TE 的峰做报告基因/CRISPR)。

结果怎么读:`*_te_subfamily_overlap.tsv` 为位点级明细(TE、class、重叠峰数、
IP/Input 读段);`*_peak_centric_te.tsv` 逐峰一行,`te_coverage_frac` 高的峰
即 TE 内峰;富集图 log2(IP/Input) > 0 的亚家族为候选被调控 TE。

## 3. 参考文献

1. bedtools — https://bedtools.readthedocs.io/ ,
   https://github.com/arq5x/bedtools2
2. Quinlan AR, Hall IM. *BEDTools: a flexible suite of utilities for
   comparing genomic features*. Bioinformatics 26(6):841–842, 2010.
   DOI: 10.1093/bioinformatics/btq033
