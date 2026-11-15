# peak_te_overlap/polygenomes — 多基因组 Peak×TE 重叠分析子模块

ChIP-seq/ATAC-seq 峰(narrowPeak)与转座子(TE)注释的重叠分析
polygenomes 变体:按基因组(`{genome}` 通配符)分别分析,输出 TE 亚家族
水平重叠/富集与 peak 水平 TE 占比,并出富集图。引擎脚本共享父目录
`modules/peak_te_overlap/bin/`(intersect_te.py / peak_centric_overlap.py /
plot_te_overlap.py),conda 环境共享父目录 `../peak_te_overlap.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `peaks_indir`(config) | 峰目录(默认 `peaks`);input 函数拼 `{peaks_indir}/{sample_id}/{sample_id}_peaks.narrowPeak` |
| `bam_indir`(config) | 可选:IP/Input BAM 目录;IP BAM 拼 `{bam_indir}/{sample_id}/{sample_id}.sorted_markdup.bam` |
| `sample_ip_input_map`(config) | IP 样本 → 对照(Input)样本映射;命中者额外取 Input BAM 做 read 计数 |
| `genome.<g>.te_gtf` | 每基因组 TE 注释 GTF(`<g>` = `{genome}` 通配符取值) |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}_peak_te_overlap.bed        # peak×TE 重叠全表(bedtools intersect)
{outdir}/{genome}/{sample_id}/{sample_id}_te_subfamily_overlap.tsv   # 位点级:重叠比例/peak 数/TE class/IP Input read 计数
{outdir}/{genome}/{sample_id}/{sample_id}_peak_centric_te.tsv        # peak 级:TE 数/class 集/覆盖碱基/覆盖比例
{outdir}/{genome}/te_subfamily_overlap_combined.<fmt>                # 富集图(<fmt> 取自 formats,默认 png)
{outdir}/{genome}/te_subfamily_overlap_combined.tsv                  # 合并 TSV
{outdir}/{genome}/.peak_te_overlap_fig.done                          # 图规则标记文件
```

(separate 模式下另按 IP:Input 对各出一图
`{ip}_vs_{input}_enrichment.<fmt>`,每格式首对图软链为
`te_subfamily_overlap_combined.<fmt>`。)收尾规则
`peak_te_overlap_result` 聚合 chart + 各样本 peak_centric TSV。

### 1.3 参数

config 位于 `Params.peak_te_overlap` 下(出图规则),工具路径在
`Procedure.bedtools`:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `method` | `plot_te_overlap.py --method` | `count` | 富集计算方法 |
| `sort_by` | `--sort-by` | `te_length` | 亚家族排序依据 |
| `top_n` | `--top-n` | 30 | 图中展示的亚家族数 |
| `combine` | `--combine` | false | true 时全部 IP:Input 对合并一张图;false 每对一图 |
| `formats` | `-f/--format`(可重复) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集) |
| `sample_ip_input_map` | `--ip-input-pairs`(格式 `ip:input;...`) | {} | IP/Input 配对(见输入) |
| `Procedure.bedtools` | — | `bedtools` | 区间交集引擎 |
| threads / mem_mb | — | 1 / 8000(每规则) | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
python modules/peak_te_overlap/bin/intersect_te.py \
  --peak S1_peaks.narrowPeak --te-gtf TE.gtf --bedtools bedtools \
  --overlap-bed S1_peak_te_overlap.bed \
  --subfamily-tsv S1_te_subfamily_overlap.tsv --sample-id S1 \
  --ip-bam S1.sorted_markdup.bam --input-bam C1.sorted_markdup.bam

python3 modules/peak_te_overlap/bin/peak_centric_overlap.py \
  --peak S1_peaks.narrowPeak --te-gtf TE.gtf \
  --output S1_peak_centric_te.tsv --bedtools bedtools --sample-id S1

python modules/peak_te_overlap/bin/plot_te_overlap.py \
  --tsvs S1_te_subfamily_overlap.tsv,S2_te_subfamily_overlap.tsv \
  --output out/ --combined-tsv out/te_subfamily_overlap_combined.tsv \
  --method count --sort-by te_length --top-n 30 \
  --ip-input-pairs "S1:C1;S2:C2" --format png --format pdf
```

(`--output` 在 separate 模式为输出目录,combined 模式(`--combine`)为
不带扩展名的输出 base 路径;每种 `--format` 各写一份 `<base>.<fmt>`。)

## 2. 算法原理与生物学意义

- **intersect_te.py**:bedtools intersect 求 peak×TE 重叠;按 TE
  亚家族(gene_id)聚合位点级数据——重叠比例、重叠 peak 数、TE class,
  并在重叠区对 IP/Input BAM 计数(read 计数供富集分析)。
- **peak_centric_overlap.py**:peak 视角——每 peak 的 TE 数、TE class
  集合(SINE/LINE/LTR 等)、合并后 TE 覆盖碱基数与 `te_coverage_frac`。
- **plot_te_overlap.py**:上图 = 各亚家族平均 TE 长度折线,下图 =
  log2(IP/Input) 富集柱状图(按 `method`/`sort_by`/`top_n` 取 top
  亚家族);`--combine` 合并所有配对,否则每对一图。
- **生物学意义**:TE 是调控序列的重要来源(TE 派生增强子/启动子),
  峰与 TE 的重叠回答"结合/开放事件有多少落在 TE 上、哪些亚家族富集";
  IP/Input 富集排除拷贝数与可及性偏差。结果读法:`*_peak_centric_te.tsv`
  用于实验验证挑峰(含 TE 的峰);`*_te_subfamily_overlap.tsv` 与富集图
  看亚家族层面的偏好。

## 3. 参考文献

1. bedtools — https://bedtools.readthedocs.io/;Quinlan AR, Hall IM.
   *BEDTools: a flexible suite of utilities for comparing genomic
   features*. Bioinformatics 26(6):841-842, 2010.
   DOI: 10.1093/bioinformatics/btq033
2. 引擎脚本 `peak_te_overlap/bin/` 为本仓库自研(共享父目录)。
