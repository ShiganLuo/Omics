# PeakCalling_report — ChIP-seq 峰调用总结报告模块

把 ChIP-seq 流程各环节产物(MACS3 峰、FRiP QC、HOMER 注释、比对/去重/质控
统计、TE 重叠、TE 热图)汇总为一份 PPT 报告 + 一个 Excel 全量数据簿。
单规则模块(`generate_report`),引擎为模块私有脚本 `bin/generate_report.py`,
环境为 SIF(`env.PeakCalling_report` → report.sif)或本地
`PeakCalling_report.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

全部按 `expand` 从目录键拼接(形态 A,目录键可由 subworkflow 覆写以适配
不同布局):

| 输入 | 路径模式 |
|---|---|
| MACS3 narrow/broad 峰、xls、summits、cutoff | `{peaks_dir}/{sample}/{sample}_peaks.narrowPeak`、`_broad_peaks.broadPeak`、`_peaks.xls`、`_summits.bed`、`_cutoff_analysis.txt`(IP 样本) |
| FRiP QC | `{qc_dir}/{sample}/{sample}.FRiP.txt`(IP 样本) |
| HOMER 注释 | `{annotation_dir}/{sample}/{sample}_peaks.annotatePeaks.txt`(IP 样本) |
| Bowtie2 日志/度量 | `{log_sample_dir}/{sample}/bowtie2_align.log`、`{metrics_dir}/{sample}/{sample}_bowtie2_metrics.txt`(IP+Input 全部) |
| MarkDuplicates 度量 | `{markdup_dir}/{sample}/{sample}.Markdup-metrics.txt`(全部) |
| TrimGalore 统计 | `{trim_dir}/{sample}/trimming_statistics_{1,2}.txt`(全部) |
| TE 重叠 | `{te_dir}/{sample}/{sample}_te_subfamily_overlap.tsv`、`_peak_centric_te.tsv`(IP);`{te_dir}/{sample}_vs_{input}_enrichment.png`(按 `sample_ip_input_map` 配对) |
| TE 热图 | `{heatmap_dir}/{sample}/{sample}_{te}_heatmap.png`(`te` 列表取 `Params.computeMatrix.regions.genes`) |

默认目录:`peaks_dir={outdir}/peaks`、`annotation_dir={outdir}/annotation`、
`qc_dir={outdir}/QC/3_frip_score`、`markdup_dir={outdir}/common/4_markdup_bam`、
`trim_dir={outdir}/common/2_trimmed_fastq`、`metrics_dir={outdir}/common/3_raw_bam`、
`te_dir={outdir}/results/te_overlap`、`heatmap_dir={outdir}/results/heatmap`。

### 1.2 输出

```
{outdir}/PeakCalling_report.pptx    # PPT 报告(QC 指标、图、分模块小结)
{outdir}/PeakCalling_report.xlsx    # Excel 全量结果数据簿
{outdir}/figures/<name>.<fmt>       # 报告插图独立交付(alignment、peak_frip、annotation、tss_distance、dup_rate、trim_stats、peak_centric_te)
{outdir}/figures/manifest.tsv       # 插图文件清单(name<TAB>filename)
```

日志:`{logdir}/PeakCalling_report.log`。

### 1.3 参数

| 参数 | CLI(generate_report.py) | 默认 | 说明 |
|---|---|---|---|
| `samples` | `--samples`(可重复) | `[]` | IP 样本名 |
| `input_samples` | `--input-samples`(可重复) | `[]` | Input 样本名 |
| `sample_ip_input_map` | `--ip-input-pair IP:Input`(可重复) | `{}` | IP→Input 配对(富集图) |
| `peaks_dir` / `annotation_dir` / `qc_dir` | `--peaks-dir` / `--annotation-dir` / `--qc-dir` | 见上 | 上游产物目录 |
| `log_sample_dir` | `--log-dir` | `logdir` | 比对日志目录 |
| `markdup_dir` / `trim_dir` / `metrics_dir` | `--markdup-dir` / `--trim-dir` / `--metrics-dir` | 见上 | 上游产物目录 |
| `te_dir` / `heatmap_dir` | `--te-dir` / `--heatmap-dir` | 见上 | TE 重叠/热图目录 |
| `Params.report.title` | `--title` | `Peak Calling Report`(json;smk 回退 `ChIP-seq Peak Calling Report`) | 报告标题 |
| `Params.report.subtitle` | `--subtitle` | `""` | 副标题 |
| `Params.report.pipeline` | `--pipeline` | `MACS3 -> deeptools -> report` | 流程说明 |
| `Params.report.genome` | `--genome` | `""` | 基因组名 |
| `Params.report.date` | `--date` | `""` | 日期 |
| `Params.report.top_n` | `--top-n` | `5`(json 示例 `20`) | 峰/图取前 N |
| `Params.report.formats` | `--format`(可重复) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba |
| `ROOT_DIR` | — | — | 拼 `modules/PeakCalling_report/bin/generate_report.py` |

脚本另支持 `--img-dir`(插图独立交付目录,smk 传 `{outdir}/figures`)、
`--excel-output`、`--lang {zh,en}`(默认 `zh`;smk 不传 `--lang`,用默认中文)。
PPT 嵌图始终用 png 副本(`formats` 不含 png 时用临时栅格副本)。环境:`PeakCalling_report.yaml` /
SIF `report.sif`。

直接调用示例(脱离流程):

```bash
python3 modules/PeakCalling_report/bin/generate_report.py \
  --samples s1IP --input-samples s1Input --ip-input-pair s1IP:s1Input \
  --peaks-dir peaks --annotation-dir annotation --qc-dir qc \
  --log-dir logs --markdup-dir markdup --trim-dir trim --metrics-dir metrics \
  --te-dir te --heatmap-dir heatmaps \
  --output PeakCalling_report.pptx --excel-output PeakCalling_report.xlsx \
  --title "Peak Calling Report" --top-n 20
```

## 2. 算法原理与生物学意义

纯工具封装(报告生成):脚本从各模块产物中解析质控与结果指标——TrimGalore
质控、Bowtie2 比对率、MarkDuplicates 重复率、FRiP(峰内读段占比,ChIP 富集
核心 QC)、MACS3 narrow/broad 峰数与 cutoff 曲线、HOMER 峰注释(启动子/
基因间分布)、峰-TE 重叠(TE 亚家族富集)、TE 信号热图——生成 PPT 叙述性
报告与 Excel 全量数据。结果怎么读:PPT 面向汇报(每模块一页小结,`top_n`
挑最强峰/图);Excel 面向复核(逐样本逐表全量);FRiP 低(<1% 量级)或
重复率高提示文库质量问题,TE 富集表给出峰与 TE 亚家族的重叠偏好。

## 3. 参考文献

1. MACS3 — https://github.com/macs3-project/MACS3(Zhang Y et al.,
   Genome Biology 9:R137, 2008. DOI: 10.1186/gb-2008-9-9-r137)
2. deepTools — https://deeptools.readthedocs.io/
3. HOMER(annotatePeaks)—— http://homer.ucsd.edu/homer/
