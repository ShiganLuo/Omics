# RIPseq_report — RIP-seq 结果汇总报告模块(PPT + Excel)

把 RIP-seq 流程各环节产物(MACS3 peaks、FRiP、HOMER 注释、比对/去重/
修剪统计、TE 重叠与热图)汇成一份 PPT 报告与一份全量 Excel(单规则
`generate_report`,引擎为模块自带 `bin/generate_report.py`)。

## 1. 模块输入输出与参数

### 1.1 输入

按 `samples`(IP 样本)与 `input_samples`(Input 样本)expand 收集上游
产物;目录可经 config 覆盖以适配不同 subworkflow 布局:

| 输入 | 目录键(默认) | 说明 |
|---|---|---|
| `{sample}_peaks.narrowPeak` / `_broad_peaks.broadPeak` / `_peaks.xls` / `_summits.bed` / `_cutoff_analysis.txt` | `peaks_dir`(`{outdir}/peaks`) | MACS3 peak 产物(每 IP 样本) |
| `{sample}.FRiP.txt` | `qc_dir`(`{outdir}/QC/3_frip_score`) | FRiP 质控(每 IP 样本) |
| `{sample}_peaks.annotatePeaks.txt` | `annotation_dir`(`{outdir}/annotation`) | HOMER peak 注释(每 IP 样本) |
| `hisat2_align.log` | `log_sample_dir`(`logdir`) | 比对日志(IP+Input) |
| `{sample}.splice.txt` | `metrics_dir`(`{outdir}/common/3_raw_bam`) | 比对锚点文件 |
| `{sample}.Markdup-metrics.txt` | `markdup_dir`(`{outdir}/common/4_markdup_bam`) | MarkDuplicates 指标 |
| `trimming_statistics_1.txt` / `_2.txt` | trim_dir(`{outdir}/common/2_trimmed_fastq`,固定) | TrimGalore 修剪统计 |
| `{sample}_te_subfamily_overlap.tsv` / `{sample}_peak_centric_te.tsv` | te_dir(`{outdir}/results/te_overlap`,固定) | TE 重叠结果(每 IP 样本) |
| `{sample}_vs_{input}_enrichment.png` | te_dir | IP vs Input TE 富集图(`sample_ip_input_map` 配对) |
| `{sample}_{te}_heatmap.png` | heatmap_dir(`{outdir}/results/heatmap`,固定) | TE 热图;TE 列表取 `Params.computeMatrix.regions.genes` |

### 1.2 输出

```
{outdir}/RIPseq_report.pptx    # PPT 报告(QC 指标、图表、各模块摘要)
{outdir}/RIPseq_report.xlsx    # 全量结果数据工作簿
{outdir}/figures/<name>.<fmt>  # 报告插图独立交付(alignment、peak_frip、annotation、tss_distance、dup_rate、trim_stats、peak_centric_te)
{outdir}/figures/manifest.tsv  # 插图文件清单(name<TAB>filename)
{outdir}/generate_report.sh    # 回执
logdir/RIPseq_report.log       # 运行日志
```

`report_result` 规则无 output,作 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `samples` | `--samples`(可重复) | [] | IP 样本列表 |
| `input_samples` | `--input-samples`(可重复) | [] | Input 样本列表 |
| `sample_ip_input_map` | `--ip-input-pair ip:input`(可重复) | {} | IP → 对照 Input 映射(富集图配对) |
| `Params.report.title` | `--title` | `RIP-seq Analysis Report` | 报告标题 |
| `Params.report.subtitle` | `--subtitle` | "" | 副标题 |
| `Params.report.pipeline` | `--pipeline` | "" | 流程描述字符串 |
| `Params.report.genome` | `--genome` | "" | 物种/基因组 |
| `Params.report.date` | `--date` | "" | 日期 |
| `Params.report.top_n` | `--top-n` | 5(json:20) | 图表取 top N |
| `Params.report.formats` | `--format`(可重复) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba |
| `peaks_dir` / `annotation_dir` / `qc_dir` / `markdup_dir` / `log_sample_dir` | 同名 `--*-dir` | 见上表 | 目录覆盖(供 subworkflow 导入) |
| `ROOT_DIR` / `outdir` / `logdir` | — | `.` / `output` / `log` | 基础键 |

线程/内存经 `rule_res("generate_report")` 配置:threads 1,mem_mb 8000。

直接调用示例(模块 bin/ 引擎脚本,docstring 原文形式):

```bash
python modules/RIPseq_report/bin/generate_report.py \
  --samples Pop5IP --input-samples Pop5Input --ip-input-pair Pop5IP:Pop5Input \
  --peaks-dir ... --annotation-dir ... --qc-dir ... --log-dir ... --markdup-dir ... \
  --trim-dir ... --metrics-dir ... --te-dir ... --heatmap-dir ... \
  --output report.pptx --excel-output report.xlsx --title "RIP-seq Report" --top-n 5
```

## 2. 算法原理与生物学意义

纯报告封装,无统计方法。RIP-seq(RNA 免疫共沉淀测序)用抗体把 RNA
结合蛋白(RBP)及其结合 RNA 共沉淀后测序,以 IP vs Input 富集找出 RBP
靶标。本模块把流程各环节证据汇成两类交付:PPT 按模块分页呈现 QC 指标、
修剪/比对/去重统计、peak 概览、TE 重叠与富集图;Excel 汇集全部结果表
供后续检索。读报告时重点看:FRiP(IP 特异信号比例)、peak 数与注释
分布、IP vs Input TE 富集是否显著。

## 3. 参考文献

1. python-pptx(PPT 生成)— https://github.com/scanny/python-pptx
2. 上游工具(MACS3、HOMER、HISAT2、deepTools、Trim Galore)出处见
   各自模块 README
