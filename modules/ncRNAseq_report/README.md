# ncRNAseq_report — ncRNA-seq 结果汇总报告模块(PPT + 文件清单)

把 ncRNA-seq 流程的 QC 与分析产物汇成多页 PPT 报告和一份文件清单
Excel(单规则 `generate_report`,引擎为模块自带 `bin/generate_report.py`);
按 `aligner` 分支只收集该分支存在的输入。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 条件(`aligner` 分支) | 说明 |
|---|---|---|
| `common/2_trimmed_dedup_fastq/final_trimmed_fastq/{sample}/trimming_statistics_1.txt` | 总是 | 修剪统计 |
| `common/4_per_gene_bam/{sample}/{sample}.bam` + `{sample}_tail.csv` | `star_3pass_gene` | per-gene 比对 BAM 与 Tailer 尾表 |
| `results/tailer/{sample}/{sample}_tail.csv` | `star_3pass`(非 per-gene) | Tailer 输出 |
| `genome/smallrna/smallrna_genes.bed` | `star_3pass*` | 小 RNA 基因注释 BED |

以上路径相对 `outdir`(`--analysis-dir`);另收集原始 FASTQ read 计数、
demux/修剪统计与 STAR 3-pass 比对日志(脚本内部按 `--analysis-dir` 检索)。

### 1.2 输出

```
{outdir}/ncRNAseq_report.pptx           # 多页 PPT 报告
{outdir}/ncRNAseq_report_files.xlsx     # 文件清单(Excel)
{outdir}/ppt_results/                   # 报告插图目录(--img-dir)
{outdir}/ncRNAseq_report.sh             # 回执
logdir/ncRNAseq_report.log              # 运行日志
```

`report_result` 规则无 output,作 DAG 终点锚点。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `samples` | `--samples`(多值) | [] | 样本列表 |
| `paired_samples` / `single_samples` | `--paired-samples` / `--single-samples` | [] | 双端/单端名单 |
| `sample_groups` | `--sample-groups`(JSON) | {} | 样本分组 |
| `short_names` | `--short-names`(JSON) | {} | 样本简称(图中显示) |
| `aligner` | `--aligner` | `star_3pass_gene` | 比对分支(star_3pass_gene/star_3pass/…),决定输入收集 |
| `Params.report.title` | `--title` | `ncRNAseq Analysis Report` | 报告标题 |
| `Params.report.subtitle` | `--subtitle` | "" | 副标题 |
| `Params.report.pipeline` | `--pipeline` | "" | 流程描述字符串 |
| `Params.report.genome` | `--genome` | "" | 物种/基因组 |
| `Params.report.date` | `--date` | "" | 日期 |
| `Params.report.lang` | `--lang` | `zh` | 报告语言 |
| `ROOT_DIR` / `outdir` / `logdir` | `--analysis-dir` 等 | `.` / `output` / `log` | 基础键 |

线程/内存经 `rule_res("generate_report")` 配置:threads 1,mem_mb 8000。

直接调用示例(模块 bin/ 引擎脚本):

```bash
python modules/ncRNAseq_report/bin/generate_report.py \
  --analysis-dir out --output out/ncRNAseq_report.pptx \
  --file-inventory out/ncRNAseq_report_files.xlsx \
  --title "ncRNAseq Analysis Report" --lang zh --img-dir out/ppt_results \
  --aligner star_3pass_gene --samples s1 s2
```

## 2. 算法原理与生物学意义

纯报告封装:脚本收集原始 FASTQ read 计数、demux/修剪统计、STAR 3-pass
比对日志、per-gene 比对 manifest 与 Tailer 3' 尾表、小 RNA 基因注释,
绘图并排版为多页 PPTX,同时输出交付文件清单 Excel。生物学意义:
ncRNA-seq(小 RNA/非编码 RNA 测序)关注 ncRNA 的表达、3' 尾修饰与
基因特异分配;报告把"数据质量(读数/修剪/比对)→ 分析产物(tail 表、
per-gene BAM)→ 交付物清单"串成一条可审阅的证据链,便于交付与复核。

## 3. 参考文献

1. python-pptx(PPT 生成)— https://github.com/scanny/python-pptx
2. Tailer(3' 尾分析,上游产物)— https://github.com/TimNicholsonShaw/tailer
