# RNAseq_report/polygenomes — 多基因组 RNA-seq 总结报告子模块

把 RNA-seq 流程各模块产物(TE 嵌合、TEtranscripts/featureCounts 定量、
GSVA、arriba 融合、DESeq2 每对比结果、GO/KEGG/GSEA)汇总为 PPT 报告 +
文件清单 Excel。C 形态子模块,共享父目录 `../RNAseq_report.yaml` 环境与
`modules/RNAseq_report/bin/generate_report.py` 脚本;环境为 SIF
(`env.RNAseq_report` → report.sif)。

## 1. 模块输入输出与参数

### 1.1 输入

全部按 genome/contrast 从 `indir`(流程根输出目录)拼接,主要来源:

| 输入 | 路径模式 |
|---|---|
| TE 嵌合 | `transcripts/{genome}/TE_chimeric/TE_chimeric_{sample_summary,group_summary,te_type_counts}.tsv` + 3 张 png |
| TE 定量 | `counts/TEtranscripts/{genome}/{genome}_TEcount{,_name}.tsv` |
| 基因定量 | `counts/featureCounts/{genome}/{genome}_featureCounts.tsv` |
| GSVA | `function/{genome}/gsva/gsva_heatmap.png`、`gsva_scores.tsv` |
| 融合(arriba) | `fusion/{genome}/arriba_report/{per_sample_summary,recurrent_fusions,high_medium_confidence_fusions,inframe_fusions}.tsv` + `figures/fig1..fig6` |
| 每对比 DESeq2 | `diff_expression/{genome}/{contrast}/TEcount/` 下 group.tsv、PCA、火山图、热图、upDown 三表、Gene/TE/Gene_TE name.tsv |
| 每对比功能富集 | `function/{genome}/{contrast}/` 下 go/kegg 背靠背图与 up/down csv、up/down genes、GSEA jpeg/csv |

contrast 列表取 `Params.DESeq2.group_pairs.<genome>` 的键;脚本运行前规则
把 `indir` 下产物软链成脚本约定的 `analysis_dir/{transcripts,fusion,counts,...}`
布局。

### 1.2 输出

```
{outdir}/{genome}/RNAseq_report.pptx          # PPT 报告
{outdir}/{genome}/RNAseq_report_files.xlsx    # 文件清单
{outdir}/{genome}/figures/<name>.<fmt>        # 报告插图独立交付(chimeric_ratio、de_summary、venn_gene_up/down、venn_te_up/down)
{outdir}/{genome}/figures/manifest.tsv        # 插图文件清单(name<TAB>filename)
```

日志:`{logdir}/RNAseq_report/{genome}/RNAseq_report.log`。

### 1.3 参数

| 参数 | CLI(generate_report.py) | 默认 | 说明 |
|---|---|---|---|
| `genome_paired_samples` / `genome_single_samples` | `--samples` / `--paired-samples` / `--single-samples` | `{}` | genome → 样本列表 |
| `Params.DESeq2.group_pairs.<genome>` | `--contrasts`(键列表) | `{}` | 对比组定义 |
| `Params.report.title` | `--title` | `RNA-seq Analysis Report` | 标题 |
| `Params.report.subtitle` | `--subtitle` | `""` | 副标题 |
| `Params.report.pipeline` | `--pipeline` | `""` | 流程说明 |
| `Params.report.genome` | `--genome` | wildcard | 基因组名 |
| `Params.report.date` | `--date` | `""` | 日期 |
| `Params.report.lang` | `--lang` | `zh` | 报告语言(zh/en) |
| `Params.report.formats` | `--format`(可重复) | `["png"]` | 图片格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba |
| 固定 | `--analysis-dir`、`--file-inventory`、`--img-dir` | — | 见上 |

环境:`conda: "../RNAseq_report.yaml"` / SIF。

直接调用示例(脱离流程):

```bash
python3 modules/RNAseq_report/bin/generate_report.py \
  --analysis-dir report_staging/GRCh38 \
  --output out/GRCh38/RNAseq_report.pptx \
  --file-inventory out/GRCh38/RNAseq_report_files.xlsx \
  --title "RNA-seq Analysis Report" --genome GRCh38 --lang zh \
  --img-dir out/GRCh38/figures \
  --samples s1 s2 --paired-samples s1 --single-samples s2 \
  --contrasts aged_vs_youth
```

## 2. 算法原理与生物学意义

纯工具封装(报告生成):按报告叙事顺序组织各模块结果——样本/文库 QC、
定量矩阵、PCA 与差异火山/热图、TE 表达与 TE 嵌合、融合基因、GSVA 通路
活性、GO/KEGG/GSEA 富集——生成 PPT(图 + 小结)与文件清单 Excel(所有
引用文件的路径清单,便于溯源交接)。结果怎么读:PPT 按对比组组织,重点
图(火山、热图、富集背靠背图)配结论位;Excel 是产物索引。

## 3. 参考文献

1. python-pptx(报告生成)—— https://github.com/scanny/python-pptx
2. 上游工具引证见对应模块:DESeq2、StringTie、arriba、GSVA、TEtranscripts
