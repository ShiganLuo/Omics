# RNAseq_report — RNA-seq 分析汇总报告模块(PPTX + 文件清单)

把 RNA-seq 流程各模块结果(TE 嵌合转录本、TE/基因定量、GSVA、融合
基因、差异表达、功能富集)汇编成 `RNAseq_report.pptx`,并输出文件清单
`RNAseq_report_files.xlsx`。自研脚本 `bin/generate_report.py`。
同目录 `polygenomes/` 为多基因组变体(不在本 README 范围)。

## 1. 模块输入输出与参数

### 1.1 输入

规则 input 为上游产物全量锚点(逐一声明,缺失即无法起跑),按节分组:

| 组 | 文件(均在 `outdir` 下) | 来源模块 |
|---|---|---|
| TE 嵌合转录本 | `transcripts/TE_chimeric/TE_chimeric_{sample_summary,group_summary,te_type_counts}.tsv` + 3 张 png | StringTie/bin TEChimeric 脚本 |
| 定量矩阵 | `counts/TEtranscripts/{genome}_TEcount.tsv`、`{genome}_TEcount_name.tsv`、`counts/featureCounts/{genome}_featureCounts.tsv` | TEtranscripts / featureCounts |
| GSVA | `function/gsva/gsva_heatmap.png`、`gsva_scores.tsv` | function |
| 融合基因 | `fusion/arriba_report/`(per_sample_summary/recurrent/high_medium/inframe + figures/fig1-6) | arriba |
| 差异表达(每 contrast) | `diff_expression/{contrast}/`(group.tsv、PCA、火山、热图、upDown、TEcount 表) | DESeq2 |
| 功能富集(每 contrast) | `function/{contrast}/`(go/kegg back-to-back png + up/down csv、up/down 基因表、GSEA jpeg+csv) | function |

| 其他输入 | 说明 |
|---|---|
| `samples` / `paired_samples` / `single_samples` / `contrasts` (config) | 样本与对比列表(进报告文本与 expand) |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 上游产物根(即 `--analysis-dir`)/日志/项目根 |

注:定量矩阵三条 input 的文件名前缀取 config 的 `Params.report.genome`
标签(缺省 `all`),规则不带 `{genome}` 通配符。

### 1.2 输出

```
{outdir}/RNAseq_report.pptx          # 汇总报告
{outdir}/RNAseq_report_files.xlsx    # 文件清单
{outdir}/ppt_results/                # 报告用图片缓存(--img-dir)
{outdir}/RNAseq_report.sh            # 回执
```

日志:`{logdir}/RNAseq_report.log`。`report_result` 为聚合规则。

### 1.3 参数

config 位于 `Params.report` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.report.title` | `--title` | `RNA-seq Analysis Report` | 标题 |
| `Params.report.subtitle` | `--subtitle` | "" | 副标题 |
| `Params.report.pipeline` | `--pipeline` | `FASTQ -> TrimGalore/Cutadapt -> STAR/HISAT2 -> TEcount/StringTie/Arriba -> DESeq2 -> RNAseq_report` | 流程名副标题 |
| `Params.report.genome` | `--genome` | "" | 基因组副标题 |
| `Params.report.date` | `--date` | "" | 日期副标题 |
| `Params.report.lang` | `--lang` | `zh` | 报告语言 |
| `samples` 等 4 个列表 | `--samples` / `--paired-samples` / `--single-samples` / `--contrasts` | [] | 非空才传 |
| (outdir) | `--analysis-dir` | config `outdir` | 上游产物根 |
| (threads) | — | 1 | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python3 modules/RNAseq_report/bin/generate_report.py \
    --analysis-dir out \
    --output out/RNAseq_report.pptx \
    --file-inventory out/RNAseq_report_files.xlsx \
    --title "RNA-seq Analysis Report" --lang zh \
    --img-dir out/ppt_results \
    --samples s1 s2 --contrasts treat_vs_ctrl
```

## 2. 算法原理与生物学意义

纯汇编型报告:按节把上游产物(表 + 图)嵌入 PPT 并生成文件清单,
不重新计算任何统计。各节读法:

- **TE 嵌合转录本**:TE-基因嵌合转录本(如 TE 插入驱动的新转录本)
  的样本/组级计数与类型分布;
- **定量矩阵**:基因与 TE 表达(TEtranscripts 命名/基因组 ID 两种
  矩阵)+ featureCounts 基因计数;
- **GSVA**:通路活性评分热图(样本级通路表型);
- **融合基因**:arriba 高/中置信、复现、in-frame 融合与统计图;
- **差异表达**:每 contrast 的 PCA(样本关系/离群)、火山图(显著
  上下调)、热图与 up/down 基因/TE 表;
- **功能富集**:GO/KEGG 上下调双向条形图、GSEA 富集曲线——把差异
  基因映射到生物学过程/通路。

`lang` 控制中英;`ppt_results/` 为出图缓存。生物学意义:一套交付物
覆盖「表达了多少(定量)、差异了什么(差异表达)、意味着什么(富集/
融合/TE 嵌合)」。

## 3. 参考文献

1. generate_report.py 为本仓库自研脚本(无外部文献);各节数据出处见
   上游模块(arriba、featureCounts、TEtranscripts、DESeq2、function)
   的 README/文档。
