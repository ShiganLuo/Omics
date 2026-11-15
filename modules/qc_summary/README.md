# qc_summary — 全流程 QC 指标汇总模块(自研 qc_summary.py)

把 Mutation 等流程各 QC 阶段(FastQC、cutadapt、flagstat、bwa-mem2、
markdup、hs_metrics、MSI、sexcheck)的产物聚合成一张样本 × 指标表
`qc_summary.tsv`,配 HTML 可视化与逐指标柱状图。指标的完整定义、
计算原理与阈值见同目录 `metrics.md`(指标登记规范亦在其中)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `samples` (config) | 样本 ID 列表 |
| `qc_dirs` (config) | stage 名 → QC 产物目录字典;subworkflow 按衔接规范传上游各模块 config 的 `outdir`。已支持 stage:`fastqc_raw` / `fastqc_trimmed` / `cutadapt` / `flagstat` / `bwa_mem2` / `markdup` / `hs_metrics` / `msi` / `sexcheck` |
| `qc_anchors` (config) | stage → 锚点文件名模板(含 `{sample}` 为逐样本锚点,不含则为 run 级单文件如 `sexcheck.tsv`);锚点作为规则 input 声明,保证 DAG 排序并在上游产物更新时触发重跑 |
| `qc_anchor_samples` (config) | 可选 stage → 样本列表覆盖(默认全部样本;如 MSI 锚点只在 tumor 样本存在) |
| `Params.qc_summary.title` (config) | 报告标题 |

未列入 `qc_anchors` 的 stage(如流程不注册的 fastqc 输出)属尽力而为:
脚本扫目录,缺失时告警跳过;缺失 stage 的单元格记 `NA`(预期行为)。

### 1.2 输出

```
{outdir}/qc_summary.tsv        # 样本 × 指标主表(列序由 build_table() core_order 决定)
{outdir}/qc_summary.html       # 可视化汇总(引用图的 png 变体)
{outdir}/qc_summary.sh         # 回执
{outdir}/qc_summary_<metric>.<fmt>   # 逐指标柱状图,按 formats 每种格式一个文件
                             # (--plot 可选 mapped_rate/duplication/total_sequences/target_coverage;无数据的指标不出图)
```

日志:`{logdir}/qc_summary/qc_summary.log`。`qc_summary_result` 为聚合
规则。

### 1.3 参数

config 位于 `Params.qc_summary` 与顶层 `qc_dirs` / `qc_anchors` /
`qc_anchor_samples` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `samples` | `--sample`(可重复) | [] | 样本列表 |
| `qc_dirs` | `--qc-dir STAGE=DIR`(可重复) | {} | stage → QC 产物目录 |
| `qc_anchors` | —(规则 input 锚点) | {} | 见输入 |
| `qc_anchor_samples` | — | {} | 见输入 |
| `Params.qc_summary.title` | `--title` | `QC Summary` | HTML 标题 |
| `Params.qc_summary.formats` | `-f` / `--format`(可重复) | `png` | 柱状图输出格式(多种同时输出);合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集);图文件 `qc_summary_<metric>.<fmt>`,HTML 引用 png 变体;`[]` = 不出图 |
| (`--plot`) | `--plot`(可重复) | `mapped_rate duplication total_sequences` | 出图指标(.smk 未暴露,脚本默认;choices 还含 `target_coverage`) |
| (`--meta`) | `--meta` | null | 可选 run 元数据 JSON 表(.smk 未暴露) |
| (threads) | — | 1 | `rule_res("qc_summary")` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python3 modules/qc_summary/bin/qc_summary.py \
    --outdir output/QC/qc_summary \
    --title "QC Summary" \
    --sample sample1 --sample sample2 \
    --qc-dir flagstat=/path/to/flagstat_outdir \
    --qc-dir markdup=/path/to/markdup_outdir \
    --format png --format pdf \
    --log log/qc_summary/qc_summary.log
```

## 2. 算法原理与生物学意义

每个 stage 一个解析函数(`parse_fastqc` / `parse_cutadapt` /
`parse_flagstat` / `parse_bwa_mem2` / `parse_markdup` / `parse_hs_metrics`
/ `parse_msi` / `parse_sexcheck`),从该阶段原始产物的固定位置提取指标
(哪一行、怎么算、怎么看,逐指标登记在 `metrics.md`),`build_table()`
按 `core_order` 拼成样本 × 指标表,缺失填 `NA`。核心指标与判读:

- **`mapped_rate`**(flagstat,`mapped/total×100`):DNA 重测序/WES 正常
  > 95%,< 80% 报告标红;偏低查污染、参考错配、接头残留;
- **`percent_duplication`**(GATK MarkDuplicates):WGS 正常 < 15%、
  WES < 30%,> 30% 报告标黄;偏高提示文库复杂度不足;
- **`mean_target_coverage` / `pct_target_bases_30x`**(CollectHsMetrics,
  仅 exome/hybrid):WES 典型 50–150×,组内差异 > 2 倍直接影响变异检出
  灵敏度可比性;
- **`msi_score`**(MSIsensor-pro,仅 tumor):≥ 2 为 MSI-H 常用阈值;
- **`sex_sex`**(sexcheck):与 meta 性别列不符提示样本混淆;
- cutadapt 保留率、FastQC 各模块 pass/warn/fail、插入片段分布等作辅助
  诊断。

生物学意义:一份表回答"这批数据能不能进下游变异分析"——覆盖度与
  重复率决定变异检出灵敏度与负荷可比性,比对率/性别判型/MSI 兜住
  样本混淆与临床分型。消费方:`Mutation_report` 对 `mapped_rate`、
  `percent_duplication` 做红/黄着色。

## 3. 参考文献

1. qc_summary.py 为本仓库自研脚本;指标收录与阈值出处见
   `modules/qc_summary/metrics.md`。
2. 各解析对象工具:FastQC(https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)、
   cutadapt(https://github.com/marcelm/cutadapt)、samtools
   (https://github.com/samtools/samtools)、GATK
   (https://github.com/broadinstitute/gatk)、MSIsensor-pro
   (https://github.com/bioinform/msisensor-pro)。
