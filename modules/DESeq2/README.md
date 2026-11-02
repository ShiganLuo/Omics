# DESeq2 模块

差异表达分析（DESeq2）+ TElocal 结果的位点级可视化。对照组和实验组至少需要两个重复。

## 规则（polygenomes/DESeq2.smk，{genome}/{contrast} 通配）

| 规则 | 输出 | 说明 |
| --- | --- | --- |
| `DESeq2_TEcount` | `diff_expression/{genome}/{contrast}/TEcount/` | TEcount（亚家族级）差异表达 |
| `DESeq2_TElocal` | `diff_expression/{genome}/{contrast}/TElocal/` | TElocal（位点级）差异表达 |
| `DESeq2_TElocal_manhattan` | `.../TElocal/manhattan/{contrast}.manhattan.png` | 全基因组 TE 位点 Manhattan 散点（爆发全局性） |
| `DESeq2_TElocal_locus_boxplot` | `.../TElocal/locus_box/{contrast}.locus_box.png` | 家族表达分布箱线/小提琴（分布上移） |
| `DESeq2_TElocal_subfamily_bar` | `.../TElocal/subfamily_bar/{contrast}.subfamily_bar.png` | 亚家族定量分组柱状（表达幅度） |

## TE 位点图（te_*.py，bin/ 下的独立 argparse CLI）

输入取流水线原始产物，不做中间预处理：DESeq2 结果 `<contrast>.TElocal_TE.tsv`、
`group.tsv`、TElocal 位点索引 `genome.references.{genome}.TEind`（.locInd，
自动缓存 parquet）；箱线/柱状另需 TElocal 计数矩阵
`counts/TEtranscripts/{genome}/{genome}_TElocal.tsv`。环境 `TEfigs.yaml`
（python + matplotlib + scipy + pyarrow + TElocal Toolkit）。

参数从 `Params.DESeq2.<rule>.<key>` 读取，缺失时回退 `Params.DESeq2.<key>`，
再回退代码默认值：

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `manhattan.label_top` | 10 | 标注 top N 位点 |
| `manhattan.ylabel` | `log2(1+\|delta TPM\|)` | y 轴标题（调用方全权控制） |
| `manhattan.highlight_families` | `{}` | 高亮家族，`{genome: {contrast: [fam, ...]}}` |
| `locus_boxplot.families` | 必填 | 作图家族，`{genome: {contrast: [pattern, ...]}}` |
| `locus_boxplot.match_field` | `subfamily` | 匹配层级 name/subfamily/family/class/any（TElocal 4 段 id 第 2 段=亚家族） |
| `locus_boxplot.match_mode` | `exact` | substring/regex/exact |
| `locus_boxplot.normalize` | `log2tpm` | log2tpm/tpm/log2cpm/cpm/none |
| `locus_boxplot.kind` | `box` | box/violin |
| `subfamily_bar.families` | 必填 | 同上 |
| `subfamily_bar.match_field` / `match_mode` | `subfamily` / `exact` | 同上 |
| `subfamily_bar.value` | `total` | total/mean（成员位点求和/均值） |
| `subfamily_bar.normalize` / `ylabel` | `tpm` / `TPM` | 标准量纲为 TPM |
| `subfamily_bar.yscale` | `linear` | linear/log |

说明：

- `manhattan` 无必需参数，随 `Params.DESeq2.enabled && Params.TElocal.enabled`
  自动进入 `outfiles`；`locus_boxplot` / `subfamily_bar` 仅当对应
  `families` 里配置了该对比时才作为目标产出（node.py 按配置追加 outfiles）。
- 选择家族按数据驱动：按显著上调位点的 delta TPM 之和排序挑选亚家族，
  宽度（位点数）× 幅度（总量）共同支撑结论，至少含一个年轻自主家族
  （IAPEz/SVA/L1HS 类）。
- 规则会生成可复现的 `.sh`（输出目录内 `DESeq2_TElocal_*.sh`），
  手动复现/调试直接跑它。
