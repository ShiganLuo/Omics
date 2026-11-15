# Mutation_report — 突变分析汇总报告模块(PPTX + XLSX)

把 Mutation 流程各模块结果(QC、体细胞突变、突变谱/签名、burden、
MSI、片段长度)汇编成 `Mutation_report.pptx` 与 `Mutation_report.xlsx`;
可选模块缺失时降级为「未运行」注记。自研脚本
`bin/generate_report.py`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `qc_summary` / `somatic_dir` / `spectrum_dir` / `burden_dir` / `msi_dir` / `fragment_dir` (config) | 上游各模块 `outdir`(subworkflow 按衔接规范传入;缺省回落 `{outdir}/QC/qc_summary`、`{outdir}/variation/somatic_snv_indel`、`{outdir}/stability/{spectrum,burden,msi}`、`{outdir}/results/fragment_size/fragment`) |
| `qc_dir/qc_summary.tsv` (规则 input) | QC 主表(硬依赖);其余节由脚本扫目录 |
| `outfiles` (config) | node.py 注册的流程产出清单;落在被消费目录下的文件成为 DAG 锚点(排序 + 上游更新时重跑),未注册的可选模块不会引发 MissingInputException |
| `samples` (config) | 样本 ID 列表 |
| `paired_design` (config) | tumor/normal 配对记录列表(`{"experimental": t, "control": n}`,与 gatk_somatic/msi 同形),转 `--pair t:n` |

### 1.2 输出

```
{outdir}/Mutation_report.pptx        # 汇总报告
{outdir}/Mutation_report.xlsx        # 结果工作簿(每模块一 sheet + 文件清单)
{outdir}/figures/somatic_calls_per_pair.<fmt>   # 报告自有图(按 formats 每种格式一个文件)
{outdir}/figures/manifest.tsv        # 图清单(name<TAB>filename)
{outdir}/generate_report.sh          # 回执
```

日志:`{logdir}/Mutation_report/generate_report.log`。`report_result` 为
聚合规则。

### 1.3 参数

config 位于 `Params.report` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `qc_summary` 等 6 个目录键 | `--qc-dir` / `--somatic-dir` / `--spectrum-dir` / `--burden-dir` / `--msi-dir` / `--fragment-dir` | 见输入 | 上游目录;qc 必填,其余可缺席 |
| `samples` | `--sample`(可重复) | [] | 样本列表 |
| `paired_design` | `--pair TUMOR:NORMAL`(可重复) | [] | 配对设计 |
| `Params.report.title` | `--title` | `Mutation 分析报告 / Mutation Report` | 报告标题 |
| `Params.report.pipeline` | `--pipeline` | `Mutation` | 流程名副标题 |
| `Params.report.genome` | `--genome` | "" | 基因组副标题 |
| `Params.report.date` | `--date` | "" | 日期副标题 |
| `Params.report.formats` | `--format`(可重复) | `png` | 报告自有图导出格式(多种同时输出),落 `figures/<name>.<fmt>` + `figures/manifest.tsv`;合法值: png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba(matplotlib 全集);仅长选项(`-f` 已被 `--fragment-dir` 占用) |
| (`--excel-output`) | `--excel-output` | `<output>.xlsx` | .smk 显式传 |
| (threads) | — | 1 | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python3 modules/Mutation_report/bin/generate_report.py \
    --qc-dir out/QC/qc_summary \
    --somatic-dir out/variation/somatic_snv_indel \
    --spectrum-dir out/stability/spectrum \
    --burden-dir out/stability/burden \
    --msi-dir out/stability/msi \
    --fragment-dir out/results/fragment_size/fragment \
    --output out/Mutation_report.pptx \
    --excel-output out/Mutation_report.xlsx \
    --title "Mutation 分析报告 / Mutation Report" \
    --pipeline Mutation --genome GRCh38 --date 2026-11-14 \
    --sample s1 --pair tumor:normal \
    --format png --format pdf \
    --log log/Mutation_report/generate_report.log
```

## 2. 算法原理与生物学意义

报告按节渲染各模块产物,缺失节显示「未运行」:

- **QC**:渲染 `qc_summary.tsv` 核心指标,`mapped_rate` < 80% 标红、
  `percent_duplication` > 30% 标黄(阈值定义在 `generate_report.py` 的
  `ALIGN_THRESH`/`DUP_WARN`,与 `qc_summary/metrics.md` 注明的消费方
  约定一致);
- **体细胞突变**:Mutect2 结果摘要(配对来自 `--pair`);
- **突变谱/签名、burden、MSI、片段长度**:谱图/表格直接引用上游
  产物(脚本按目录重发现文件,上游刷新后重跑报告即可更新);
- **XLSX**:每模块一个 sheet,附文件清单表,便于逐项核查。

生物学意义:一页式交付物回答「这批样本的测序质量、突变负荷、致突变
过程、MSI 状态分别如何」;锚点机制保证报告永远反映最新上游结果,
可选模块不阻塞出报告。

## 3. 参考文献

1. generate_report.py 为本仓库自研脚本(无外部文献);QC 指标定义见
   `modules/qc_summary/metrics.md`,其余节引用各上游模块 README。
