# DR_report — DR(DNA×RNA)联合分析汇总报告模块(PPTX + 文件清单)

把 DR 工作流两条链与 DR 独特分析的结果汇编成一份 `DR_report.pptx`,并输出
文件清单 `DR_report_files.xlsx`。覆盖范围:

- DNA 链(Mutation 同款分析集):QC / 体细胞 SNV/InDel / MSI / 突变谱与签名 /
  片段长度;
- RNA 链(RNAseq 同款分析集):QC / TE 与嵌合转录本 / 差异表达 / 功能富集 /
  融合基因 / GSVA;
- DR 独特分析:基因组稳定性负荷(burden)、DRclass D×R 联合证据分类、
  TransEffect 转录效应。

自研脚本 `bin/generate_report.py`。本模块是 DR 工作流唯一的报告模块
(modules.md「报告模块(_report)」规范:每个 Workflow 一个自己的 _report
模块,不复用 Mutation_report / RNAseq_report)。

## 1. 模块输入输出与参数

### 1.1 输入

硬输入(缺失即无法起跑):

| 文件 | 来源 |
|---|---|
| `{dna_dir}/QC/qc_summary/qc_summary.tsv` | qc_summary(DNA 链) |
| `{rna_dir}/QC/qc_summary/qc_summary.tsv` | qc_summary(RNA 链) |

其余全部为**锚点输入**(Mutation_report 同款锚定模式):subworkflow 把
node.py 注册的 `outfiles` 传入,凡落在消费根(`dna_dir`/`rna_dir`/
`drclass_dir`/`teffect_dir`/`burden_dir`)下的注册文件都成为该规则的
input —— 报告排在其全部生产者之后,任一节产物刷新即重跑;未启用的可选
阶段没有注册文件,不产生 MissingInputException。脚本运行期按目录扫描
消费(缺失节降级为"未运行"注记,绝不报错)。

| 消费根(默认派生) | 主要消费内容 | 来源模块 |
|---|---|---|
| `dna_dir`(`{outdir}/DNA`) | QC/sexcheck、`variation/somatic_snv_indel/{N_vs_T}/{N_vs_T}.filtered.vcf.gz`、`stability/msi/msi_input.txt`、`stability/spectrum/`(spectrum/signature 图与 csv)、`results/fragment_size/fragment/FragmentSize.png` | qc_summary、gatk_somatic、msi、somatic_spectrum、fragment_size |
| `rna_dir`(`{outdir}/RNA`) | QC/sexcheck、`transcripts/{genome}/TE_chimeric/`、`counts/`、`diff_expression/{genome}/{contrast}/`、`function/{genome}/{contrast}/`、`fusion/{genome}/arriba_report/` | qc_summary、StringTie、TEtranscripts、featureCounts、DESeq2、function、arriba |
| `burden_dir`(`{dna_dir}/stability/burden`) | `burden_summary.tsv` + 每样本 `*.snv_indel/sv/cnv_burden.tsv` | burden |
| `drclass_dir`(`{outdir}/DRclass`) | `DRclass_summary.tsv` + `{sample}/{sample}.{dna,rna}_anchored.tsv` | DRclass |
| `teffect_dir`(`{outdir}/TransEffect`) | `rna_support_sites/vaf_ratios/variant_gene_delta/cnv_dosage_summary/cnv_segments/program_overlap.tsv` + 4 张 fig | TransEffect |

| 其他输入 | 说明 |
|---|---|
| `samples` / `dna_samples` / `rna_samples` / `paired_samples` / `single_samples` (config) | 样本列表(进报告文本) |
| `paired_design` (config) | `{experimental, control}` 肿瘤:正常设计对(体细胞计数) |
| `contrasts` (config) | DE 对比键(`sample_params.comparisons[<genome>]` 的键) |
| `genome_label` (config) | RNA 树 genome 目录标签(如 `GRCm39`);缺省取 `Params.report.generate_report.genome` |
| `outfiles` (config) | node.py 注册输出(锚点来源) |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 交付物根/日志/项目根 |

### 1.2 输出

```
{outdir}/DR_report.pptx          # 联合汇总报告(15 页)
{outdir}/DR_report_files.xlsx    # 文件清单(逐产物多 sheet)
{outdir}/figures/<name>.<fmt>    # 报告自有图独立交付(somatic_calls_per_pair、
                                 # burden_overview、drclass_label_composition、
                                 # de_updown_per_contrast;--img-dir)
{outdir}/figures/manifest.tsv    # 插图文件清单(name<TAB>filename)
{outdir}/DR_report.sh            # 回执(RECEIPT_DIR 另有带时间戳历史副本)
```

日志:`{logdir}/DR_report.log`。`report_result` 为聚合规则。

### 1.3 参数

config 位于 `Params.report.generate_report` 下(规则级精确读点):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `title` | `-T/--title` | `DNA×RNA 联合分析报告 / DR Joint Analysis Report` | 标题 |
| `subtitle` | `-S/--subtitle` | "" | 副标题 |
| `pipeline` | `-L/--pipeline` | 模块内默认链描述 | 流程标注 |
| `genome` | `-g/--genome` | node.py 写回全部样本 organism | 副标题基因组标注 |
| (config `genome_label`) | `-G/--genome-label` | `genome` 兜底 | RNA 树目录标签 |
| `date` | `-y/--date` | node.py 写回运行当天 | 日期 |
| `lang` | `-l/--lang` | `zh` | 报告语言(zh/en) |
| `formats` | `-f/--format` | `png` | 自有图导出格式(多值;[] 不出图) |

脚本 CLI(脱离流程直跑)其余参数:`-a/--analysis-dir`(必填,workflow 输出根)、
`-d/-r/-k/-t/-b`(五个消费根覆盖)、`-o/--output`(必填 PPTX)、
`-e/--file-inventory`(XLSX,缺省 `<output stem>_files.xlsx`)、`-i/--img-dir`
(缺省 `<output dir>/figures`)、`-s/-D/-R/-p/-u`(样本列表,repeatable)、
`-c/--contrast`(对比,repeatable)、`-P/--pair TUMOR:NORMAL`(体细胞设计对,
repeatable)、`-m/--log`。

## 2. 报告怎么读

15 页 PPTX,标题带 ★ 的是 DR 独特分析:

| 页 | 内容 | 数据源 |
|---|---|---|
| 1-2 | 标题、流程总览(两链 + 独特分析卡片) | 元信息 |
| 3 | DNA 质控概览(qc_summary + sexcheck 表) | DNA/QC |
| 4 | DNA 变异概览(体细胞计数表+柱状图、MSI 表) | gatk_somatic、msi |
| 5 | 突变谱与突变签名(嵌入 spectrum/signature 图) | somatic_spectrum |
| 6 | 片段长度分布(FragmentSize 图) | fragment_size |
| 7 ★ | 基因组稳定性负荷(burden_summary 表 + SNV/SV/CNV 三联柱状图) | burden |
| 8 | RNA 质控概览(qc_summary + sexcheck 表) | RNA/QC |
| 9 | TE 与嵌合转录本(sample_summary 表 + 类型 top 图) | StringTie TEChimeric |
| 10 | 差异表达概览(每对比 up/down 表 + 柱状图) | function 阶段基因表 |
| 11 | 功能富集概览(GO/KEGG top 词条 + GSVA 热图) | function |
| 12 | 融合基因概览(per_sample_summary 表 + fig1) | arriba |
| 13 ★ | DRclass D×R 联合证据分类(summary 表 + 标签堆叠柱状图) | DRclass |
| 14 ★ | TransEffect 转录效应(vaf_ratios 表 + 4 图) | TransEffect |
| 15 | 结论(样本/对比/独特分析计数) | 汇总 |

XLSX 逐产物多 sheet(dna_qc、somatic_calls、msi_scores、burden_summary、
drclass_summary、drclass_<sample>_{dna,rna}_anchored、teffect_*、rna_qc、
te_chimeric_*、de_<contrast>_*、func_<contrast>_*、fusion_*、file_inventory);
多样本合并表首列带 sample_id;行数截断 5000。

## 3. 怎么算(算法与口径)

- **体细胞计数**:对每个设计对 `{normal}_vs_{tumor}` 流式统计 Mutect2
  `*.filtered.vcf.gz` 的非 header 行数(不重算、不过滤,口径 = 过滤后 VCF
  行数);VCF 缺失记 NA。
- **MSI**:读 `msi/msi_input.txt`(`sample<TAB>score`)。
- **burden 页**:直接呈现 `burden_summary.tsv`(外连接宽表:snv/indel/tmb
  计数、SV 分型计数、CNV gain/loss 段计数),柱状图为每样本
  `snv_count` / `sv_total` / `cnv_total_segments` 三联。口径见
  modules/burden/README.md(VAF 阈值、CNV altered 判定)。
- **DRclass 页**:`DRclass_summary.tsv` 每样本一行,12 个 D×R 标签计数
  (D1R1…D0R1、R1-ANCHOR-REJ)+ 穿透率 `penetrance=(R1+R2)/(R1+R2+R6)` +
  `r4_fraction=D1R4/(D1R1+D1R4)` + QC(WARN 当 r4_fraction > max_r4)。
  堆叠柱状图为每样本标签构成。标签语义与判定状态机见
  modules/DRclass/README.md。
- **TransEffect 页**:直接呈现 `vaf_ratios.tsv` 与 4 张成图
  (RNA 支持位点 VAF、变异基因 Δ'、CNV 剂量一致性、稳定性程序基因重叠)。
- **差异表达计数**:每对比 `function/{genome}/{contrast}/{up,down}_genes.txt`
  的非空行数(即功能分析消费的 DE 基因表口径,与富集分析一致)。
- **表格列宽**:excel 自适应宽度对全 NaN 列安全(pandas 3.x 下
  `astype(str)` 保留缺失,宽度计算用 `pd.isna` 守卫)。
- **降级行为**:任何节缺输入/未启用 → 该页"未运行"注记 + xlsx 空 sheet,
  脚本退出码 0;仅当两个硬输入 qc_summary.tsv 缺失时规则无法起跑
  (DAG 拦截)。

## 4. 参考文献

- modules/modules.md「报告模块(_report)」(目录结构、两规则、subworkflow 集成规范)
- modules/burden/README.md(burden 口径)
- modules/DRclass/README.md(D×R 标签空间与状态机;设计文档 `3.2.2.3节写作构思.md`)
- modules/TransEffect/README.md(转录效应指标定义)
- modules/Mutation_report、modules/RNAseq_report(报告房式 pptx 助手与数据收集器蓝本)
