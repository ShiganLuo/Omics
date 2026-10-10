# Fiber-seq 下游分析方案（供审阅）

> 目标：在现有上游流程（m6A/核小体/MSP/FIRE calling + 比对 + 提取 + QC）之上，
> 增加下游解读层。分析项均来自三篇指定文献中实际使用的方法 [1][2][3]。
> 每一项标注：来源文献 → 分析内容 → 拟用工具 → 落点（复用模块 / 新增模块）→ 优先级。

## 实施状态（已按审阅决定实现）

审阅决定：1 保留 `ft fire`（补 FDR/haps）；2 LDSC 需要、保留接口；3 单倍型要做；
4 着丝粒要做（参考基因组可配）；5 motif 库择优（选 MEME/JASPAR 自扫描，任意物种）；
6 组间比较按 [3] 设计。

已落地（2026-10-09）：
- **上游修正**：`ft call-peaks` 默认 FDR 模式（bedtools-shuffle 等价的纯 Python null +
  `--shuffled`，输出 `fire_fdr_table.tsv`），`--haps` 可选；`ft fire --extract` 产出
  per-fiber FIRE BED；HP 拆分 + 分相提取规则。
- **新模块**（初版为 fiberseq_analysis 聚合模块，后按规范拆分，见下）：annot/motif/
  footprint/nucpos/quant/diff/coactuation/haplotype/censat/track。
- **新模块** `ldsc`（接口：annot/munge/h2-cts）；富集分析复用 `modules/function`。
- **接线**：Fiberseq.smk Step 4–7、node.py 按 `Params.*.enabled` 注册输出、
  config/schema 扩展；待补文件清单见 `docs/fiberseq-downstream-files.md`。
- **验证**：13 个 CLI 合成数据冒烟测试通过（含 Enrichr 真实 API 路径）；
  snakemake dry-run 全链 DAG 构建通过。
- **规范化重构（modules.md 合规）**：enrichment 重复实现已删除、复用 modules/function
  （function_go_kegg/function_gsea，diff 直接产出 go-kegg.r 兼容基因表）；
  下游分析按布局 C 组织为 `modules/fiberseq/` 父模块 + 10 个子模块
  （annotate/motif/footprint/nucpos/quant/diff/coactuation/haplotype/censat/track），
  共享父目录 fiberseq.yaml（一个 SIF），共享库 src/common/util/FiberlibUtil.py；
  规则体全部为 modules.md 标准模板；输出目录 `{outdir}/<module>/{sample_id}/`。
  复测：冒烟 47/47、wiring dry-run 全部入 DAG。

---

## 1. 现状

当前 `subworkflow/Fiberseq.smk` 三步：

1. `ft predict-m6a`（SPRQ 数据自动跳过）→ `ft add-nucleosomes` → `ft fire`
2. `pbmm2` 比对
3. `ft extract --all` + `ft call-peaks` + `ft qc`

输出止于 `results/<sample>/{m6a,nuc,msp,cpg}.bed.gz`、`fire_peaks.bed`、`qc.tsv`，
即"数据加工 + peak calling + QC"。**peak 之后的注释、motif/footprint、核小体
定位、单倍型、共活化、重复/着丝粒专题等解读层全部缺失。**

---

## 2. 三篇文献的下游分析方法盘点

### 2.1 人脑细胞类型特异 Fiber-seq（Cell Rep Methods 2024）[1]

该文 STAR Methods 中列出的下游分析（均已复现于 pyft/pandas 脚本）：

| 分析 | 方法细节 |
|---|---|
| FIRE 注释与 peak calling | FIRE 半监督流程给每个 MSP 打 precision 分，≥0.9 视为调控元件；基于打乱坐标的 null 分布按 FDR 5% 调 peak，过滤 <10 fibers 的 peak [1] |
| Peak 聚类与注释 | pseudo-bulk 信号 k-means (k=4) 聚类，deepTools 可视化，ChIPseeker 注释到 promoter/enhancer/intron 等 [1] |
| Repeat 注释 | bedtools overlap，peak ≥25% 长度落在 repeat 即标注重复序列 [1] |
| Motif 富集 | HOMER known motif 流程按 peak cluster 跑 motif 富集（NFY/SP1/CTCF/SOX9/NEUROD1 等）[1] |
| TF footprint 评分 | motif 中心计算 m6A/A 聚合谱，footprint score = 两侧最大值之和 − motif 区 m6A 总量（GC 校正），KS 检验比较分组 [1] |
| TF footprinting（单纤维） | `ft footprint` 在 CTCF 35 bp motif 上做 footprint；fiber 分三类（MSP+footprint / MSP 无 footprint / 核小体覆盖），±1 kb 聚合 m6A/A 谱；观察到核小体 phasing 随距 motif 距离衰减 [1][2] |
| 核小体定位 | TSS/调控元件 ±5 个核小体中心 offset 聚合（随机参考 fiber 相对化），Wilcoxon 检验 -1/+1 核小体 [1] |
| Co-actuation（共活化） | 同一 fiber 上两个位点均含 precision≥0.9 的 FIRE 元件；Fisher 检验 + actual vs expected 共活化比例（显著对 0.149 vs 0.036）[1] |
| 单倍型分析 | DeepVariant + pbsv (+sniffles) 调变异，hiphase 分相（meryl/merqury k-mer 辅助）；peak 上 H1/H2 可及性差异 + Fisher 检验 [1][2] |
| 功能富集 | 非启动子 peak 按可及性排序取 top1000 基因，Enrichr GO-BP 富集（p<0.05）[1] |
| GWAS 富集（可选） | LD score regression（去 repeat peak、去 MHC 区）[1] |
| 与短读长整合 | ATAC-seq/ChIP-seq 信号差异与 Fiber-seq percent accessibility 相关性分析 [1] |

### 2.2 fibertools 工具文（Genome Research 2024）[2]

| 分析 | 方法细节 |
|---|---|
| m6A 精度分级 | ML tag 阈值控制 precision（>95% 用 244/256），可按精度要求调阈值 [2] |
| 核小体 calling 理解 | HMM（A/T 碱基）+ 简单 >85 bp 无修饰段细化；核小体/MSP 定义影响下游统计口径 [2] |
| 遗传-表观整合 | 单分子遗传+表观共处理、molecule↔reference 坐标互转、SV 区域内分析；分相信息写入 BAM haplotag 字段后可直接做 haplotype-aware 分析 [2] |
| ONT 支持 | Dorado 6mA/5mC calling 走同一 ft 流程（ONT 数据加 `Params.fibertools.ont=true`）[2] |

### 2.3 着丝粒 DNA 甲基化（Nat Genet 2025）[3]

| 分析 | 方法细节 |
|---|---|
| 区域聚焦提取 | 按 chm13v2.0_censat_v2.1.bed（α-sat HOR / CDR）用 samtools view + bedtools intersect 提取 core/noncore reads，分组对比 [3] |
| 修饰定量 | modkit call-mods 动态阈值；modkit pileup 按位点计数 5mC/5hmC（阈值 0.8），bedGraph→bigWig 滚动 mCpG/CpG 谱 [3] |
| 可及性定量 | `ft add-nucleosomes` + `ft extract`，计算 >50 bp MSP 占比作为可及性指标（自定义脚本）[3] |
| 重复区无参考甲基化 | 着丝粒 reads 转 fasta 重新比对 / mapping-independent 甲基化统计（适配高度重复区）[3] |
| 蛋白定位整合 | DiMeLo-seq（CENP-A、H3K9me3）与 Fiber-seq 联合解读染色质状态 [3] |
| 亚硫酸盐验证 | COBRA 实验验证（湿实验，不入流程）[3] |

---

## 3. 下游分析模块方案

原则（遵循 modules.md/subworkflow.md 规范）：一个模块一个原子分析步骤；逻辑在
`bin/*.py`；`node.py` 按 config 开关注册输出文件，用输出决定流程走向；
新脚本遵循 omics-script-engineering 规范（argparse 全参数化、LogUtil、README）。

### Tier 1 — 核心解读层（P0，覆盖 90% 常规需求）

**A. FIRE peak 注释与分类**（来源 [1]）
- 内容：ChIPseeker/bedtools 注释 promoter/enhancer/intron/intergenic；
  repeat overlap（≥25% 规则，复用 `peak_te_overlap` 模块思路）；
  基于 percent accessibility 的 peak 聚类 + 汇总表。
- 工具：ChIPseeker（R）或 bedtools + GENCODE GTF；deepTools 聚类图。
- 落点：新模块 `fire_annotate`（bin/fire_annotate.py）；聚类图复用 `deeptools`。
- 输出：`results/<sample>/annotate/{peaks_annotated.tsv, peaks_repeat_flag.tsv, cluster/}`。

**B. Motif 富集 + motif 中心可及性谱**（来源 [1]）
- 内容：peak cluster 上 HOMER known-motif 富集；motif 位点中心 ±N bp 的
  m6A/A 聚合谱（footprint 初筛）。
- 工具：HOMER `findMotifs.pl`（现 `homer` 模块只有 annotatepeaks，需补
  findMotifs 规则）+ 新聚合脚本 `motif_profile.py`。
- 落点：`homer` 模块补规则；新模块 `fiberseq_footprint` 的 profile 部分。
- 输出：`results/<sample>/motif/{homer/, motif_m6a_profile.{tsv,pdf}}`。

**C. 核小体定位分析**（来源 [1]）
- 内容：以 TSS / FIRE peak 中心为锚点，取 ±5 个核小体中心 offset 聚合
  （随机参考 fiber 相对化），核小体 phasing 强度随距离衰减曲线，Wilcoxon 统计。
- 工具：新脚本 `nucleosome_positioning.py`（pandas/numpy，纯计算）。
- 落点：新模块 `fiberseq_nucpos`。
- 输出：`results/<sample>/nucpos/{nuc_offset_boxplot.pdf, nuc_phasing_decay.pdf, nuc_offset.tsv}`。

**D. 基因组轨道与可视化**（来源 [1][3]）
- 内容：BED12 → bigBed/bigWig（percent accessibility、m6A density、
  nucleosome occupancy 信号），UCSC/IGV 轨道。
- 工具：UCSC bedToBigBed/bedGraphToBigWig；复用 `track` 模块。
- 落点：`track` 模块扩展 fiberseq 轨道规则。
- 输出：`results/<sample>/track/*.big{Bed,Wig}` + UCSC trackDb。

**E. 可及性定量矩阵**（来源 [1][3]）
- 内容：每个 peak/区域的 percent accessibility（可及 fiber 比例）、
  MSP>50 bp 占比 [3]、m6A/A 比例 → 样本×区域矩阵，供分组比较。
- 工具：新脚本 `accessibility_matrix.py`；组间差异用 Wilcoxon/Fisher。
- 落点：`fiberseq_nucpos` 同模块或独立 `fire_quant`。
- 输出：`results/matrix/{accessibility_matrix.tsv, diff_peaks.tsv}`。

### Tier 2 — 高级分析（P1）

**F. TF footprinting（单纤维三分类）**（来源 [1][2]）
- 内容：`ft footprint` 在指定 motif（如 CTCF 35 bp）上调用 footprint；
  fiber 三分类（MSP+footprint / MSP 无 footprint / 核小体覆盖）±1 kb 聚合谱；
  footprint score（GC 校正）+ KS 检验。
- 工具：`ft footprint`（ft ≥0.13 是否已含该子命令待实测确认）+ 聚合脚本。
- 落点：`fibertools` 模块补 `ft_footprint` 规则 + 新模块 `fiberseq_footprint`。
- 依赖：Tier 1-B 的 motif 位点集合 + ChIP-seq 验证集（可选输入）。

**G. 单倍型特异染色质**（来源 [1][2]）
- 内容：DeepVariant + pbsv 调变异 → hiphase 分相 → BAM 写入 haplotag →
  按 haplotype 分别统计 FIRE/可及性差异（Fisher 检验）。
- 工具：复用现有 `deepvariant`、`pbsv`、`hiphase` 模块。
- 落点：`Fiberseq.smk` 串联上述模块 + 新脚本 `haplotype_fire_diff.py`。
- 注意：成本高（全基因组变异检测），默认关闭、按需开启。

**H. Co-actuation 共活化分析**（来源 [1]）
- 内容：同一 fiber 上 FIRE 元件对的共现统计（Fisher + actual vs expected
  差值，过滤覆盖 <10 fibers 的对）。
- 工具：新脚本 `coactuation.py`（读 BAM fire 标签或 fire BED12）。
- 落点：新模块 `fiberseq_coactuation`；输入可接受 promoter–enhancer 候选对
  列表（来自外部 ABC/Hi-C 数据）或全 peak 对（限量）。
- 输出：`results/coactuation/{coactuation_pairs.tsv, coactuation_stats.pdf}`。

**I. 功能富集**（来源 [1]）
- 内容：非启动子 peak 按可及性排序取 top N 基因 → GO-BP/KEGG 富集
  （Enrichr 或本地 clusterProfiler）。
- 工具：现 `enrichment` 模块是空壳（`# need test`），需实现。
- 落点：补全 `enrichment` 模块。
- 备注：GWAS LDSC [1] 计算与依赖（GWAS summary statistics、LD reference）
  重，列为可选插件，不在常规流程。

### Tier 3 — 专题分析（P2，按项目需要）

**J. 着丝粒/重复区域专题**（来源 [3]）
- 内容：censat BED（α-sat HOR/CDR）区域聚焦提取、core vs noncore 分组；
  modkit pileup 5mC/5hmC 定量与滚动谱；无参考/重比对甲基化统计；
  CDR 边界 CENP-A/H3K9me3（DiMeLo-seq）整合位点图。
- 工具：samtools/bedtools/modkit/minimap2 + 新脚本 `censat_profile.py`。
- 落点：新模块 `fiberseq_censat`；依赖 T2T-CHM13v2.0 参考与
  chm13v2.0_censat_v2.1.bed 注释（config genome references 增加
  `censat_bed` 字段）。
- 附带：ONT 数据用户走 modkit 路径 [2][3]；PacBio CpG 用 `ft extract` 的
  cpg.bed.gz。

**K. 多组学整合**（来源 [1][3]）
- 内容：与 ATAC-seq/H3K4me3/H3K27ac/CENP-A 等信号的聚类整合、差异谱；
  与 scRNA-seq 细胞类型数据联合解读 [1]。
- 工具：deepTools + 复用 `PeakCalling` 相关模块产出。
- 落点：依赖外部数据登记（data_registry），按项目配置。

---

## 4. 工程落点汇总

| 项 | 复用 | 新增 |
|---|---|---|
| A 注释 | `peak_te_overlap`（repeat）, `deeptools`（聚类图） | `fire_annotate` 模块 |
| B motif | `homer`（补 findMotifs 规则） | `motif_profile.py` |
| C 核小体定位 | — | `fiberseq_nucpos` |
| D 轨道 | `track` | fiberseq 轨道规则 |
| E 定量矩阵 | — | `fire_quant` |
| F footprint | `fibertools`（补 ft_footprint） | `fiberseq_footprint` |
| G 单倍型 | `deepvariant` `pbsv` `hiphase` | `haplotype_fire_diff.py` |
| H 共活化 | — | `fiberseq_coactuation` |
| I 富集 | `enrichment`（需实现，现为空壳） | — |
| J 着丝粒 | `samtools` `bedtools` | `fiberseq_censat`（modkit 容器） |

其他改动：
- `node.py runFiberseq`：按 `Params.*` 开关注册对应 outfiles
  （输出决定流程走向，开关逻辑只在 node.py 层）。
- `config/Fiberseq.schema.json`：增加 `Params` 下的下游分析块
  （annotate/motif/nucpos/track/quant/footprint/haplotype/coactuation/
  enrichment/censat 各 bool + 参数）；`genome.references` 增加
  `gtf`、`censat_bed`、`motif_db` 等可选字段。
- `Fiberseq.smk`：Step 4+ 按 use rule 引入新模块，严格按 subworkflow.md
  的 `*_config` 字典模式。
- 新脚本统一：argparse 全参数化、短选项、LogUtil 日志、`--dry-run`、
  skip-existing、README（omics-script-engineering 规范）。
- 容器：modkit/HOMER/ChIPseeker 各自独立 SIF（一个模块一个环境）。

---

## 5. 分期实施

| 阶期 | 内容 | 验收标准 |
|---|---|---|
| Phase 1 | A 注释 + D 轨道 + E 定量矩阵 | fire_peaks.bed 全部注释归类；bigWig 可在 IGV 打开；矩阵含全部样本 |
| Phase 2 | B motif + C 核小体定位 + I 富集 | HOMER 富集表 + motif 中心谱图；TSS ±5 核小体 offset 图；GO 富集表 |
| Phase 3 | F footprint + G 单倍型 + H 共活化 | CTCF 三分类谱复现文献形态；haplotype 差异表；共活化统计表 |
| Phase 4 | J 着丝粒专题 + K 多组学 | censat profile 图；整合聚类图 |

Phase 1–2 覆盖文献 [1] 的常规图表主体，建议先做。

---

## 6. 待确认问题（审阅时请定夺）

1. **FIRE 实现**：继续用 `ft fire`（已有），还是引入官方 FIRE 流水线 + pyft
   （文献 [1] 所用，precision 分级、haplotype-aware aggregate score 更完整）[4]？
   建议：先保持 `ft fire`，Phase 3 评估官方 FIRE。
2. **GWAS LDSC** 是否需要？依赖外部 GWAS 汇总统计与 LD panel，建议默认不做。
3. **单倍型分析**默认开还是按需开？（变异检测 + 分相成本高）
4. **着丝粒专题**是否要求 T2T-CHM13v2.0 参考？现有 GRCh38 项目无法做
   α-sat HOR 级别分析 [3]。
5. **motif 数据库**：HOMER known motifs 默认即可，还是接入 JASPAR/HOCOMOCO？
6. 比较分析（处理组 vs 对照，如文献 [3] 的 demethylation 时间序列）是否需要
   在流程内支持成对比较设计，还是下游用 R 即可？

---

## Sources

[1] https://doi.org/10.1016/j.crmeth.2024.100911 — Single chromatin fiber profiling and nucleosome position mapping in the human brain (Cell Rep Methods 2024)
[2] https://doi.org/10.1101/gr.279095.124 — DNA-m6A calling and integrated long-read epigenetic and genetic analysis with fibertools (Genome Res 2024)
[3] https://doi.org/10.1038/s41588-025-02324-w — DNA methylation influences human centromere positioning and function (Nat Genet 2025)
[4] https://github.com/fiberseq/FIRE/blob/main/docs/README.md — FIRE analysis pipeline documentation
