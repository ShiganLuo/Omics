# qc_summary 指标收录文档

本模块聚合 Mutation 流程各 QC 阶段的产物，产出 `qc_summary.tsv`（样本 × 指标表）、
`qc_summary.html`（可视化）与逐指标柱状图。本文档**稳定收录表中出现的每一个指标**：
含义、计算原理（原始文件哪一行、怎么算）、反映了什么、典型范围与告警阈值、异常时
该查什么。新增指标时必须同步在本文档登记（含义/计算/解读三要素齐全才算收录完成）。

指标来源 stage（subworkflow 经 `--qc-dir STAGE=DIR` 传入，STAGE 名与下文一致）：

| stage | 原始文件 | 解析函数 |
|---|---|---|
| `fastqc_raw` / `fastqc_trimmed` | `<sample>/*_fastqc/fastqc_data.txt`（FastQC） | `parse_fastqc` |
| `cutadapt` | `<sample>/<sample>.cutadapt_report.txt`（单端为 `.single.`） | `parse_cutadapt` |
| `flagstat` | `<sample>/<sample>_flagstat.txt`（samtools flagstat） | `parse_flagstat` |
| `bwa_mem2` | `<sample>/bwa_mem2-alignment.txt`（比对规则日志） | `parse_bwa_mem2` |
| `markdup` | `<sample>/<sample>.Markdup-metrics.txt`（GATK MarkDuplicates） | `parse_markdup` |
| `hs_metrics` | `<sample>/<sample>.hs_metrics.txt`（GATK CollectHsMetrics） | `parse_hs_metrics` |
| `msi` | `<sample>/<sample>_msi_pro.tsv`（MSIsensor-pro） | `parse_msi` |
| `sexcheck` | `sexcheck.tsv`（run 级整表，每样本一行） | `parse_sexcheck` |

表内列顺序由 `build_table()` 的 `core_order` 决定；缺失 stage/文件对应单元格保持
`NA`（可选 stage 由 node.py 的 outfiles 注册门控，属预期行为而非错误）。

消费方：`Mutation_report` 报告对 `mapped_rate`（< 80% WARN）与
`percent_duplication`（> 30% WARN）做红/黄着色（`generate_report.py` 的
`ALIGN_THRESH`/`DUP_WARN`）。

---

## 1. 测序质量（fastqc_raw / fastqc_trimmed）

raw 为原始下机数据，trimmed 为 cutadapt 去接头后；两 stage 指标同名，对比可评估
 trimming 收益。FastQC 模块状态行 `>>模块名\t状态` 收录为 `fastqc_status_<模块>`
 （`pass`/`warn`/`fail`），下表只列数值型核心项。

| 指标 | 含义 | 计算原理 | 反映了什么 / 怎么看 |
|---|---|---|---|
| `total_sequences` | 该端 reads 总数（条） | `>>Basic Statistics` 模块 `Total Sequences` 行原值 | 测序通量。PE 两端分别统计，`_1` 端文件优先收录；与 `cutadapt_read_pairs` 对照可发现文库量异常 |
| `sequence_length` | reads 长度（bp，范围如 `35-76`） | `Sequence length` 行原值 | 测序仪读长设定。trimmed 后显著变短 → 过度修剪 |
| `percent_gc` | 全reads 平均 GC 含量（%） | `%GC` 行原值 | 物种基因组 GC 的期望值（人 ~41%，小鼠 ~42%）。明显偏离 → 污染或建库偏好；fastqc `Per sequence GC content` warn 时以状态列为准 |
| `fastqc_status_<module>` | FastQC 各质检模块的判定 | `>>Per base sequence quality\tpass` 等行，状态列原样收录 | 每个模块一个布尔告警：`per_base_sequence_quality`、`per_sequence_quality`、`per_sequence_gc_content`、`sequence_duplication_levels`、`overrepresented_sequences`、`adapter_content` 等。任一非 pass 即需人工复核对应 fastqc 报告 |

## 2. 去接头（cutadapt）

来自 cutadapt 摘要（Summary）块；PE 取自 `<sample>.cutadapt_report.txt`，SE 取自
`<sample>.single.cutadapt_report.txt`。百分比取自行内 `(...%)`。

| 指标 | 含义 | 计算原理 | 反映了什么 / 怎么看 |
|---|---|---|---|
| `cutadapt_read_pairs` | 输入 read pairs 总数（PE） | `Total read pairs processed:` 行数值（去逗号） | trimming 入口通量，应 ≈ fastqc_raw `total_sequences`/2（PE） |
| `cutadapt_reads` | 输入 reads 总数（SE） | `Total reads processed:` 行数值 | 同上，SE 版 |
| `cutadapt_pairs_written_pct` | PE 保留的 pairs 百分比（%） | `Pairs written (passing filters):` 行 `(...%)` | 修剪后仍满足最小长度过滤的 pairs 比例。< 90% 常见原因：太短（插入片段 < 读长时 adapter read-through 严重，正常）、`--minimum-length` 过严 |
| `cutadapt_reads_written_pct` | SE 保留 reads 百分比（%） | `Reads written (passing filters):` 行 `(...%)` | 同上，SE 版 |
| `cutadapt_quality_trimmed_pct` | 低质量碱基调低质量修剪的 reads 占比（%） | `Quality-trimmed:` 行 `(...%)` | 测序质量健康度。大幅偏高 → 测序仪该轮质量下滑 |
| `cutadapt_bases_written_pct` | 输出碱基占输入碱基百分比（%） | `Total written (filtered):` 行 `(...%)` | 总信息保留率，综合了长度过滤与修剪 |

## 3. 比对（flagstat，作用于排序 BAM）

`parse_flagstat` 解析每行 `n + n label (pct%)`，取每个 label 的**第一个** n（QC-passed
reads 侧；`primary`/`secondary` 前缀行自动跳过——同一 label 只留首个）。计数均含
supplementary/secondary（除非指明 primary，本表未收录 primary 细分）。

| 指标 | 含义 | 计算原理 | 反映了什么 / 怎么看 |
|---|---|---|---|
| `total` | flagstat 统计总 reads 数 | `in total (QC-passed reads + QC-failed reads)` 行 | 进入比对的 reads 分母。应 ≈ trimmed 后 reads 数；偏大 → 上游多写或重复比对 |
| `mapped` | 比对上参考基因组的 reads 数 | `mapped (pct%)` 行 | 比对上的绝对量 |
| `mapped_rate` | 比对率（%） | `mapped / total × 100`，保留两位 | **核心指标**。DNA 重测序/WES 正常 > 95%；< 80% 报告标红（`ALIGN_THRESH`）。偏低查：污染（查 fastqc GC/overrepresented）、物种/参考版本错配、adapter 残留、reads 质量 |
| `properly_paired` | 双端正确配对比对 reads 数 | `properly paired (pct%)` 行 | 比对质量与文库插入分布健康度。PE 数据正常 > 90%（占 paired）；占比低 → 文库短插入/嵌合、参考不合适 |
| `paired` | 属于 pair 的 reads 数 | `paired in sequencing` 行（无百分号行，经别名归一） | PE 分母校验：应 = read1 + read2 |
| `read1` / `read2` | R1/R2 各自 reads 数 | 对应行原值（无百分号行） | 两侧应严格相等，不等 → 上游数据丢失/重复 |
| `with_itself` | mate 也比对上且同一对同时在的 reads 数 | `with itself and mate mapped` 行（别名归一） | 与 properly_paired 差值 = 配到但方向/距离异常的对 |
| `singletons` | 自身比对上但 mate 未比上 | `singletons (pct%)` 行 | 占比应 < 1%；偏高 → mate 质量差或参考不完整（剪切本 vs 全基因组） |
| `secondary` | 多重比对 secondary hits 数 | `secondary` 行（无百分号行） | 重复序列含量指示，bwa MEM 默认输出少量 |
| `supplementary` | 分段比对 supplementary 数 | `supplementary` 行（无百分号行） | 嵌合/长读取 split hits；WGS 正常几个百分点，爆高 → 污染或结构变异-heavy 样本 |
| `duplicates` | flag 标记的重复 reads 数 | `duplicates` 行 | 比对前 flag 非零即有（本流程 markdup 在下游，此处通常 0） |
| `mate_mapped_diff_chr` | mate 比到不同染色体的 reads 数 | `with mate mapped to a different chr` 行（别名归一） | 高于 ~5% 提示污染或比对错误率过高 |
| `mate_mapped_diff_chr_mapq5` | 同上但 mate MAPQ≥5 | `with mate mapped to a different chr (mapQ>=5)` 行 | 排除低质量错配后仍跨染色体 → 真 event（易位/嵌合）或污染 |

## 4. 比对日志推断（bwa_mem2，来自比对规则日志）

bwa-mem2 运行日志中 `[PE]` 统计块的推断值（日志只有首个 ~100 万 pairs 的推断，
非全数据精确值）。仅 PE 样本有值，SE 为 NA。

| 指标 | 含义 | 反映了什么 / 怎么看 |
|---|---|---|
| `insert_size_mean` / `insert_size_std` / `insert_size_median` | 文库插入片段长度分布的均值/标准差/中位数（bp） | 文库构建质量。DNA 文库中位数典型 150–400 bp；std 大 → 文库分布弥散（酶切/机械片段化不均）；中位数偏小 → 建库降解或短插入过量 |
| `bwa_unique_pairs_fr` | 唯一比对 pair 候选数（FR 方向，第二个分量） | 高质量比对 pair 的量；与 mapped 对照可估多重比对比例 |

计算原理：日志行 `[PE] # candidate unique pairs for (FF, FR, RF, RR): (n_FF, n_FR, n_RF, n_RR)` 取 FR 分量（首个 chunk 的推断值）；插入片段分布取首个 `analyzing insert size distribution for orientation FR` 块内的 `(25, 50, 75) percentile: (...)` 值元组第 2 个（中位数）与 `[PE] mean and std.dev: (mean, std)` 行。

## 5. 去重（markdup，GATK MarkDuplicates metrics）

metrics 文件 `LIBRARY` 行（表头在 `## METRICS CLASS` 之后）。

| 指标 | 含义 | 计算原理 | 反映了什么 / 怎么看 |
|---|---|---|---|
| `read_pairs_examined` | 参与去重判定的 read pairs 数（primary, mapped, 双端均比上） | LIBRARY 行 `READ_PAIRS_EXAMINED` | 去重分母；≈ mapped primary pairs |
| `read_pair_duplicates` | 被判重复的 pairs 数 | `READ_PAIR_DUPLICATES` | 与下一项同为重复来源分析 |
| `read_pair_optical_duplicates` | 光学/测序簇重复 pairs 数 | `READ_PAIR_OPTICAL_DUPLICATES` | 占 read_pair_duplicates 比例高 → 测序载玻片簇密度过高；PCR 重复则 library 复杂度低 |
| `percent_duplication` | 综合重复率（0–1 小数） | GATK 自身计算：duplicates / examined（含 single-end 折算） | **核心指标**。WGS 正常 < 15%，WES < 30%；> 30% 报告标黄（`DUP_WARN`）。偏高 → 文库复杂度不足（PCR cycles 过多、起始量不足）或数据加 Downsampling 前 |
| `estimated_library_size` | 估算的文库独特分子数 | GATK 由 duplicates–depth 关系外推（`ESTIMATED_LIBRARY_SIZE`） | 文库复杂度上限指示，越大越好；与 read_pairs_examined 差距越大 → 饱和度越低 |


## 6. 捕获质控（hs_metrics，GATK CollectHsMetrics）

仅 exome/hybrid 运行有值（node.py 在解析出 targets BED 时注册 collectHsMetrics
输出并加入本 stage 的 qc_anchors）；WGS 无捕获设计，整列为 NA，属预期。运行对象
是 analysis-ready 的 `{sample}.sorted_markdup.bam`（与变异调用同一文件），BAM 中
已标记的 duplicate 默认被排除，指标描述的是**有效可调用深度**。注意 GATK
4.6.2.0 的 CollectHsMetrics 无多线程参数（help 无 thread/parallel 类选项），
单线程运行（约 10–20 分钟/样本，40M pairs）。

| 指标 | 含义 | 计算原理 | 反映了什么 / 怎么看 |
|---|---|---|---|
| `mean_target_coverage` | 靶区平均覆盖深度（X） | `MEAN_TARGET_COVERAGE` 列原值 | **核心指标**。WES 典型 50–150×；组内样本间差异 > 2 倍会直接影响变异检出灵敏度与负荷可比性，需在下机前均衡 |
| `median_target_coverage` | 靶区覆盖深度中位数（X） | `MEDIAN_TARGET_COVERAGE` 列原值 | 比均值更抗离群（off-target 高深度区域拉高均值）；与均值差距大 → 覆盖分布右偏/存在超深热点 |
| `pct_target_bases_10x` / `_20x` / `_30x` / `_50x` / `_100x` | 覆盖深度 ≥ N× 的靶区碱基比例（0–1） | `PCT_TARGET_BASES_10X` 等列原值 | 覆盖均一性与可调用性。体细胞 SNV 检出通常要求 ≥ 30X；`pct_target_bases_30x` < 0.8 时低深度靶区变异漏检风险上升。40X 亦收录（不进 core 列序，追加在表尾） |
| `fold_enrichment` | 捕获富集倍数 | `FOLD_ENRICHMENT` 列原值 = on-near-bait 碱基比例 / 靶区占基因组比例 | 捕获效率。WES 典型 20–70×；< 10 → 探针效率低或文库降解，伴随 off-bait 升高 |
| `pct_off_bait` | off-bait 碱基占已比对碱基比例（0–1） | `PCT_OFF_BAIT` 列原值 | 落在捕获区外的测序占比。WES 正常 < 15–20%；偏高 → 捕获效率低、探头交叉杂交或文库污染 |
| `pct_selected_bases` | on+near-bait 占已比对碱基比例（0–1） | `PCT_SELECTED_BASES` 列原值 | 与 `pct_off_bait` 互补（和 ≈ 1），解读同上 |
| `target_territory` | 靶区总碱基数（bp） | `TARGET_TERRITORY` 列原值 | 捕获面板大小。同 run 内各样本应一致（本流程共用同一 targets.bed）；突变率 per-Mb 分母与此对应（98.83 Mb） |
| `bait_territory` | bait 区总碱基数（bp） | `BAIT_TERRITORY` 列原值 | 本流程 bait=target（区间由覆盖推导，无独立 bait 设计），应与 `target_territory` 相等 |
| `hs_library_size` | 杂交捕获文库复杂度估计（独特分子数） | `HS_LIBRARY_SIZE` 列原值 | 越大越好；与 `estimated_library_size`（markdup）互证，偏小 → 文库复杂度不足 |
| `total_reads` / `pf_reads` / `pf_unique_reads` | 进入 CollectHsMetrics 的总/PF/去重后 reads 数 | 对应列原值 | 与 flagstat/markdup 计数对照可核口径（PF=passing filter） |
| `on_target_bases` / `on_bait_bases` / `near_bait_bases` / `off_bait_bases` | 各类区域已比对碱基数 | 对应列原值 | 深度指标的分子，备查 |
| `mean_bait_coverage` | bait 区平均深度（X） | `MEAN_BAIT_COVERAGE` 列原值 | bait=target 时 ≈ `mean_target_coverage` |
| `at_dropout` / `gc_dropout` | AT/GC 偏好性遗漏率（%） | `AT_DROPOUT` / `GC_DROPOUT` 列原值 | 覆盖均一性的极端区域度量，WES 典型 < 5%；偏高 → GC 极端区域系统性漏检 |
| `pct_exc_dupe` / `pct_exc_off_target` / `pct_exc_mapq` | 被 duplicate/off-target/MAPQ 排除的碱基比例 | `PCT_EXC_DUPE` 等列原值 | 排除项分解，定位深度损失原因时用 |

## 7. MSI（msi，MSIsensor-pro）

仅 tumor 样本有值（node.py 按 paired_design 的 experimental 样本注册 msi 输出），
normal 及未配对样本为 NA。

| 指标 | 含义 | 计算原理 | 反映了什么 / 怎么看 |
|---|---|---|---|
| `msi_score` | 微卫星不稳定分数 | MSIsensor-pro `msi` 输出 tsv 数据行第 3 列（总位点数、可评估位点数、score） | ≥ 2 为 MSI-H 常用阈值（MSIsensor 系惯例，具体以临床指南为准）；tumor-only 无配对 normal 时可信度降低，解读需谨慎 |
| `sex_sex` | 样本性别判定 | sexcheck 模块判型（DNA 模式按性染色体覆盖度比值，RNA 模式按标记基因表达，见 §8）；XY/XX/unknown | 与 meta 表性别列不符提示样本混淆/标签错误；unknown 不判好坏，先查参考是否含 chrY / 标记基因 |
| `sex_ratio_x_to_y` | X/Y 归一化覆盖度比值 | `x_reads_per_mb / y_reads_per_mb`（inf 表示 Y 覆盖为 0） | 判型依据本身；XY 预期 <1 区间、XX 预期 inf 或极大；落中间带输出 unknown |
| `sex_xist_rpkm` / `sex_y_rpkm_sum` | RNA 模式判型指标 | XIST 与 Y 连锁标记基因的 RPKM（reads / 基因长度kb / 总 mapped 百万） | XY=y 高 xist 低、XX=xist 高 y 低；两者皆高/皆低输出 unknown |

---

## 8. 性别判定（sexcheck，run 级整表）

sexcheck.tsv 为全流程一张表（每样本一行），解析时整表读入、按 sample_id 取行，全部列以 `sex_` 前缀收录。核心列 `sex_sex` 进 core 列序；其余明细列追加在表尾备查。

sexcheck 支持两种判型逻辑，由 `Params.sexcheck.data_type` 选择：

- **DNA（默认）**：性染色体长度归一化覆盖度（reads/Mb）比值判型；XY=Y 有覆盖且比值低于阈值、XX=Y 近零或比值高、unknown=参考无 chrY 或落在模糊带。明细列 `sex_ratio_x_to_y`、`sex_x_reads_per_mb`、`sex_y_reads_per_mb`、`sex_x_mapped` 等。注意该指标对转录组数据无效（Y 基因表达稀疏，雄性样本也会落进 XX 区）。
- **RNA**（`data_type: "RNA"`，LRtranscriptome 流程默认启用）：XIST + Y 连锁标记基因（默认 Ddx3y/Eif2s3y/Kdm5d/Uty/Usp9y，按 GTF `gene_name` 大小写不敏感匹配，`xist_gene`/`y_genes` 可覆盖）表达量判型。XY=`y_rpkm_sum >= y_min_rpkm` 且 `xist_rpkm < xist_min_rpkm`、XX=镜像条件、unknown=两者皆有/皆无、标记基因缺失或样本无 BAM。明细列 `sex_xist_rpkm`、`sex_y_rpkm_sum`、`sex_y_rpkm_max`、`sex_y_rpkm_detail`、`sex_total_mapped` 等；阈值 `y_min_rpkm`（默认 0.5）、`xist_min_rpkm`（默认 0.1）可调。


## 指标登记规范（新增指标时）

1. 在对应解析函数提取该指标，并加入 `build_table()` 的 `core_order`（数值核心项
   排前，状态/计数排后）。
2. 在本文档对应 stage 章节登记一行：指标名（=tsv 列名）、含义、计算原理（原始
   文件哪一行/哪个字段、如何换算）、反映什么/怎么看（含典型范围与阈值）。
3. 需要出图的指标加进 `render_plots()` 的 `--plot` choices 并实现分支。
4. 报告着色阈值变动需同步 `Mutation_report/generate_report.py` 顶部常量并在
   本文档注明阈值出处。
