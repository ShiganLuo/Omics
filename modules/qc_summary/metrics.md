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
| `msi` | `<sample>/<sample>_msi_pro.tsv`（MSIsensor-pro） | `parse_msi` |

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

## 5b. samtools idxstats 补充（sexcheck 等）——非本表收录范围，见对应模块文档。

## 6. MSI（msi，MSIsensor-pro）

仅 tumor 样本有值（node.py 按 paired_design 的 experimental 样本注册 msi 输出），
normal 及未配对样本为 NA。

| 指标 | 含义 | 计算原理 | 反映了什么 / 怎么看 |
|---|---|---|---|
| `msi_score` | 微卫星不稳定分数 | MSIsensor-pro `msi` 输出 tsv 数据行第 3 列（总位点数、可评估位点数、score） | ≥ 2 为 MSI-H 常用阈值（MSIsensor 系惯例，具体以临床指南为准）；tumor-only 无配对 normal 时可信度降低，解读需谨慎 |

---

## 指标登记规范（新增指标时）

1. 在对应解析函数提取该指标，并加入 `build_table()` 的 `core_order`（数值核心项
   排前，状态/计数排后）。
2. 在本文档对应 stage 章节登记一行：指标名（=tsv 列名）、含义、计算原理（原始
   文件哪一行/哪个字段、如何换算）、反映什么/怎么看（含典型范围与阈值）。
3. 需要出图的指标加进 `render_plots()` 的 `--plot` choices 并实现分支。
4. 报告着色阈值变动需同步 `Mutation_report/generate_report.py` 顶部常量并在
   本文档注明阈值出处。
