# TRA / TE-TE 重组与年轻插入验证 —— 方法文档

日期: 2026-11-05 | 文档位置: `modules/sv/docs/`(与代码同库管理)
前置文档: `MEI_METHODS.md`(MEI 管线总览) / `README_polyA_TSD.md`(polyA+TSD 原理)
本文档只记录**方法**: 输入契约、算法步骤、判定规则、统计设计、参数依据、
质控点、方法学教训、局限。不含具体数据结果与运行用法; 结果见输出文件,
运行方式见代码 argparse 帮助。

---

## 0. 问题定义与证据框架

**问题**:
- Q1 (MEI 线): INS 中哪些是转座事件? 是否有新发(体细胞)转座?
- Q2 (TRA 线): 跨染色体易位是否由 TE 驱动(NAHR)? 依据什么分级判定?

**TRA 线证据漏斗(四级, 逐级收紧)**:

| 证据级 | 判定 | TE 驱动? |
|---|---|---|
| T0 事件 | 跨染色体 TRA(镜像去重后) | 未判定 |
| T1 相关 | ≥1 断点落 TE 内 | 否(单侧巧合可能) |
| T2 配对 | 双断点均落 TE 内, 且富集检验显著 | 候选 |
| T3 驱动 | 双拷贝序列同源达 NAHR 机制要求 | **是(序列级)** |

判定口径: 以跨染色体断点(inter-chromosomal TRA)为转座(易位)事件入口;
断点落 TE 注释为"TE 相关"; 两拷贝序列同源达到 NAHR 机制前提才判
"TE 驱动"。单侧 TE、双侧无同源、仅富集未过同源的事件, 各级分列,
不越级归因。TE 驱动的确定性: 双拷贝高同源是 NAHR 的机制前提本身
(重组必须配对同源序列), 统计富集 + 序列同源构成完整证据链;
直接实锤为接合序列本地组装(未系统化, 局限之一)。

---

## 1. 输入契约

| 输入 | 关键字段/约定 | 用途 |
|---|---|---|
| SURVIVOR 合并 VCF | INFO: SVTYPE, CHR2, END, SUPP_VEC, STRANDS, CIPOS/CIEND | TRA/INS 事件来源 |
| 单样本 pbsv VCF | FORMAT=GT:AD:DP; BND 记录含 MATEID、ALT 括号方向 | VAF 验证、链方向恢复 |
| MEI 注释表 (mei_events.tsv) | family/polyA/tsd/context/gene/te_repeat/in_te_region | 事件分类输入 |
| 定年表 (te_divergence_per_ins.tsv) | w_div(Kimura), jc_div, mu 档年龄 | 候选筛选 |
| rmsk TE BED | chrom/start/end/name | 断点 TE 注释 |
| 参考 FASTA(.fai) | — | TE 拷贝序列提取 |
| 基因注释 GTF | gene 层, gene_name/gene_id | 断点基因注释 |
| fiber BAM | HiFi | 接合读长抽查 |

### 1.1 坐标与 ID 约定(易错点)

- VCF POS: **1-based**; INS 插入位于 POS 与 POS+1 之间; BND 的 POS 是
  锚定碱基。
- BED: **0-based half-open**。断点判定: `bisect_right(starts, pos)-1`
  取候选区间, 验 `ends[i] > pos`(TE 长度 >> 1bp 边界误差, 1-based/
  0-based 边界差对判定无实质影响)。
- samtools faidx: **1-based inclusive**; BED→faidx 坐标转换: start+1。
- MEI 表与 ins.fa 的 id 映射: ins.fa 第 i 条序列 = 事件表中超过长度
  阈值的第 i 条(VCF 记录顺序); 需校验 best-hit repeat 名一致(0 错配
  才能继续)。
- SUPP_VEC: 位数=样本数, 对照组在前。读取时 '01' 丢前导零变 '1',
  **必须补零对齐**(已踩过的 bug, 导致一组统计消失)。

### 1.2 工具环境要点

- blastn: BLAST+ 参数是 `-num_threads`, **不是 `-threads`**(传错则
  进程报错退出、stdout 空);
- 容器内临时文件必须写绝对路径(apptainer 只 bind 数据盘, 容器内相对
  路径不可见);
- 外部工具调用必须检查 returncode 与输出非空。

---

## 2. INS 构成分类(转座 vs 非转座)

**方法**: 事件表 family 列计数。family 由 RepeatMasker 分类产出:
插入 ALT 序列超阈值者写 FASTA → RepeatMasker `-s -no_is -lib <te_lib>`
→ 每条取**覆盖跨度最大**的命中(非 SW 最高)→ 按 repeat 名前缀归桶:

| family 桶 | 规则 |
|---|---|
| Alu | 名以 Alu 开头 (SINE) |
| L1 | 名以 L1 开头 (LINE) |
| SVA | 名以 SVA 开头 |
| LTR | 名含 ERV 或以 LTR/MLT/MER/THE 开头 |
| DNA | 名以 Tigger/hAT/Charlie/Zaphod 等前缀开头 |
| simple_repeat / microsatellite | 名含 (X)n 式括号模式 |
| short_unclassified | 长度低于阈值且无命中 |
| other / unclassified | 卫星等其余; 未命中且超阈值 |

**判定逻辑三层区分**(详见 MEI_METHODS.md 第三节): 序列来源层
(RepeatMasker 同源)≠ 移动机制层(TPRT / 整合酶 / 剪切-粘贴)≠
事件层(新生插入 / 历史残留 / 伪影, 需指纹+定年证据分级)。
simple/microsatellite = 串连重复扩张, 非转座。
未分类处理**声明为假设**(假设非 TE → 下界; 全 TE → 上界),
结论以区间表达, 不给单点值。

---

## 3. VAF 验证方法(生殖系多态 vs 体细胞嵌合)

### 3.1 候选筛选(确定性规则)

- 事件表与定年表按 (chrom, pos) 等值 merge(同源 VCF 坐标, 无偏移);
- 候选 = (family ∈ {Alu,L1,SVA} 且 w_div < 1.0) ∪ (context='exon' 且
  w_div < 2.5);
- (chrom, pos) 去重。

### 3.2 原始支持数回查

合并 VCF 丢失 FORMAT(只剩 SUPP_VEC), 必须回查单样本 pbsv VCF:

1. 全量扫描建索引 `CHROM_POS -> FORMAT 字段原文`;
2. 候选 POS 直接查表(合并 POS 与单样本 POS 同源; 若发现偏移改用
   ±1bp 容差);
3. 解析 GT:AD:DP, AD = (ref, alt)。

### 3.3 判定规则

- VAF = alt / (ref + alt);
- 生殖系杂合: GT=0/1 且 VAF≈0.5(抽样噪声下不设硬阈值, 报中位数与
  范围); 纯合: GT=1/1 且 VAF=1.0;
- **体细胞嵌合判据**: VAF 显著低于 0.5(经验 <0.25)且深度足够;
- 质控: 回查的组别归属须与 SUPP_VEC 自洽(三方交叉验证 merge 正确性);
  深度需足以区分 0.5 与低 VAF。

---

## 4. TRA 记录解析与独立事件判定

### 4.1 解析

- 合并 VCF 取 SVTYPE=TRA; INFO 正则 `(\w+)=([^;]+)` 抽 CHR2/END(第二
  断点)/STRANDS/SUPP_VEC;
- 术语: pbsv 单样本报 BND, SURVIVOR 合并转为 TRA——同一事件在管线
  不同阶段的名字。TRA 标签 = inter-chromosomal BND 的合并态;
  intra-BND 属染色体内重排, 不在转座(易位)分析范围。

### 4.2 镜像去重(记录级→事件级)

- 同一接合从两端各报一条(A→B 与 B→A), 去重 key =
  `sorted([(chrom1,pos1),(chrom2,pos2)])`, drop_duplicates;
- 质控: 记录数必须 = 2 × 事件数(每接合两端各一条, 无残缺);
- 局限: 仅合并坐标完全相同的镜像; 近邻断点仍各自计数(见 4.4)。

### 4.3 链信息(SURVIVOR 有损压缩 —— 已实验定性)

迷你 VCF 实验(4 条覆盖全部括号方向 N[X[ / N]X] / ]X]N / [X[N 的 BND
记录跑 SURVIVOR merge): 输出**全部**被改写为 ALT=N[chr2:pos[ +
STRANDS=++; 遍历合并参数组合(type=0/1 × strand=0/1)均如此。

**结论**: SURVIVOR merge **没有任何参数可以保住链信息**。第 5 参数
(strand 0/1)只控制合并时是否按链分组, 不控制输出保真度; TRA 的
ALT/STRANDS 在写出时被硬编码重写, 方向信息确定性丢失。
**可恢复**: 单样本 pbsv VCF 保留原始括号方向(4 种模式齐全), merged
TRA 可按坐标回连单样本 BND, 真实链方向回查原始 ALT 即可恢复;
平衡/不平衡判定需回单样本 VCF + BAM 读深。

### 4.4 复发断点簇(独立性评估)

- 全部断点按 (chrom, pos) 排序, 同染色体相邻间距 >2000bp 处切簇
  (cumsum(diff>2000));
- 簇内 ≥2 断点为复发簇;
- 判读: 簇横跨两组 → 群体复发热点或比对伪影, 非组特异获得性重排;
  近邻簇未合并时, 事件数可能高估。

---

## 5. 断点 TE 富集检验(统计设计)

### 5.1 断点-TE 判定

BED 按 (chrom, start) 排序建 per-chrom 数组, 二分查找判定(见 1.1)。
BED 列有混合类型脏数据, 数值化+去空后使用。

### 5.2 三重对照设计(只报倍数没有说服力, 必须有零分布)

1. **基因组背景**: TE 总碱基 / 基因组总碱基(染色体长度用 max(end)
   近似, T2T 无 gap 时误差 <1%);
2. **随机置换-单断点**: N(默认 1000)次随机断点(染色体按长度加权,
   位置均匀), TE 命中率 → 用于 sanity check(应≈背景, 证明采样无偏);
3. **随机置换-双断点独立期望**: M(默认 2000)对独立随机断点对, 双
   命中率 → 零分布(≈单断点率², 符合独立假设);
4. **二项检验**: H0: 双 TE 率 = 独立期望; H1: >; 单侧 greater。
   置换与二项两种设计独立给出同一结论, 互为交叉验证。

### 5.3 参数表

| 参数 | 默认 | 依据 |
|---|---|---|
| 置换次数(单断点) | 1000 | SE≈1.4%, 标准经验值 |
| 置换次数(双断点对) | 2000 | 同上 |
| 随机种子 | 42 | 可复现 |
| 显著性 | α=0.05 | 常规 |

---

## 6. TE 家族配对与 shuffle 期望

### 6.1 家族归类(前缀规则)

Alu / L1 / SVA / LTR/ERV / LINE2 / MIR / DNA / satellite-simple / other,
按 repeat 名前缀映射。**审核重点**: 手工映射表, 误归类直接影响期望值。
注意点: MacERV* 归 LTR/ERV; Alu 前缀(含 AluYRa1 等)归 Alu;
L1 开头(含 L1_RS2 等)归 L1; 名含 ERV 或 LTR/MLT/MER 前缀归 LTR/ERV。

### 6.2 配对观测

both-TE 事件, 每事件两断点家族组成 pair = sorted(fam1, fam2) 拼接, 计数。

### 6.3 shuffle 期望公式

设 both-TE 断点总数 N(= 2×事件数), 家族 f 频数 c_f。从 N 个断点随机
不放回抽两个为一对, 重复事件数次:

- 跨家族 (a≠b): E = 2·c_a·c_b / (N(N-1)) × 事件数
- 同家族 (a=b): E = c_a(c_a-1) / (N(N-1)) × 事件数

### 6.4 判读规则

- ratio = obs/exp > 1.5 → 家族配对偏好(富集);
- ratio < 0.7 → 回避(低于随机);
- exp < 1 的小类别倍数不稳, 仅参考不作结论。

---

## 7. 序列级同源验证(blastn)

### 7.1 流程(逐对执行)

1. rmsk BED 取两断点各自命中区间;
2. samtools faidx 提取参考序列(坐标转换见 1.1; 临时文件写绝对路径);
3. blastn 双链比对(outfmt 6, pident/length/坐标列);
4. 取最优 HSP: max(pident × length)(匹配碱基数最多的高同源块);
5. 记录 ident(%)、hsp_len(bp)、hom_index = pident × hsp_len / 100。

### 7.2 NAHR 判定阈值与依据

- **NAHR 级**: ident ≥ 90% 且 hsp_len ≥ 100bp。依据: 哺乳动物非等位
  同源重组效率要求 ~95%+ 同源、数百 bp 以上同源臂(经典 L1 近期拷贝对
  值); 90%/100bp 为保守下限, 留测序/比对噪声余量;
- **无同源**: blastn 无任何 HSP(默认 e-value)→ 排除 NAHR, 归
  NHEJ/巧合落点(断点碰巧落 TE——基因组 TE 覆盖高, 随机断裂常落 TE);
- 中间层(80-90%)单列, 不计入 NAHR。

### 7.3 方法学教训(两次错误迭代, 引以为戒)

- v1 difflib SequenceMatcher 一致块: 只匹配正向相同序列, **不搜反向
  互补**——负链 TE 拷贝被误判不同源(错误结论: 同名对多数"不同源");
- v2 blastn 传参 `-threads`(正确为 `-num_threads`), 进程报错退出、
  stdout 空, 全部记 0(错误结论: "全部无同源");
- v3 修正后结果与生物学预期一致。
- 教训: (a) 双序列同源必须用真正比对器且搜双链; (b) 外部调用必须检查
  returncode 与输出非空; (c) **"全部为 0"是系统性错误信号, 不是生物
  学结论**; (d) 已知阳性对照必须在流水线中复现(质控点)。

---

## 8. TRA 全景刻画方法

### 8.1 染色体对网络

cpair = sorted(chrom1, chrom2) 计数; 热点判定按计数排序。

### 8.2 基因注释

- GTF gene 层, per-chrom 按 gene start 排序建 bisect 数组;
- 断点判定: 候选 gene 附近 ±2 行内逐一验 `gene_start <= pos <=
  gene_end`(处理重叠基因); 记 gene_name(无则 gene_id);
- "断点在基因内" = 落 gene body(**不含**启动子窗口, 0 容忍);
- **融合候选** = 双断点均在 gene body 内。

### 8.3 融合候选判读规则

仅断点层面注释; 无 RNA 数据时**不做转录层结论**; 小样本(n<10)下
组间计数失衡不作方向性结论(二项检验确认不显著后仅列候选)。

---

## 9. 统计设计汇总

| 检验 | 位置 | H0 | 统计量 |
|---|---|---|---|
| 置换(单断点) | Stage1 | 断点随机 | N 次 TE 命中率(≈背景的 sanity check) |
| 置换(双断点) | Stage1 | 两断点独立 | M 对双命中率(零分布) |
| 二项(greater) | Stage1 | 双 TE 率 = 独立期望 | k/n 单侧检验 |
| shuffle 期望 | Stage2 | 配对随机 | 家族频数乘积期望 |
| 二项(小样本事后) | 全景 | 两组计数 1:1 | 失衡显著性 |
| MWU/KS | te_age | 组间年龄分布相同 | div 分布差异 |

---

## 10. 质控点清单(复现时逐项核对)

1. 事件表 ↔ ins.fa 一致性: 长度过滤后条数相等, best-hit repeat 名
   0 错配;
2. VAF 回查组别 vs SUPP_VEC 自洽(三方一致);
3. 置换单断点率 ≈ 基因组背景(采样无偏 sanity check);
4. blastn 已知阳性对照(手工验证过的高同源对)在流水线 stage3 复现;
5. 镜像去重计数: 记录数 = 2 × 事件数;
6. SUPP_VEC 补零后三组计数与 VCF 直接 grep 一致。

---

## 11. 局限汇总(按影响排序)

1. **1v1 设计**: "组特异"= 个体私有, 无年龄效应统计力; 一切方向性
   结论仅候选级;
2. **接合序列未组装**: NAHR 判定靠拷贝同源度 + 读长比对签名, 非接合
   序列直接证据; 严格实锤需逐事件本地组装;
3. **STRANDS 信息丢失(可恢复)**: SURVIVOR 写出时确定性丢链信息(无
   参数可解, 见 4.3 实测); 可按坐标回查单样本 VCF 恢复;
4. **独立事件数可能高估**: 近邻断点簇未合并(可能是同一复杂重排的
   多次拆分);
5. VAF 取自 caller 的 AD, 未回 BAM 独立复核跨断点读数;
6. 融合候选无转录层验证;
7. TE 侧翼同源未纳入 NAHR 同源臂计算(hom_index 低估)。
