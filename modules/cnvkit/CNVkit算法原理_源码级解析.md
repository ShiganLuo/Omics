# CNVkit 算法原理:源码级解析

> **配套资源**:`demo_cnvkit_pipeline.py`(可运行演示)、`cnvkit_demo_figure.png`(三面板示意图)、`run_output.txt`(完整运行输出)
> **源码版本**:[etal/cnvkit](https://github.com/etal/cnvkit) master 分支(commit `7539ce3`)
> **定位**:本文是《NGS CNV 检测技术调研与算法分析》报告第 3 章(检测原理)的源码级展开——报告讲"为什么",本文讲"具体怎么做",每个公式都给出对应的函数、文件与行号。所有行号均已对照本机源码核实。

---

## 0. 信号使用总览:深度是主干,BAF 是可选增强,断裂点不使用

报告 3.1 列出了三类信号(RD、BAF、PE/SR)。CNVkit 对三者的使用并不对称,这是理解它在工具图谱中位置的关键:

| 信号 | 是否使用 | 接入环节 | 开启条件 |
|---|---|---|---|
| 深度 RD/log2 | **必需,唯一主干** | coverage→reference→fix→segment→call 全程 | 无(默认即深度驱动) |
| 等位基因 BAF | 可选增强 | ① 分割(联合建模或二次分割) ② call(等位基因特异分解) | `--variants sample.vcf` 且 VCF 含 `alt_freq`(AF)字段 |
| 断裂点 PE/SR | **不使用** | — | —(CNVkit 无 read-pair/split-read 分析) |

BAF 通道的四个接入点(§4.5、§5.3 逐个展开):

1. **CBS/none + VCF** → 分割后对每个 segment 内的镜像 BAF 跑一个 2 态 Viterbi HMM 做二次分割(`variants_in_segment`),BAF 指示的亚段断点被插回 .cns;
2. **haar + VCF** → log2 与 BAF 各跑一次 haar,**断点并集**重建 segment(一次到位,无需二次分割);
3. **hmm-tumor/hmm-germline + VCF** → BAF 以 beta-binomial 发射概率与 log2 高斯**联合建模**,状态空间展开为 (总拷贝数, minor 拷贝数)——三种方式中信息利用最充分的;
4. **call + VCF** → segment 级中位 BAF 参与纯度还原与 `cn1/cn2`(major/minor)分解,使 cnLOH(total CN=2、minor=0)可被显式调用。

不提供 VCF 时,1-4 全部退化为纯深度流程,输出不含 `baf/cn1/cn2` 列——即"CNVkit 只用深度"的说法对默认配置成立,对完整配置不成立。

---

## 0.5 总览:一条命令流水线,五个数学阶段

CNVkit 的肿瘤分析主流程(`cnvkit.py batch` / `import-rna` 之外的经典路径)是:

```text
BAM ──coverage──> .cnn ──reference(多样本)──> reference.cnn
                    │                              │
                    └──────────fix─────────────────┴──> .cnr (bin 级 log2 + weight)
                                                              │
                                                          segment ──> .cns (片段级 log2)
                                                              │
                                                    call (+ VCF 的 BAF) ──> 绝对拷贝数
```

对应到源码,每个箭头是一个独立模块:

| 阶段 | 核心文件 | 输入→输出 | 数学本质 |
|---|---|---|---|
| coverage | `cnvlib/coverage.py` | BAM → `.cnn` | bin 内碱基数 / bin 长度 |
| reference | `cnvlib/reference.py` | N 个正常 `.cnn` → `reference.cnn` | 稳健中心 + 稳健离散度 |
| fix | `cnvlib/fix.py` | `.cnn` + reference → `.cnr` | 偏差回归剔除 + log2 相减 + 加权 |
| segment | `cnvlib/segmentation/` | `.cnr` → `.cns` | 变化点检测 / 隐状态解码 |
| call | `cnvlib/call.py` | `.cns` (+VCF) → 绝对 cn | log2 反变换 + 纯度重缩放 |
| bintest | `cnvlib/bintest.py` | `.cnr` (+`.cns`) → p 值 | z 检验 + BH 多重校正 |

演示脚本 `demo_cnvkit_pipeline.py` 用 400 个合成 bin 的模拟肿瘤样本(p=0.7)走完了全过程,输出显示: CBS 分割将 5 个真实台阶(2→1→2→5→2)全部精确还原,片段 log2 误差 ≤0.06;HMM 自动网格搜索估计出 purity=0.65(真值 0.70)、ploidy=2.0(真值 2.0);纯度重缩放后绝对拷贝数调用命中真实值 2/1/5。

### 0.7 数据结构:CopyNumArray 与 VariantArray

理解源码前先认识两个载体(`skgenome/gary.py::GenomicArray` 是共同基类,包装 pandas DataFrame + 区间运算):

- **`CopyNumArray`(`cnary.py` L91+)**——`.cnr/.cns/.cnn` 的内存形态:必需列 `chromosome/start/end/gene/log2`,可选列 `depth/weight/spread/gc/rmask/baf/cn/cn1/cn2`。它封装了本章反复用到的操作:`center_all`(两层中位数中心化)、`residuals`(对 segment 求残差)、`drop_low_coverage`、`by_arm/by_chromosome/by_gene`(区间分组迭代);
- **`VariantArray`(`vary.py` L13+)**——VCF 的内存形态:必需列 `chromosome/start/end`(SNV 的位置),核心列 `alt_freq`(B 等位基因频率,即 BAF),可选 `alt_count/depth/n_alt_freq`(配对正常样本的 AF)、`zygosity`。

两者按基因组区间对齐的桥接方法是 `vary.py::baf_by_ranges`(L28-73):把落在每个 bin/segment 内的 het SNVs 的 BAF **先镜像到同一侧**(`_mirrored_baf`:0.5+|baf−0.5|),再取中位数——聚合为每区间一个数。镜像这步小而关键:由于 phasing 缺失,某位点观测到 BAF=0.4 还是 0.6 是同一份证据(哪个等位基因为" B"是任意的),不镜像的话聚合中位数会被系统性拉向 0.5,稀释 LOH 信号。配套的 `tumor_boost`(L187+,Bengtsson 2010)在有配对正常 VCF 时进一步把肿瘤 BAF 归一化,放大 LOH 对比度。

---

## 1. 信号从哪里来:coverage(报告 3.1)

### 1.0 深度信号为什么是 CNVkit 的唯一主干

三个事实决定了这一架构选择(源头是 2016 年 PLoS Comput Biol 的 CNVkit 论文,设计目标是 panel/WES 的 RD 检测):

1. **RD 是 panel/WES 唯一可用的信号**:外显子组捕获读长小(2×150bp)、插入片段短,断裂点 read-pair 证据在捕获边界外几乎不存在;BAF 需要 VCF 步骤(突变 calling),而胚系 het SNP 密度在 panel 上往往只有几十~几百个,稀疏;
2. **RD 的统计性质好**:Poisson 计数→log2 近正态、方差稳定,整个加性统计框架(reference 池化、分割、z 检验)都建立在它上面;
3. **antitarget 补齐了分辨率**:CNVkit 用"靶内密 bin + 靶外疏 bin"让 WGS 也获得近均匀的深度分辨率,不必依赖 SV 类证据(对比:WGS 专用工具才主打 PE/SR)。

代价同样明确:纯深度无法区分 cnLOH(cn=2 但杂合性丢失)与正常二倍体——深度完全相同;也无力解析断裂点位置。这正是 §4.5/§5.3 的 BAF 通道存在的意义:BAF 是等位基因层的信息,与深度正交,二者联合才能拆出 major/minor 拷贝。

### 1.1 深度信号的原始定义

CNVkit 对"某个 bin 的深度"的定义朴素而精确(`coverage.py::region_depth_count`,L510-594):

```python
# 只统计通过过滤的 reads (重复/二次比对/未比对/QC失败/MAPQ不足, L536-548)
def filter_read(read) -> bool:
    return not (read.is_duplicate or read.is_secondary
                or read.is_unmapped or read.is_qcfail
                or read.mapq < min_mapq)

# 逐 read 累加"落在区间内的比对碱基数", 最后:
depth = bases / (end - start) if end > start else 0
```

三个值得注意的工程细节,都会直接影响 CNV 信号质量:

1. **过滤旗标**与 `samtools depth` 默认一致(0x704:dup/secondary/unmapped/qcfail),再加 MAPQ 下限(`--min-mapq`,默认 0,生产建议 ≥20);
2. **`--no-overlap` 选项**(L524-533):短插入片段(典型如 FFPE 样本)的两条 mates 会重叠,按 read 数计会把同一片段算两遍。开启后按 template 合并比对块(区间合并算法,L575-586),即 `samtools depth -s` 语义;
3. 默认走 **`depth`(碱基数)而非 `count`(reads 数)**:`by_count=False`。reads 数对 bin 长度敏感(长 bin 天然捕到更多 reads),碱基数才是"每碱基期望深度"的无偏估计。

### 1.2 log2 比值:depth 信号的最终形态

coverage 命令输出的 `.cnn` 中,每个 bin 带有 `log2` 列。对肿瘤样本,它此刻还是"相对预期深度的对数";真正的"比值"要等 fix 阶段减去参考后才成立。log2 是整个 CNVkit 的统一信号尺度——原因有二:

- 拷贝数是倍乘关系(cn=2→1 是"减半",2→4 是"加倍"),对数尺度下 gain/loss 对称(log2=+1 与 -1);
- 深度噪声近似 Poisson→正态,取对数后方差近似恒定,后续的加性模型(分割、检验)才成立。

### 1.3 bin 的设计:target 与 antitarget(这是 CNVkit 的招牌)

`target.py::do_target`(L36-105)做三件事:丢弃零宽区间、**合并重叠/相邻 bait**(L67,`merge(bp=1)`),可选把大区间均分为 ~267 bp 的小 bin(`subdivide`,默认 `avg_size=200/0.75`)。

`antitarget.py` 则在靶区之外的"安静区"均匀撒 bin(避开高 GC、高重复区)。**这套 target+antitarget 双尺度 bin 是 CNVkit 区别于纯 WES 工具的核心设计**:WES 数据只有外显子区的均匀覆盖,CNVkit 用 antitarget bins 让 WGS 也能低成本获得全基因组均匀分辨率——靶内 bin 密、靶外 bin 疏,深度分辨率随测序设计自适应。

---

## 2. 参考基线:reference(报告 3.2 "PON 是生命线")

### 2.1 为什么需要 reference:乘法模型与可辨识性

单个 bin 的观测深度 = f(拷贝数) × g(GC 含量) × h(靶区长度/可映射性) × m(该样本总测序量) + 噪声。这个乘法结构不是建模约定,而是计数过程的代数:一条 read 对应一个模板分子,分子数 = locus 模板数(∝ cn × 捕获效率 × 留存率)× 抽样份额(∝ 测序量)。每一级都是"比例",比例复合只能是乘法。三个可证伪的定量预言均被本仓库测试数据证实(`test_multiplicative.py`/`test_ratio.py`,CNVkit 官方测试集):

1. 同 panel 两样本 log2 相减后 GC 因子抵消(样本对 Δlog2 vs GC 残差 robust sd 仅 0.20);
2. 女性chrX(2拷贝)/男性(1拷贝)的 log2 差 = log2(2) = **+1.000**——无参数硬预言,实测 +1.009 与 +0.998;
3. 逐 bin 深度比 = 2^Δlog2(33,840 bins 相关 = 1.0000)。

取对数后乘法变加法,每个 bin 的观测是三成分之和:

```text
y_i  =  tv(cn_i)  +  b(GC_i)  +  e_i
        ↑拷贝数信号   ↑系统偏差    ↑噪声
```

基线校正的目标是扣掉 b、留下 tv。**核心定理:当且仅当基线的估计来自"不含 tv 的数据"(或 tv 与协变量独立)时,扣除才干净。**这一可辨识性条件决定了 CNVkit 的整个架构:必须有一批正常样本来充当"信号-free 的基线源"——这就是 reference 存在的数学理由,报告 3.2 "PON 是生命线"说的就是它(违反该条件的后果与三种路线的定量对照见 §2.6)。

CNVkit 的思路:**用一批正常样本构建每个 bin 的"期望深度 + 离散度",之后任何样本都相对这个基线取比值**。

### 2.2 池化的数学:不用均值,用 biweight

`reference.py::combine_probes`(L385-494)→ `load_sample_block`(L497-606)把 N 个正常样本的 log2 堆成矩阵 `(bins × samples)`,然后 `summarize_info`(L717-741)逐 bin 计算:

```python
cvg_centers = np.apply_along_axis(descriptives.biweight_location, 0, all_logr)
spreads = np.array([
    descriptives.biweight_midvariance(a, initial=i)
    for a, i in zip(all_logr.T, cvg_centers, strict=True)])
```

**biweight(双权)估计量**是这里的关键选型。以 location 为例(`descriptives.py::biweight_location`,L96-142):

```python
d = a - initial                     # 残差
scale = c * mad                     # c=6, MAD=中位绝对偏差
w = (1 - (d / scale) ** 2) ** 2     # 离中心越远权重越小, |d|>scale 的点权重=0
initial += (d[mask] * w[mask]).sum() / weightsum   # 迭代至收敛 (max 5 轮)
```

对比算术均值:均值对离群值是线性敏感的——假如 20 个"正常对照"里混入 1 个真实携带 3 Mb 缺失的样本(临床队列常见),算术均值会把该 bin 的基线整体拉低,下游所有样本在该 bin 出现假性 gain。biweight 给这个离群样本的权重趋近 0,基线几乎不动。scale 估计 `biweight_midvariance`(L188-232,c=9)同理稳健,产出的 `spread` 列后续直接决定每个 bin 的权重(§3.4)。

### 2.3 reference 中还存了什么

`load_sample_block` 同时为每个 bin 计算并缓存(`ref_columns`,L546-561):

- `gc`:该 bin 区间在参考基因组上的 GC 含量(`get_fasta_stats`,FASTA 逐区间提取);
- `rmask`:RepeatMasker 注释的重复序列占比;
- `depth`:池化的深度中心;
- `spread`:池化离散度(≈bin 间 sd)。

`gc` 与 `rmask` 是 fix 阶段偏差校正的"排序键"(§3.3),`spread` 是加权依据(§3.4)。**reference.cnn 因此不只是"平均深度表",而是把后续校正所需的全部协变量一并封装**——这是"reference 一旦建好,单样本处理只跑 fix"高效性的来源。

另有 `--cluster` 选项(`create_clusters`,L744+):对参考样本按 log2 相关性做层次聚类,为每个簇单独存 `log2_i/spread_i` 列;fix 时(`fix.py` L179-208)按样本与各簇的相关系数选最优簇。应用场景:参考池混有不同plate/建库方法批次时,避免跨批次的系统差。

### 2.4 性别染色体的预处理

多样本池化前,男女样本的 chrX 深度天然差一倍,必须统一。`shift_sex_chroms`(L670-714)按样本性别把 X/Y 的 log2 平移到参考性别基准(代码注释里的真值表:XY 参考 chrX=-1,XX 参考 chrX=0),性别由 `verify_sex`/chrX 相对残差推断(`_guess_sexes`),未知样本默认按女性处理并告警(L701-707)。PAR(拟常染色体区)按二倍体处理,坐标来自 `skgenome.genomebuild` 的硬编码表。

### 2.5 PON 构建自身的校正链

正常样本在池化前也要走一遍偏差校正(`reference.py::bias_correct_logr`,L625-667):`center_all` → 性别平移 → `center_by_window` 按 GC/rmask/edge 逐层扣除。**这一步用 M2(滑窗中位数)是安全的,因为正常样本没有 CNV**——"窗内信号均值为零"的假设自动成立(可辨识性条件在这里天然满足)。正常样本完成校正后才堆叠池化。换言之:reference 的干净,是"信号-free 数据 + 非参数滑窗 + 稳健池化"三层叠加的结果,任何一层缺失都会让偏差泄漏进基线。

### 2.6 反面教材:三条基线路线的定量对照

如果不用 PON,改用"样本自身回归掉 GC"会怎样?`test_baseline_regression.py`(2000 bins,GC 偏差 2.3 log2,真实台阶 cn=2/1/4、p=0.7)给出定量答案。两种情境:情境 A(cn 与 GC 独立)、情境 B(缺失区故意放在 GC 最高端——真实对应 GC 富集启动子/脆性位点更易断裂,肿瘤 CNV 与 GC 并不独立):

| 路线 | 基线来源 | 情境 A 缺失区误差 | 情境 B 缺失区误差 |
|---|---|---|---|
| 不校正 | — | −0.02 | −1.35 |
| M1 参数回归 y~poly(GC) | 样本自身 | −0.03 | **+0.29(吃掉 47% 信号)** |
| M2 滑窗中位数(CNVkit 式) | 样本自身 | −0.00 | **+0.57(吃掉 92% 信号)** |
| M3 PON 相减 | 30 个正常样本 | +0.01 | **−0.00 ✓** |

- **M1 必错(共线性)**:回归 `y ~ poly(GC)` 时,y 中与 GC 共线的成分是 `b(GC) + tv(cn)·corr(cn, GC)`——偏差与信号混在同一协变量方向(情境 B 中 corr=+0.21),最小二乘无从区分,只能一起回归掉。CNVkit 全程刻意避免对样本自身 log2 做参数化 GC 回归。
- **M2 是受控妥协**:滑窗中位数隐含"窗内 CNV 均值为零"。§3.3 的"先随机打散再按 GC 排序"正是为了在 cn⊥GC 时维持它;情境 B 中假设破产,窗内均值≠0 → 偏差被高估、信号被部分吸收。所以 CNVkit 只把 M2 用在信号-free 或信号已稀释的位置(正常样本校正 §2.5、PON 之后的残余偏差 §3.3)。
- **M3 是金标准**:正常样本 tv=0,biweight 中位数是 b 的无偏且信号-free 估计,相减后 b 严格抵消。**PON 不是经验主义,它是可辨识性的唯一来源。**

**实战推论**:PON 必须用匹配的正常人样本,**不能用肿瘤样本混池**——混池平均携带大量 CNV,用它当基线 = 基线里混入信号,所有样本的 calling 都会向"混池的平均肿瘤图谱"收缩(系统性漏检混池中常见的事件)。PON 样本量与代表性比算法选择更影响最终灵敏度:这就是报告 4.1 质控要点中"PON 质量是 panel 分析第一质控项"的统计学根据。

---

## 3. 核心校正:fix(报告 3.2 标准化的源码主体)

`fix.py::do_fix`(L20-215)是 CNVkit 单样本处理的中枢,四个子步骤顺序执行。

### 3.1 坏 bin 过滤:`mask_bad_bins`

```python
mask = (cnarr["log2"] < params.MIN_REF_COVERAGE        # -5.0: 参考中近零覆盖
      | (cnarr["spread"] > params.MAX_REF_SPREAD)       # 1.0: 参考样本间分歧过大
      | (cnarr["gc"] > 0.7) | (cnarr["gc"] < 0.3))      # GC 极端区
```
(参数见 `params.py` L4-10;`fix.py` L295-321)

这一步在"减参考"**之前**执行:参考中不可靠的 bin(如 para-unique 区、极端 GC 区)不管样本测得多好,其 log2 都没有可比性,直接丢弃比加权保留更干净。demo 中该规则丢弃了 31/400 个 GC 极端的模拟 bin。

### 3.2 总量中心化:`center_all`

`cnary.py::center_all`(L282-342)把全基因组 log2 的中心平移到 0。默认估计器是 median,且 `by_chrom=True`:先对每条染色体取 median,再对这些 per-chromosome median 取 median。两层 median 的目的:**高度非整倍体的肿瘤(如三倍体、chr7 扩增)若直接全基因组 median,会被真实 CNV 拖偏**;两层结构假设"多数染色体处于同一倍性水平",对真实肿瘤更成立。`skip_low=True` 时剔除近零覆盖 bin 再估计,防止 0 深度 bin(对应 `NULL_LOG2_COVERAGE=-20` 哨兵值)污染中心。

### 3.3 偏差校正:`center_by_window` —— CNVkit 最"巧"的一段

GC/重复序列/edge-effect 偏差都不是简单线性关系(GC 通常呈倒 U 形),CNVkit 不拟合参数化曲线,而是**非参数的"沿协变量排序 + 滑动中位数 + 相减"**(`fix.py::center_by_window`,L373-427):

```python
rng = np.random.default_rng(0xA5EED)
shuffle_order = rng.permutation(df.index)   # (1) 先随机打乱所有 bins
df = df.iloc[shuffle_order]
sort_key = sort_key[shuffle_order]
order = np.argsort(sort_key, kind="mergesort")
df = df.iloc[order]                          # (2) 按协变量(GC 等)排序
biases = smoothing.rolling_median(df["log2"], fraction)   # (3) 滑窗中位数
df["log2"] -= biases                         # (4) 减去局部趋势
```

逐步拆解:

- **为什么先打乱?** 真实 CNV 在基因组上是连续的区段(相邻 bin 拷贝数相关)。若直接按基因组顺序滑窗,窗口内 CNV 信号占主导;先随机打散再按 GC 排序,窗口内的 bins 来自基因组各处,**它们共享的是 GC 偏差(协变量),而 CNV 信号在窗口内近似均值为零**——滑动中位数估计的就是纯偏差曲线;
- **窗口多大?** `fraction = max(0.01, len(cnarr) ** -0.5)`(L263),400 bins → 20 bins。窗口是"偏差拟合精度 vs CNV 信号保留"的调节旋钮:窗口越小,GC 曲线贴得越紧,但一个恰好落进窗口的宽 CNV 会被部分吸收;窗口越大反之。CNVkit 用 √n 规则让窗口随数据量自动缩放;
- **为什么用中位数?** 滑窗内即使混入少数 CNV bins,中位数(而非均值)不被 ≤50% 的污染拉动——又一处稳健统计选型;
- 同一框架复用于三种偏差:GC(`sort_key=ref['gc']`)、RepeatMasker(`ref['rmask']`,仅 antitarget)、edge effect(`get_edge_bias` 计算的解析量)。**顺序**:GC → edge → rmask 依次执行,每次减去一层偏差。

**edge effect 的解析模型**(`fix.py::edge_losses/edge_gains`,L462-512)值得一提:杂交捕获的探针只能"拉"到插入片段长度范围内的片段,靶区两端各损失约 `i/(2t)`(i=insert size=250,t=靶长)的覆盖;相邻靶区间距 <i 时,彼此的侧翼 reads 又会互相"溢出"增益 `(i-g)²/(4it)`(g=间距)。CNVkit 把损失与增益解析算出、合成一个数值作为该靶点的 edge 排序键——这不是经验回归,而是从片段-探针物理模型推出的期望覆盖形状。

demo 输出中,GC 十分位组均值轨迹极差从 0.245(校正前)压到 0.173(校正后),同时 5 个真实台阶保持完好——这正是"窗口内 CNV 均值为零"假设成立时的预期行为。

### 3.4 参考相减与加权:`apply_weights`

偏差校正后,`do_fix` 的最后两步(L212-214):

```python
cnarr.data["log2"] -= ref_matched[log2_key]        # log2比值 = 样本 - 参考
cnarr = apply_weights(cnarr, ref_matched, ...)
cnarr.center_all(skip_low=True, ...)               # 最终中心化
```

至此每个 bin 的 log2 才是报告 3.1 所说的"深度比值信号"。

`apply_weights`(L515-617)为每个 bin 计算 `weight ∈ (0, 1]`,思路是 **weight ≈ 1 − 方差**:

- **样本内项**(`simple_wt`):分别对 target/antitarget bins 计算全基因组残差方差(对 segment 中位数的残差做 `biweight_midvariance`),再按 bin 大小归一:`w = 1 − var / (bin_sz / mean_sz)`——大 bin 深度更稳、权重更高;
- **参考项**(`fancy_wt`,仅池化参考时):`w = 1 − spread²`,参考样本间分歧大的 bin(如存在隐性结构变异的区段)降权;
- 混合:`weights = 0.9*fancy + 0.1*simple`(L605-606),NaN 兜底为 `epsilon=1e-4` 并 `clip(1e-4, 1.0)`。

这些权重随后贯穿整个下游:分割时加权平均(CBS 的 `weights=`、haar 的 weight 通道、HMM 的 sd 估计)、`bintest` 的 z 分母、`.cns` 的 segment 均值加权。**"每个 bin 有一个质量分数"是 CNVkit 输出格式的关键设计**——不是所有 bins 生而平等,低质量 bin 参与统计但话语权小。

> **本节方法为什么这样设计**:§3.2/§3.3 用样本自身做基线(中心化、滑窗中位数),其合法性来自 §2.6 的可辨识性条件——它们只被用在"PON 已扣除主偏差、残余偏差且信号已稀释"的位置;真正的信号-free 基线是 §2.5 的 PON。三条路线的定量对照见 §2.6。

---

## 4. 分割:从 bin 级噪声到片段级台阶(报告 3.3)

输入:369 个 bins 的加权 log2(demo 实际有效 bin 数)。输出:`.cns`,每行一个 segment(区间 + 加权均值 log2 + probes 数)。CNVkit 提供 5 种方法:`cbs / haar / hmm / hmm-tumor / hmm-germline`(`segmentation/__init__.py::do_segmentation`,L37-201;CBS/haar 按染色体臂并行,HMM 全基因组单模型)。公共预处理(L272-332):丢非有限 log2 → 可选 `skip_low` → `drop_outliers`(滑动分位数 10×MAD 外的散点)→ 丢零权重 bins。

### 4.1 CBS:统计检验驱动的循环二分(默认方法)

`segmentation/cbs.py` 内嵌完整 R 脚本(L1-78)调用 Bioconductor **DNAcopy** 包,核心一行:

```r
fit = segment(cna, weights=tbl$weight, alpha=threshold)   # threshold 默认 1e-4
```

算法骨架:对当前区间,枚举所有二元切分点 t,统计量 `max |Σ_{i≤t} w_i(x_i−x̄)| / sqrt(Σw²·p(1−p))`(标准化加权 CUSUM);置换检验(重排 bin 顺序)给出 p 值;p < α 则在极值点切开、递归两侧;否则停止。**每个新 fragment 的诞生都经过一次显式假设检验**,α=1e-4 直接控制假片段率。R 脚本收尾处(L64-73)还有个常被忽视的细节:用原始 bin 边界重写 segment 端点、用加权均值重算 `seg.mean`——保证 .cns 与 .cnr 完全对齐、均值可复现。

demo 中 CBS 恰好输出 5 个 segment,与真实结构 2/1/2/5/2 一一对应,端点误差 ≤1 bin。

### 4.2 HaarSeg:小波多尺度 + FDR(纯 Python、免 R 依赖)

`segmentation/haar.py`(移植自 HaarSeg R/C 原实现)。核心思想:**拷贝数台阶在不同尺度下都表现为"阶跃",Haar 小波对阶跃的响应是局部极值**——用一族宽度指数递增的 Haar 小波(2^1 到 2^5 个 bins)扫过信号,每个尺度独立找峰,再跨尺度合并(`UnifyLevels`,L615-663:细尺度断点若落在粗尺度断点 ±window 内则不重复记录)。

显著性控制(`FDRThres`,L374-402):对每个尺度上的峰振幅做双侧正态检验,p 值经 **Benjamini-Hochberg FDR** 校正(`p <= m*q` 的最大秩),q=`fdr_q` 默认 1e-4。噪声 σ 用最细尺度响应的 **MAD×1.4826** 稳健估计(L308-333,含 σ=0 时的 SIGMA_FLOOR 兜底)。最终 `SegmentByPeaks`(L405-451)按断点分组计算**加权均值**。

提供 VCF 时,`one_chrom_baf`(L188+)对 log2 与逐 bin BAF **各跑一次 haar、取断点并集**后再重建 segment——专门捕捉"深度平坦、BAF 偏移"的拷贝中性 LOH(cnLOH),这是纯深度分割在原理上不可能发现的事件。

一个反直觉的实现事实(`_haar_signal`,L101-114 的长注释):haar 输入必须是**未平滑的原始 log2**——FDR 阈值以"最细尺度 MAD≈真实噪声"为前提,若先做 Savitzky-Golay 平滑,噪声被人为压低,阈值塌陷导致 ~8× 过分割。demo 中 haar 出现 1 例过分割(170-180kb 从中性段切出),是 FDR 1e-4 下的单峰误报,恰好演示了"更灵敏的方法需要更严的 q 值"。

### 4.3 HMM:拷贝数作为隐状态(源码自带纯 NumPy 实现)

`segmentation/hmm.py`(968 行,无 pomegranate 依赖)。建模五要素:

**状态空间**(`enumerate_states`,L36-60):`(总拷贝数 cn, minor 等位拷贝)` 二元组,`max_copies=5` 且启用 BAF 时为 (0,0),(1,0),(2,0),(2,1),(3,0),(3,1),(4,0)…——**正是报告 3.4 所说"把每段拷贝拆为 major/minor 两条染色体"的离散化**;不提供 VCF 时 BAF 分量缺失,状态退化为纯 cn。

**发射概率**(`log_emission_probs`,L135-196):log2 分量是高斯 `N(obs | E[log2], σ)`;BAF 分量是 beta-binomial(聚合 allele counts,ρ=200 吸收过离散)。两个期望值函数就是报告里的混合公式:

```python
# _expected_log2 (L66-83): E[log2] = log2( (p·cn + (1−p)·2) / (p·ploidy + (1−p)·2) )
# _expected_minor_freq (L86-104): E[BAF] = (p·minor + (1−p)·1) / (p·cn + (1−p)·2)
```

注意 E[BAF] 分子里的 `1`:正常细胞每条染色体恰 1 份 B 等位;cn=0 时约定返回 0.5(无等位信息)。**联合发射的威力**在于二维似然可分离纯深度无法区分的状态:cn=2/minor=1(正常)与 cn=2/minor=0(cnLOH)的 log2 期望完全相同,但 BAF 期望分别为 0.5 与 1.0——Viterbi 依据 BAF 观测即可分开二者;这也是报告 3.1 "BAF 是 cnLOH 的关键证据"的模型化表述。`_switch_weights`(L202-226)按 |Δcn| 软约束转移(见下),同 cn 不同 minor 的状态间切换代价低,符合"等位失衡常后于拷贝数改变发生"的生物学直觉。

**转移概率**(`log_transition_matrices`,L229-271):与距离相关——`p_stay = exp(−d/L)`,L=5 Mb:相邻 bin 越远越容易发生状态切换(着丝粒/间隙处断点先验更高);切换幅度另有惩罚 `exp(−8·|Δcn|)`:cn2→3 比 cn2→0 便宜得多,编码"CNV 是局部事件、大跳罕见"的先验;

**解码**(`viterbi_log`,L320-360):标准 log 空间 Viterbi,`δ_t(j) = e_t(j) + max_i [δ_{t−1}(i) + a_ij]`,回溯得最优状态路径;

**纯度/倍性自动估计**(`grid_search_purity_ploidy`,L575-668):对 (purity, ploidy) 网格逐点跑 forward 算法算边际似然,取 argmax。demo 中估计出 purity=0.65 / ploidy=2.0(真值 0.70/2.0)——**报告 3.4 说的"在 LRR 平面上联合拟合纯度与倍性",在这里是一个显式的网格最大似然**。

三种方法输出对比(demo 实测):

| 方法 | 段数 | 断点 | 备注 |
|---|---|---|---|
| CBS | 5 | 全部正确 | 每段经置换检验,最保守 |
| haar | 6 | 4 正确 + 1 过分割 | 快(~ms 级),免 R |
| HMM | 8 | 结构正确,扩增段两端有锯齿 | 直接输出整数状态,自动估纯度 |

选择建议(与 CNVkit 文档一致):肿瘤 WES/panel 默认 **CBS**;大批量 WGS 或无 R 环境用 **haar**;需要等位基因特异输出(major/minor)或联合 BAF 分割时用 **hmm-tumor**(胚系用 `hmm-germline`,纯度固定 1.0)。

### 4.4 bintest:bin 级显著性(报告 3.3 的"小事件"检验)

分割解决"台阶在哪",但**台阶内可能还有更小的独立事件**(如 1 个外显子的缺失被大片段增益淹没)。`bintest.py::do_bintest`(L20-95):残差 = bin log2 − 所属 segment 均值(`cnary.py::residuals`,L637-679),标准差从权重反推 `sd = sqrt(1 − weight)`(L98-105),z = 残差/sd,双侧正态 p 值经 **BH 多重校正**(L116-123)。

demo 的对比实验很有教学意义:**同一份数据,零假设取"所属 segment 均值"时 0/369 显著(分割已解释一切),取"log2=0(全基因组中性)"时 80/369 显著且精确聚焦在缺失/扩增区内**。这演示了假设检验的一个通用原则:零假设的选取决定了"什么算异常"。

### 4.5 BAF 进入分割的三条路径(按方法分派)

`segment --variants sample.vcf` 时,BAF 如何接入取决于所选分割方法(`segmentation/__init__.py` L391-416 的分派逻辑)——三条路径的"联合程度"递增:

**路径 A:CBS / none → 分割后二次分割(后验修正)**

```python
# segmentation/__init__.py L397-404
newsegs = [hmm.variants_in_segment(subvarr, segment)
           for segment, subvarr in variants.by_ranges(segarr)]
```

先纯深度跑完 CBS,再对每个 segment 内部用**镜像 BAF 上的 2 态 HMM**(`hmm.py::variants_in_segment`,L864+;状态=中性 BAF~0.5 vs 偏离 BAF~0.67,p_stay=0.99,最少需 50 个 SNVs)检测等位失衡亚段;亚段断点取相邻 SNV 的中点插回 .cns。逻辑:深度方法对"拷贝数不变、只有等位比改变"的区段天然失明,BAF 二次分割负责补救——这是 cnLOH 在 CBS 路径下被发现的唯一机会。

**路径 B:haar → 断点并集(一步联合)**

`one_chrom_baf`(haar.py L188-226)对 log2 与逐 bin 镜像 BAF **各跑一次 haar,取断点并集**重建 segment。比路径 A 简单,且断点在第一轮就同时受两种信号约束;代价是两种信号没有各自的权重/噪声模型,联合是几何的而非概率的。

**路径 C:hmm-tumor / hmm-germline → 完全联合建模**

如 §4.3,BAF 的 beta-binomial 发射与 log2 高斯发射在**同一状态的似然中相加**,Viterbi 同时用两种观测解码一条状态路径。状态空间直接以 (cn, minor) 为坐标——不需要"二次分割"这种补丁式设计,也不存在"先深度后 BAF"的顺序偏置。信息利用最充分,计算也最重(要网格搜索 purity/ploidy)。

| 路径 | 方法 | 联合时机 | cnLOH 检出 | 适用 |
|---|---|---|---|---|
| A | cbs/none | 分割后(2态HMM 修正) | ✓(亚段级) | 默认选择,稳健 |
| B | haar | 断点并集 | ✓(区段级) | 快、免 R |
| C | hmm-tumor/germline | 发射概率级 | ✓(状态级,最细) | WGS + 足量 SNVs |

三条路径殊途同归的架构含义:**CNVkit 把"深度主干 + BAF 增强"做成了正交的可选层**——不供 VCF 时它是最纯粹的 RD 工具,供 VCF 时获得等位基因特异能力,但始终不碰 read-level 的断裂点证据(那是 SV 工具的领地,报告 3.1 也建议与其联用)。

---

## 5. 绝对拷贝数与纯度:call(报告 3.4 混合去卷积)

### 5.1 log2 → 绝对拷贝数:一个反解 + 一个重缩放

`call.py::do_call`(L33+)支持三种方法。核心公式在 `_log2_ratio_to_absolute`(L532-579),注释里写全了推导:

```text
v = log2( p·n/r + (1−p)·x/r )        # 观测: 纯度p, 肿瘤拷贝n, 参考倍性r, 期望拷贝x
2^v = p·n/r + (1−p)·x/r
n = ( r·2^v − x·(1−p) ) / p          # 纯度重缩放 (p<1 时)
n = r·2^v                            # 纯样本/胚系 (p=1)
```

注意这与报告 3.1 的期望信号公式是同一个混合模型的正反两面:**fix 之前它用来"预测观测",call 之中它用来"反解真相"**。demo 的数值验证:cn=1 片段在 p=0.7 下的观测 v=−0.621,代入公式还原 n=1.00——同时 `do_call` 对 n<0 做了 clamp(纯合缺失在正常细胞污染下推到负值,L564-576,issue #503)。

BAF 侧对称地有 `rescale_baf`(L594-624):`t_baf = (obs_baf − 0.5·(1−p)) / p`,把正常细胞把 BAF 拉向 0.5 的"稀释"反推回去;还原后超出 [0,1] 的值说明纯度估计与观测矛盾(典型原因:VCF 里混入了 somatic 而非 germline het 位点),clamp 并告警。

### 5.2 两种调用哲学:threshold vs clonal

- **`threshold` 法**(默认,`absolute_threshold`,L297-353):4 个 log2 阈值 `(−1.1, −0.25, +0.2, +0.7)` 把连续 log2 切成 5 档 cn=0/1/2/3/4+。阈值的来历写在 docstring 里:假设 cn∈{0..4}、50% 肿瘤细胞,期望 log2 = log2(2:6/4) ≈ (−1.0, −0.415, 0, +0.322, +0.585),再放 ±0.1 噪声余量取中点。**它隐含假设了 p=0.5**;
- **`clonal` 法**:显式输入 purity(可来自 `--diagram --tumor-boost` 的估计或外部工具如 ASCAT/FACETS),按 §5.1 公式重缩放,输出连续绝对拷贝数(可 `--round-result` 取整)。

demo 中同一份 .cns:threshold 法把扩增段调为 cn=4(因为观测 log2=+0.97 > +0.7 只够到"AMP"档),clonal 法(p=0.7)调为 cn=5——与真实值一致。**这就是报告 3.4 的核心论点的源码体现:不校正纯度,拷贝数系统性低估;纯度越低,高拷贝事件被压得越扁**。

### 5.3 等位基因特异输出:从 VCF 到 cn1/cn2 的完整链路

`do_call(variants=VCF)` 触发一条四步流水(call.py L162-236),这是 BAF 通道的最终汇合点:

**① 聚合**:每个 segment 内 het SNVs 的 BAF 镜像后取中位数 → `baf` 列;同时数出支持 SNP 个数 → `nbaf` 列(L162-172)。`nbaf` 不是冗余信息:`export vcf` 会把它写成 BAFN INFO 字段,让下游按证据量加权 LOH 判断(一个 nbaf=3 的 BAF 与 nbaf=60 的 BAF 不该同权重)。

**② 纯度还原**:`rescale_baf`(L594-624)把正常细胞稀释效应反推掉:

```text
t_baf = (obs_baf − 0.5·(1−p)) / p
```

推导与 §5.1 的 log2 重缩放同源:观测 BAF 是肿瘤 BAF 与正常 BAF(=0.5)按纯度的线性混合。还原后落在 [0,1] 之外的值是"纯度估计与观测矛盾"的哨兵(典型原因:VCF 混入了 somatic 变异而非胚系 het),会被 clamp 并告警——**这也是 purity 参数的质量校验器**。

**③ major/minor 分解**(L226-234),三行代码:

```python
upper_baf = (outarr["baf"] − 0.5).abs() + 0.5          # BAF 折回 ≥0.5 半边
outarr["cn1"] = (absolutes * upper_baf).round().clip(0, outarr["cn"])
outarr["cn2"] = outarr["cn"] − outarr["cn1"]
```

原理:肿瘤 BAF(还原后)≈ major 拷贝占比,故 `cn1 = round(cn × upper_baf)`、`cn2 = cn − cn1`。数值验证:cn=3,baf=0.67 → cn1=2, cn2=1(AAB ✓);cn=2,baf=0.99 → cn1=2, cn2=0(**cnLOH 显式现形**——纯深度看 cn=2 是"正常",BAF 揭示两条同源);cn=1,baf=0 → cn1=1, cn2=0(半合缺失)。`baf` 为 NaN(该段无 het SNVs)时 cn1/cn2 置 NaN 而非猜值。

**④ BAF 在 segmentation 阶段的去向**:若走路径 A/B(§4.5),.cns 已带 BAF 修正过的断点;若走路径 C,`baf` 直接来自联合模型。call 不再重算分割,只在既定 segment 上聚合。

**能力边界**:这条链路只利用**胚系 het SNPs**(配对正常样本的 `n_alt_freq` 可选启用 TumorBoost 增强)。肿瘤纯合突变、体细胞 SNV 的 VAF 不进入模型——体细胞等位比异常(如早期克隆的 LOH 残留)需要 ASCAT/FACETS 类工具或专门的 subclone 分析。panel 数据 het SNP 太稀(nbaf 普遍 <10)时,cn1/cn2 的置信度有限,报告时应谨慎使用。

---

## 6. 与报告第 3 章论断的逐条对照

| 报告论断 | 源码落点 | demo 验证 |
|---|---|---|
| 3.1 "期望信号为 log2[(p·C_t+(1−p)·2)/2]" | `hmm.py::_expected_log2` L66-83;`call.py::_log2_ratio_to_absolute` L532(反解) | p=0.7 下 cn=1 期望 −0.621,实测 −0.622;反解还原 n=1.00 |
| 3.2 "按 GC/靶区特征的样内校正" | `fix.py::center_by_window` L373-427;`edge_losses/gains` L462-512 | GC 分箱轨迹极差 0.245→0.173,台阶保留 |
| 3.2 "基于对照样本集的基线建模" | `reference.py::summarize_info` L717-741(biweight) | 5 正常样本池化,离群免疫 |
| 3.3 "CBS:统计检验支持下迭代寻找均值变化点" | `cbs.py` L60(`DNAcopy::segment(alpha=1e-4)`) | 5/5 台阶精确还原 |
| 3.3 "置换/自助检验评估分段显著性" | CBS 置换检验;`haar.py::FDRThres` L374-402(BH-FDR) | haar 过分割 1 例,直观展示 q 值的作用 |
| 3.3 "z-score 检验输出事件置信度" | `bintest.py::z_prob` L98-113 + `p_adjust_bh` L116-123 | 零假设取中性时 80/369 显著,聚焦真实事件 |
| 3.4 "联合拟合肿瘤纯度 p 与宿主倍性" | `hmm.py::grid_search_purity_ploidy` L575-668;`call.py` clonal 法 | HMM 自估 purity=0.65(真值 0.70) |
| 3.4 "把每段拷贝拆为 major/minor" | `hmm.py::enumerate_states` L36-60;`do_call` 的 cn1/cn2 | 状态空间即 (cn, minor) 二元组 |
| 3.1 "BAF 与深度信号联合分析" | §4.5 三条路径;`do_call` cn1/cn2 | (本 demo 未含 VCF,见 §7 扩展建议) |
| 3.1 "PE/SR 断裂点证据" | **不使用**(无 read-pair 分析) | — (建议与 Delly/Manta 联用) |

---

## 6.5 信号使用对比:CNVkit 在工具图谱中的位置

报告 3.5 的工具表格按数据类型分类;若按"用了哪些信号"重新排列,CNVkit 的定位一目了然:

| 工具 | RD | BAF | PE/SR | 备注 |
|---|---|---|---|---|
| CNVkit(默认) | ✓ | — | — | 纯深度;antitarget 双尺度 bin 是其特色 |
| CNVkit(+VCF) | ✓ | ✓(分割+call) | — | 等位基因特异;仍无断裂点 |
| XHMM / GATK-gCNV / CoNIFER | ✓ | — | — | 纯深度,WES 胚系,PCA/HMM 去噪路线 |
| CONTRA / ExomeDepth | ✓ | — | — | 纯深度,WES,case-control 比值路线 |
| FACETS / Sequenza / ASCAT | ✓ | ✓(全程联合) | — | 深度+BAF 联合拟合 (purity, ploidy, major/minor),CNVkit 的 BAF 通道是它的子集 |
| TITAN | ✓ | ✓ | — | 贝叶斯混合模型,亚克隆分辨率 |
| Delly / Manta / LUMPY | — | — | ✓ | 纯断裂点,SV 工具;CNV 非主要目标 |
| GRIDSS / SvABA | (✓) | — | ✓ | 断裂点 assembly 为主,深度为辅 |

两个推论:

1. **"CNVkit 只用深度"只在默认配置下成立**;但即使加上 BAF,它的联合深度(发射概率级或断点并集)也不及 ASCAT/FACETS 的网格联合拟合彻底——后者的目标函数直接在 (LRR, BAF) 平面上搜索 (purity, ploidy, CN_A, CN_B) 的全局最优。CNVkit 的定位是"稳健、快速、可及的 CN 检测",等位基因特异分析是辅助能力;
2. **与 SV 工具联用是互补而非冗余**(报告 3.1 的 Ho 2020 / Kosugi 2019 论断):断裂点证据补的是 CNVkit 最弱的两项——断点精读级定位、复杂重排解析;深度证据补的是 SV 工具的大片段 CN 检出缺失。临床流程中 `CNVkit(call) ∩ Delly(CNV)` 取交集可显著降假阳性。

---

## 7. 复现与扩展

```bash
# 环境 (一次): conda create -n cnvdemo -c conda-forge python=3.12 numpy pandas scipy \
#              biopython bioframe pysam pyfaidx statsmodels matplotlib reportlab fsspec
#              R 用户另需: Rscript -e 'install.packages("DNAcopy", repos="https://bioconductor.org/packages/3.18/bioc")'
cd /data/16T/luosg/ngs_cnv_report_work/cnvkit_demo
/home/luosg/miniforge3/envs/cnvdemo/bin/python demo_cnvkit_pipeline.py
```

产出:`sim_tumor.cnr`(bin 级)、`sim_tumor.haar.cns`(分割)、`cnvkit_demo_figure.png`(三面板图:原始信号→校正→分割)、`run_output.txt`(全过程数值)。

**扩展方向**(按投入产出排序):

1. **BAF 通道(推荐首选,补齐本文未演示的部分)**:为模拟加 germline het SNPs——每个 bin 2-5 个,B 等位计数按所属区段的 minor 拷贝数二项采样(中性区 p=0.5,cn=1 区 p=1.0,cnLOH 区深度不变但 p=1.0),写成最小 VCF 后:
   - 跑 `segment --method hmm-tumor --variants`,观察状态空间展开后 cnLOH 区段(BAF→1、深度平坦)被独立标出;
   - 跑 `segment --method cbs --variants`,观察 `variants_in_segment` 在 cnLOH 段内插出 BAF 断点(深度方法本身看不到它);
   - 跑 `do_call(variants=..., purity=0.7)`,核对 `cn1/cn2` 与构造时的 minor 拷贝;
2. **过分割/欠分割扫描**:`fdr_q ∈ {1e-2..1e-5}` × `alpha ∈ {1e-2..1e-5}` 网格,统计检出的段数与断点误差,复现报告中"灵敏度-特异度权衡";
3. **antitarget 演示**:加一组 4 kb 的 antitarget bins(深度 15×),体会双尺度分辨率;
4. **对照真实 CNVkit 输出**:把本仓库任一真实样本的 `.cnr/.cns` 替换模拟数据,`tabio.read` 后复用 Step 4-6。

**阅读源码的建议顺序**(按信息密度):

`fix.py::do_fix` → `reference.py::summarize_info` + `descriptives.py`(biweight) → `segmentation/haar.py::haarSeg`(短,250 行核心) → `call.py::_log2_ratio_to_absolute` + `rescale_baf` → `segmentation/hmm.py`(长,但注释极佳) → `cnary.py::CopyNumArray`(数据结构,GenomicArray 之上的 log2/weight/gene 封装) → 想追 BAF 通道再读 `vary.py::VariantArray` 与 `hmm.py::variants_in_segment`。

---

## 附:关键参数速查(`params.py` + 各模块默认值)

| 参数 | 默认值 | 含义 | 位置 |
|---|---|---|---|
| `MIN_REF_COVERAGE` | −5.0 | 参考 log2 低于此 → 丢 bin | params L4 |
| `MAX_REF_SPREAD` | 1.0 | 参考离散度高于此 → 丢 bin | params L5 |
| `GC_MIN/MAX_FRACTION` | 0.3 / 0.7 | GC 越界 → 丢 bin | params L9-10 |
| `INSERT_SIZE` | 250 | edge-effect 的片段长度假设 | params L13 |
| `fdr_q`(haar) | 1e-4 | 断点 FDR q 值 | segmentation/__init__ L131 |
| `alpha`(CBS) | 1e-4 | 置换检验显著性 | 同上 L130 |
| `alpha`(bintest) | 0.005 | bin 级 z 检验阈值(FDR 校正后) | bintest L23 |
| `thresholds`(call) | (−1.1,−0.25,+0.2,+0.7) | threshold 法的 cn 分档 | call L44-49 |
| `DEFAULT_LENGTH_SCALE` | 5e6 bp | HMM p_stay=exp(−d/L) | hmm L26 |
| `TRANSITION_PENALTY` | 8.0 | HMM 大幅切换惩罚 | hmm L27 |
| `BETABINOM_RHO` | 200 | HMM BAF 发射的集中度 | hmm L29 |
| `LOG2_STDEV_FLOOR` | 0.05 | HMM log2 噪声下限 | hmm L28 |
