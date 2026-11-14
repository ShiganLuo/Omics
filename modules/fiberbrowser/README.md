# fiberbrowser — Fiber-seq 单分子染色质图模块(Fig 2H 风格)

把带 Ma/Aq/MM 标签的比对文件渲染为 Peter et al., Cell Rep Methods 2024,
100911 Fig 2H 风格的单分子染色质图。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `sample_sams`(config)/ `--sample-sams`(CLI) | 样本 → 比对文件:SAM 文本或 BAM/CRAM(BAM 需有索引),需带 Ma/Aq/MM 标签。形态 B 字典,路径由 subworkflow 按上游 config outdir 拼接,保证 DAG 依赖 |
| `loci` | 位点文件,每行一个位点:BED 行(chrom<TAB>start<TAB>end[<TAB>name])或 `[NAME=]chrom:start-end`。由 node.py 解析写回 `loci_map`(净化名 → `chrom:start-end`,模块消费,勿手填) |
| `sample_group_dict` | 样本 → 分组标签,渲染为 "s1 (aged)";手填优先,node.py 自动从 meta 分组列补齐 |

解析规则(引擎内部):只取 primary 比对(flag & 0x900 == 0);CIGAR 的
M/=/X 合并为比对段(D/N 为间隙),黑线与覆盖度只按比对段计算;Ma 块
(nuc/msp/fire)按比对段裁剪;MM 标签解出修饰位点,C+m/C+h 与 G-m/G-h
计入 m5CpG(A+a/T-a 不画)。

### 1.2 输出

```
{outdir}/{locus}/{locus}.png            # plain 变体
{outdir}/{locus}/{locus}_zoom.png       # zoom 变体
{outdir}/{locus}/fiberbrowser.sh        # 回执(单条引擎命令,含线程限流导出)
```

### 1.3 参数

config 位于 `Params.fiberbrowser` 下,与引擎 CLI 一一对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `enabled` | — | false | 出图交付总开关;关闭时 node.py 不注册输出,模块不进 DAG |
| `loci` | `--locus` | — | 位点文件(config)/ 单个位点 `chrom:start-end\|name`(CLI);启用时必填 |
| `loci_map` | — | — | node.py 生成:净化名 → `chrom:start-end`;勿手填 |
| `sample_sams` | `--sample-sams` | — | 见输入;CLI 形态为 JSON 字典 `{"s1":"path.sam",...}` |
| `sample_group_dict` | `--sample-group-dict` | {} | 见输入;CLI 形态为 JSON 字典 |
| `variants` | `--variant`(可重复) | plain, zoom | `plain` 出 `<name>.<fmt>`;`zoom` 出 `<name>_zoom.<fmt>`(浅蓝框选 + 灰影重点区 + 虚线 + 150 bp 放大行) |
| `max_fibers` | `--max-fibers` | 12 | 每样本纤维行数上限,按窗口内比对碱基数降序取 |
| `zoom_bp` | `--zoom-bp` | 2100 | zoom 蓝框宽度(bp);灰影重点区为其 1/3 |
| `formats` | `--format`(可重复) | png, svg | 输出格式(png pdf svg ps eps tif tiff jpg jpeg pgf raw rgba) |
| `dpi` | `--dpi` | 200 | 光栅格式分辨率;矢量格式忽略 |

直接调用示例(脱离流程,config 原生 JSON 直传):

```bash
python bin/fiberbrowser.py \
  --locus "NC_133414.1:68551797-68560579|MSTA (ERVL-MaLR)" \
  --sample-sams '{"s1":"region_s1.sam","s2":"s2.fiberseq.fire.sorted.bam"}' \
  --sample-group-dict '{"s1":"aged","s2":"youth"}' \
  --out-dir out --format png --variant plain --variant zoom
```

日志打印每样本 primary 纤维数与 FIRE 调用数(实测 44/188、51/218),
可作数据健全性快检。

## 2. 算法原理与生物学意义

### 2.1 生物学意义

Fiber-seq 长读长测序在单分子分辨率上同时读出核小体定位(足迹)、连接子
DNA、CpG 甲基化,并把甲基化敏感区段(MSP)中的调控元件判别出来。本图把
一个位点上每个样本的 DNA 分子逐条铺开:黑线裸露段 = 单核小体(n),紫块 =
连接子(L),橙/红/棕红 = FIRE 判定的调控元件(启动子 P / 其他 RE),顶部
m5CpG 尖峰给出 CpG 甲基化密度,覆盖度按元素类型分层堆叠 —— 从而直观比较
不同样本(如 aged / youth)在该位点的染色质纤维组织、核小体定位与调控
元件占用的差异。

### 2.2 算法原理

- **核小体 / MSP**:fibertools-rs v0.4.2 `predict-m6a` + `add-nucleosome`
  (最短 75 bp);MSP = 每条纤维上不被核小体足迹占据的区段。
- **m5CpG**:SMRTLink/jasmine 的 5mC 调用(MM 标签 C+m/G-m);本图按
  "胞嘧啶修饰合计"池化(C+h/G-h 含 5hmC),严格 5mC 可只取 m 型代码。
- **FIRE**:semi-supervised 机器学习(Mokapot + XGBoost)把 MSP 分类为
  调控元件或 linker;每个 FIRE 元件带 precision(真调控元件概率,训练
  目标 ≥ 0.90,上限 0.99),q = 1 − precision。
- **FIRE score 轨道公式**:

      S_g = -(50 / C_g) * Σ_i log10(q_i)

  C_g 为覆盖位点 g 的纤维数(真实比对深度),q_i = 1 − precision_i,
  S ∈ [0, 100];覆盖该位点的 FIRE 元件 < 3 个记 −1(画 0),不做平滑。
- **Aq 标签**:`Ma:Z` 中 `fire.Q` 是质量规格而非分数;每调用的
  precision 存于 `Aq` B-array 标签(precision = Aq/255,uint8),与 Ma 串
  的 fire 块严格 1:1 配对;q = 1 − Aq/255。
- **FIRE q 三档着色**(`fire_tier_color`):
  - `aq >= 252` → q ≤ 0.01,棕红 `#7D3219`;
  - `243 <= aq < 252` → 0.01 < q ≤ 0.05,亮红 `#E01B24`;
  - `aq <= 242` → q > 0.05,橙 `#F0A246`。
- **三档着色的生物学意义**:三档对应 MSP 被 FIRE 判为**真实调控元件**的
  统计置信度(q = FDR,即 1 − precision):棕红(q≤0.01)是高置信调控
  元件,亮红(0.01<q≤0.05)置信度次之,橙(0.10≥q>0.05)是边界性提示
  调用,可能部分是假阳性。颜色越深 = 该段无核小体足迹的开放染色质越可能
  承担调控功能(启动子 P / 其他 RE),而非仅仅是连接子 DNA;看图时先盯
  深色块,橙色块作为候选线索。这与图注的 FDR 三档一致(橙 0.10≥q>0.05、
  亮红 0.05≥q>0.01、棕红 q≤0.01)。
- **覆盖度三层实心堆叠**(总高 = 纤维覆盖,先画大层再画内层):
  红 `#C0392B`(FIRE 覆盖)→ 紫 `#8A78B6`(连接子 = MSP − FIRE)→
  灰 `#B0B1B5`(核小体覆盖)。零覆盖缺口 ≥500 bp 标注 `deletion X.X kb`。
- **zoom 变体**:自动选 FIRE 调用最多的代表性纤维,以 q 最小的 FIRE 块
  为中心,蓝框宽 `zoom_bp`、灰影为其 1/3,放大行画该分子真实调用
  (n/L/RE 标注 + 150 bp 比例尺)。

## 3. 参考文献

1. Peter CJ et al. *Single chromatin fiber profiling and nucleosome position
   mapping in the human brain*. Cell Reports Methods 4(12):100911, 2024.
   DOI: 10.1016/j.crmeth.2024.100911 —— 图注与方法原文,Fig 2H / 1D 绘图
   风格来源。
2. FIRE analysis pipeline — https://github.com/fiberseq/FIRE
   (docs/README.md)—— FIRE score 公式、precision/q 约定、
   element_coverages(FIRE / 核小体 / 连接子三种覆盖度)。
3. fibertools-rs(Rust 实现,CLI 缩写 `ft`,文档 docs.rs)v0.4.2 ——
   `predict-m6a` / `add-nucleosome` / `extract`。
4. molecular-annotation crate(docs.rs)—— `Ma:Z` / `Aq` 标签规范
   (`fire.Q` 质量规格、Aq B-array 逐调用质量)。
5. Revio / SMRTLink jasmine —— 5mC(m5CpG)自动调用(方法原文:
   "Methylated CpGs are automatically identified using Revio's
   on-board jasmine")。
