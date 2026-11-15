# fiberbrowser — Fiber-seq 单分子染色质图模块(Fig 2H 风格)

把带 Ma/Aq/MM 标签的比对文件渲染为 Peter et al., Cell Rep Methods 2024,
100911 Fig 2H 风格的单分子染色质图。

本文按"怎么跑 → 怎么读 → 怎么算"组织,三节各司其职、互不混合:
**§1 输入输出与参数**(接口契约)、**§2 图像结构与元件含义**(纯视觉语义:
图上每个元件长什么样、表示什么)、**§3 算法原理**(纯数据推导:每个量怎么
算出来,不含颜色/绘图描述)。只想读图看 §2;想核对数值来源看 §3。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `sample_sams`(config)/ `--sample-sams`(CLI) | 样本 → 比对文件:SAM 文本或 BAM/CRAM(BAM 需有索引),需带 Ma/Aq/MM 标签。形态 B 字典,路径由 subworkflow 按上游 config outdir 拼接,保证 DAG 依赖 |
| `loci` | 位点文件,每行一个位点:BED 行(chrom<TAB>start<TAB>end[<TAB>name])或 `[NAME=]chrom:start-end`。由 node.py 解析写回 `loci_map`(净化名 → `chrom:start-end`,模块消费,勿手填) |
| `sample_group_dict` | 样本 → 分组标签,渲染为 "s1 (aged)";手填优先,node.py 自动从 meta 分组列补齐 |

(比对解析的引擎内部规则——primary 过滤、CIGAR 分段、缺口桥接、Ma/MM
标签消费——统一见 §3.1,不在本节重复。)

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
| `variants` | `--variant`(可重复) | plain, zoom | `plain` 出 `<name>.<fmt>`;`zoom` 出 `<name>_zoom.<fmt>`(框选 + 放大行,元件见 §2.5) |
| `max_fibers` | `--max-fibers` | 12 | 每样本纤维行数上限,按窗口内比对碱基数降序取 |
| `zoom_bp` | `--zoom-bp` | 2100 | zoom 蓝框宽度(bp);灰影重点区为其 1/3 |
| `gtf` | `--gtf` | — | 可选 GTF/GTF.gz 注释;设置后每张图底部加基因结构轨道,zoom 灰影延伸贯穿该轨道(元件见 §2.4) |
| `zoom_center` | `--zoom-center` | — | 蓝框(与灰影)中心的基因组坐标;显式指定时优先于一切自动规则。缺省自动选择(规则见 §3.7);越界值收敛到位点内,落在无读段覆盖处会告警 |
| `zoom_center_mode` | `--zoom-center-mode` | — | `zoom_center` 为空时的自动选心规则:`diff-fire` = 两组间 FIRE score 差值最大处;`diff-cpg` = 两组间 m5CpG 差值最大处(恰好需要两个分组;算法见 §3.7) |
| `sample_order` | `--sample-order` | — | 逗号分隔 sample_id,固定样本自上而下显示顺序(如 `s_youth,s_aged`);未列出的样本按输入序排在所列样本之后 |
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

## 2. 图像结构与元件含义(怎么读图)

本节只回答"图上画的是什么":每个元件的视觉呈现与含义、看图时怎么解读。
所有量的计算方式(数值从哪来)一律在 §3,本节不重复。全图无坐标轴
(UCSC 截图风格),横轴即位点线性坐标。

### 2.1 整体布局(自上而下)

```
标题:<name>  (chrom:start-end)
图例:一行 4 个色块(L / 橙 / 亮红 / 棕红,见 §2.2)
─────────────────────────────────────────────
m5CpG 轨道      × N 样本(每样本一行,堆叠在最顶部)
─────────────────────────────────────────────
[样本 1]  FIRE score 轨道      ┐
          coverage 轨道        ├ 每样本一组,三行固定顺序
          Single Fibers 纤维行 ┘(左侧绿色标签盒 + 方括号)
[样本 2]  同上三行
...
─────────────────────────────────────────────
gene 轨道(可选,--gtf 时出现,每个基因模型一行)
zoom 放大行(可选,--variant zoom 时出现,见 §2.5)
```

### 2.2 单纤维行的元件

每行 = 一条 DNA 分子(单个 read),从左到右按位点坐标铺开:

| 元件 | 视觉呈现 | 含义 |
|---|---|---|
| 核小体(n) | 黑色横线 `#1A1A1A` 的可见段 | 单核小体足迹——黑色线本身是纤维骨架,两个元件块之间露出的线段即核小体 |
| 连接子(L) | 紫色矩形 `#8A78B6` | 连接子 DNA = MSP 中未被判为调控元件的部分 |
| 调控元件(RE/P) | 橙 `#F0A246` / 亮红 `#E01B24` / 棕红 `#7D3219` 矩形 | FIRE 判定的调控元件(P = 启动子 / RE = 其他),按置信度三档着色(见下表) |
| 断口 | 黑线断开处 | >150 bp 的真实缺失/大缺口(小缺口已被桥接,见 §3.1) |

FIRE 三档着色与置信度(q = FDR = 1 − precision):

| 颜色 | q 区间 | Aq 值(u8) | 含义 |
|---|---|---|---|
| 棕红 `#7D3219` | q ≤ 0.01 | ≥ 252 | 高置信调控元件,基本可信 |
| 亮红 `#E01B24` | 0.01 < q ≤ 0.05 | 243–251 | 置信度次之 |
| 橙 `#F0A246` | q > 0.05 | ≤ 242 | 边界性提示调用,可能部分是假阳性 |

(与论文图注的 FDR 三档一致:橙 0.10≥q>0.05、亮红 0.05≥q>0.01、
棕红 q≤0.01。)

### 2.3 三条数据轨道

每样本固定三条轨道,自上而下:

| 轨道 | 视觉呈现 | 读数含义 |
|---|---|---|
| m5CpG | 棕红 `#7D3219` 细尖峰(离散直方,非填充轮廓) | 每根尖峰 = 该位置被修饰的胞嘧啶计数,峰高 = 覆盖该位点的纤维中 m5CpG 个数——即 CpG 甲基化密度 |
| FIRE score | 黑色 `#1A1A1A` 实心块,分段常数 | 每个位置的 FIRE 评分 S(公式见 §3.5),块越高 = 该区域调控元件证据越强;<3 个 FIRE 覆盖的位置记 0 |
| coverage | 三层实心堆叠:底红 `#C0392B` → 中紫 `#8A78B6` → 顶灰 `#B0B1B5` | 总高 = 该位置纤维覆盖度(测序深度);红层 = FIRE 元件覆盖,紫层 = 连接子覆盖,灰层 = 核小体覆盖(三层定义见 §3.6) |

coverage 轨道内的灰色小字 `deletion X.X kb`:内部零覆盖缺口 ≥500 bp 时
标注,提示该区间可能是纯合缺失或组装缺口。

### 2.4 标注与辅助元素

| 元件 | 视觉呈现 | 含义 |
|---|---|---|
| 标题 | 顶部居中 | 位点名 + 坐标 `chrom:start-end` |
| 图例 | 标题下方一行 4 色块 | L 连接子 / 三档 FIRE(标签文字与 §2.2 表一致) |
| 样本标签盒 | 纤维栈左侧浅绿盒 `#D8ECD8` + 黑色方括号 | 样本名(含分组,如 "s1 (aged)")+ "Single Fibers",括号圈出该样本的纤维栈范围 |
| gene 轨道 | 深蓝 `#2F4F6F` 外显子盒 + 黑色骨架线 + 线上箭头 | 基因结构:盒 = 外显子,线 = 内含子骨架,箭头 = 转录方向(>/< = 正/负链),基因名标于模型上方;每个基因模型一行 |

### 2.5 zoom 变体特有元件

zoom 变体 = plain 全图 + 一处选区 + 一行放大,用于看清单分子细节:

| 元件 | 视觉呈现 | 含义 |
|---|---|---|
| 蓝框 | 浅蓝 `#C9E4F2` 矩形(边 `#7FB3D0`) | 框住代表纤维上被放大的区间,宽 `zoom_bp`;整个图只有这一处选区 |
| 灰影重点区 | 半透明灰 `#B9B9B9` 竖带,连续贯穿全部数据轨道 | 蓝框中央 1/3 = 放大行实际放大的范围;延伸贯穿 gene 轨道,便于与基因结构对位 |
| 灰虚线 | 两条 `#9A9A9A` 虚线 | 连接灰影底边两端与放大行内灰影边缘——引导视线,把主图与放大行对上 |
| 放大行 | 独立一行,代表纤维的真实调用放大后重画 | 分子细节:黑线(核小体)/ 紫块(L)/ 三档色块(RE)同 §2.2;行内灰影延续,提示当前放大范围 |
| n / L / RE 标签 | 放大行元件上方黑色文字 | 三类结构各标一次(标在该类最宽的实例上;标签冲突时垂直错开)——n = 核小体、L = 连接子、RE = 调控元件 |
| 150 bp 比例尺 | 放大行右侧黑线 + "150 bp" 文字 | 放大行的长度基准 |

### 2.6 看图要点

- **颜色越深 = 越可能是真调控元件**:棕红块是高置信开放染色质区(RE/P),
  看图先盯棕红,亮红次之,橙色块只作候选线索(假阳性率较高)。
- **样本间比较**:同一位点上下对齐比较各组的纤维行组织(核小体定位
  滑动/丢失)、FIRE 块数量与档位、coverage 分层比例、m5CpG 峰密度——
  这正是 Fig 2H 要展示的 aged / youth 染色质纤维组织差异。
- **FIRE score 轨道是区域级汇总**,适合先用它锁定差异区段,再下钻到
  纤维行看单分子证据;zoom 变体自动(或按 `zoom_center`)放大到最值得
  细看的区段。

## 3. 算法原理

本节只讲每个量是怎么从数据推导出来的,不描述任何绘图/颜色;视觉语义
一律见 §2。

### 3.1 比对解析(引擎内部)

- 只取 primary 比对(flag & 0x900 == 0)。
- CIGAR 的 M/=/X 合并为比对段,D/N 为间隙;覆盖度只按比对段计算。
- 画纤维黑线时 ≤150 bp 的比对间隙自动桥接:分子是连续 DNA,zoom 下
  ~40 bp 的比对缺口会在元件链上显示为断裂;更大间隙 = 真实缺失,保留
  断口。
- Ma 块(nuc/msp/fire)按比对段裁剪。
- MM 标签解出修饰位点:C+m/C+h 与 G-m/G-h 计入 m5CpG(A+a/T-a 不画)。

### 3.2 核小体与 MSP

fibertools-rs v0.4.2 `predict-m6a` + `add-nucleosome`(最短 75 bp);
MSP = 每条纤维上不被核小体足迹占据的区段。渲染时 MSP 按 FIRE 调用切分:
落在 FIRE 调用内的子区间 = 调控元件,剩余子区间 = 连接子(§2.2 的
紫块与三档色块由此而来)。

### 3.3 m5CpG

SMRTLink/jasmine 的 5mC 调用(MM 标签 C+m/G-m);本图按
"胞嘧啶修饰合计"池化(C+h/G-h 含 5hmC),严格 5mC 可只取 m 型代码。

### 3.4 FIRE 调用与 q 值

- **FIRE**:semi-supervised 机器学习(Mokapot + XGBoost)把 MSP 分类为
  调控元件或 linker;每个 FIRE 元件带 precision(真调控元件概率,训练
  目标 ≥ 0.90,上限 0.99),q = 1 − precision(q 即 FDR)。
- **Aq 标签**:`Ma:Z` 中 `fire.Q` 是质量规格而非分数;每调用的
  precision 存于 `Aq` B-array 标签(precision = Aq/255,uint8),与 Ma 串
  的 fire 块严格 1:1 配对;q = 1 − Aq/255。precision 上限 0.99 在 u8
  量化下取 252,故 Aq ≥ 252 即 q ≤ 0.01 档(三档 → 颜色的映射见 §2.2
  表,此处不重复)。

### 3.5 FIRE score 轨道公式

      S_g = -(50 / C_g) * Σ_i log10(q_i)

C_g 为覆盖位点 g 的纤维数(真实比对深度),q_i = 1 − precision_i,
S ∈ [0, 100];覆盖该位点的 FIRE 元件 < 3 个记 −1(画 0),不做平滑
(真实 S 为分段常数)。

### 3.6 覆盖度轨道

三类逐碱基覆盖分别累计后堆叠:FIRE 覆盖、连接子覆盖(= MSP − FIRE)、
核小体覆盖;总高 = 三者之和 = 纤维覆盖度。零覆盖缺口 ≥500 bp 触发
deletion 标注(阈值在此,视觉见 §2.3)。

### 3.7 纤维选取与 zoom 选区

- **纤维行选取**:每样本按窗口内比对碱基数降序取前 `max_fibers` 条。
- **代表纤维**:窗口内 FIRE 调用数最多的纤维(平手取比对碱基多者)。
- **放大中心**,按优先级:
  1. `zoom_center` 显式坐标(收敛到位点内;落在任何样本都无覆盖处时
     告警,放大行可能为空);
  2. `zoom_center_mode = diff-fire` / `diff-cpg`:恰好两个分组时,取
     两组间 FIRE score / m5CpG 逐碱基差值绝对值最大的位置(单侧无
     覆盖的位置不参与比较,差值按 200 bp 移动平均平滑后取极值);
  3. 缺省:代表纤维上 q 最小的 FIRE 块中点;该位点无 FIRE 则取位点
     中点。
- **选区尺寸**:蓝框宽 = min(`zoom_bp`, 位点宽度);灰影重点区 = 蓝框
  中央 1/3(下限 50 bp)。放大行重画代表纤维在该区间内的真实调用,
  n/L/RE 标签各标一次(取该类最宽实例)。

### 3.8 生物学意义

Fiber-seq 长读长测序在单分子分辨率上同时读出核小体定位(足迹)、连接子
DNA、CpG 甲基化,并把甲基化敏感区段(MSP)中的调控元件判别出来。本图把
一个位点上每个样本的 DNA 分子逐条铺开,从而直观比较不同样本(如
aged / youth)在该位点的染色质纤维组织、核小体定位与调控元件占用的差异。

## 4. 参考文献

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