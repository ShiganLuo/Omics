# fibercluster — Fiber-seq 单纤维染色质模式聚类模块

对目标区域（如 marker 基因启动子窗口）提取各样本的单纤维可及性模式并聚类，
输出每类纤维的模式剖面、各样本中的类别比例，以及每区域一张单分子聚类图。
聚类单元 = 1 条纤维 = 1 个细胞的一条等位基因分子；类别比例近似处于该染色质
状态的（等位基因）细胞比例。bulk Fiber-seq 无细胞条码，聚类结果是「染色质
状态亚群」，不是细胞类型归属（见 §3.7）。

本文按"怎么跑 → 怎么读 → 怎么算"组织，各节互不混合：
**§1 输入输出与参数**（接口契约）、**§2 图像结构与元件含义**（纯视觉语义）、
**§3 算法原理**（纯数据推导，不含颜色/绘图描述）。只想看图读 §2；想核对
数值来源读 §3。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `-m/--msp sample=path` | 必填，可重复。ft extract MSP BED12（bgzip）：一行一条纤维，col4 = read name，block = MSP（开放片段）调用。聚类默认通道 |
| `-u/--nuc sample=path` | 可选，可重复。nuc BED12：block = 核小体调用。仅作图上覆盖轨道 |
| `-f/--fire sample=path` | 可选，可重复。`ft fire --extract` 的 fire_peaks.bed 峰表：col5 = FDR。FIRE 通道/覆盖轨道按 `--fire-fdr` 过滤后使用。**不要传 fire.bed.gz**——它承载全部 MSP 区间（已在本项目数据实测核对），显著性只能在峰层过滤 |
| `-r/--regions` / `-c/--locus` | 目标区域：BED（≥4 列，col4 = 区域名）或行内位点 `chr:start-end\|name`（可重复），二者可叠加 |
| `-S/--sample-order` | 样本显示顺序（逗号分隔），首样本在图上方；缺省为 `-m` 传入顺序 |

所有 `sample=path` 的样本名在三个通道间必须一致（同名匹配 read）。BED12
哨兵 block 解析、FDR 过滤等引擎内部规则统一见 §3，本节不重复。


```
{outdir}/                        # 一次 fiber_cluster.py 运行（一个通道/k 设置）
  {region}/
    fiber_matrix.tsv          # 纤维 × bin 二值矩阵（聚类输入的落地）
    cluster_assignments.tsv   # 每条纤维的类别归属
    cluster_patterns.tsv      # 每类的逐 bin 模式（长表）
    cluster_summary.tsv       # 每类 × 每样本的纤维数与比例
    {region}.clusters.png     # 单分子纤维轨道图（结构见 §2）
  regions_failed.txt          # 仅当有区域失败时出现（region<TAB>reason）
```

多次运行（不同通道/k）的结果由 `bin/fiber_report.py` 汇编为报告包
（见 §1.4 的调用示例）：

```
{report_outdir}/
  figures/{region}.clusters.png     # 定稿图：canonical 运行中 ≥2 类别的区域
  summary/loci_summary.tsv          # 全部位点（含单类别位点，只有表不出图）
  summary/cluster_proportions.png   # 跨位点类别比例总览（结论主图）
  sensitivity/sensitivity_table.tsv # 各设置下的类别结构对照
```

表格列含义：

| 文件 | 列 | 含义 |
|---|---|---|
| fiber_matrix.tsv | read_name / sample / cluster / span_frac | 纤维 ID / 样本 / 类别（1 起）/ 纤维跨度覆盖区域的比例（§3.1） |
| | bin_{start} | 各 bin 的二值状态：1 = 开放（判定见 §3.2） |
| cluster_assignments.tsv | n_open_bins | 该纤维开放 bin 总数（Σ bin 列） |
| cluster_patterns.tsv | bin_start / bin_end | bin 的基因组坐标（末 bin 可能短于 bin_size） |
| | open_frac / nuc_frac / fire_frac | 该类内 bin 开放/核小体占位/FIRE 的纤维比例（§3.5） |
| | n_fibers | 该类纤维数 |
| cluster_summary.tsv | n_fibers / frac_of_sample | 该类中该样本的纤维数 / 占该样本全部纤维的比例（§3.5） |

**运行语义**：默认 skip-existing（已存在 `cluster_assignments.tsv` 的区域
跳过），`-w/--overwrite` 重算；`-D/--dry-run` 只打印计划不落盘；区域失败
（如纤维数 < `-n`）收集进 `regions_failed.txt` 并使进程退出码为 1，全部成功
为 0，参数/输入错误为 2。

### 1.3 参数

CLI 即引擎 argparse 全量（含默认值与选择集）：

| 参数 | 默认 | 说明 |
|---|---|---|
| `-r/--regions PATH` | — | 区域 BED；与 `-c` 至少给一个 |
| `-c/--locus STR` | [] | 行内位点 `chr:start-end\|name`，可重复 |
| `-m/--msp S=P` | 必填 | MSP BED12，可重复 |
| `-u/--nuc S=P` | [] | nuc BED12，可重复 |
| `-f/--fire S=P` | [] | fire_peaks.bed，可重复 |
| `-q/--fire-fdr F` | 0.1 | FIRE 峰显著性阈值：保留 FDR ≤ q 的峰 |
| `-o/--outdir DIR` | 必填 | 输出目录 |
| `-b/--bin-size INT` | 100 | bin 宽度（bp） |
| `-B/--bin-mode {any,frac}` | any | bin 开放判定模式（§3.2） |
| `-P/--min-msp-bp INT` | 10 | `--bin-mode any` 时 bin 内最小 MSP 重叠 bp |
| `-k/--clusters INT` | 3 | 目标类别数 |
| `-e/--method {hierarchical,kmeans}` | hierarchical | 聚类算法（§3.4） |
| `-n/--min-fibers INT` | 10 | 区域最少纤维数，低于则跳过并记入失败清单 |
| `-s/--min-span-frac F` | 0.8 | 纤维入选的最小跨度覆盖率（§3.1） |
| `-t/--bin-threshold F` | 0.5 | `--bin-mode frac` 时 bin 开放的覆盖率阈值 |
| `-S/--sample-order STR` | — | 显示顺序，逗号分隔 |
| `-x/--seed INT` | 0 | k-means 随机种子 |
| `-a/--max-fibers INT` | 300 | 单图显示抽样上限（仅显示抽样，统计永远用全部纤维） |
| `-T/--right-tracks STR` | fire,lnk,nuc | 右子图覆盖度轨道的层（逗号分隔，均可选，自下而上 fire→lnk→nuc）：`fire`=FIRE 覆盖、`lnk`=linker(MSP) 覆盖、`nuc`=核小体覆盖；`none` 全部关闭 |
| `-y/--cache-dir DIR` | — | 提取缓存目录：按区域集（md5 校验）缓存元素子集，且在此/源文件旁查找 tabix 索引；重跑秒级 |
| `-w/--overwrite` | 关 | 重算已存在区域 |
| `-D/--dry-run` | 关 | 只打印计划 |
| `-l/--log PATH` | — | 日志文件（缺省仅 stdout） |

### 1.4 示例命令（实际执行过的命令，路径即真实输入）

| 输入文件 | 绝对路径 |
|---|---|
| MSP（青年） | `/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-21310/luanchao-21310.msp.bed.gz` |
| MSP（衰老） | `/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-11236/luanchao-11236.msp.bed.gz` |
| nuc ×2 | 同目录 `luanchao-{21310,11236}.nuc.bed.gz` |
| FIRE 峰 ×2 | 同目录 `luanchao-{21310,11236}.fire_peaks.bed` |
| 区域 BED | `/data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/regions.bed`（marker 基因 TSS±2.5kb，坐标取自 `/data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/GCF_049350105.2_T2T-MMU8v2.0_genomic.gtf`） |
| 输出 | `/data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/{k2,k3,fire_k2,fire_k3}/` |

FIRE 通道 k=2 的单条命令（k=3 换 `-k 3 -o …/fire_k3`；MSP 通道换
`-C msp -o …/k2`）：

```bash
apptainer exec --bind /data/16T:/data/16T /data/16T/Database/env/fiberseq/fiberseq.sif \
  python /data/16T/luosg/workflow/Omics/modules/fibercluster/bin/fiber_cluster.py \
  -r /data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/regions.bed \
  -m luanchao-21310=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-21310/luanchao-21310.msp.bed.gz \
  -m luanchao-11236=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-11236/luanchao-11236.msp.bed.gz \
  -u luanchao-21310=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-21310/luanchao-21310.nuc.bed.gz \
  -u luanchao-11236=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-11236/luanchao-11236.nuc.bed.gz \
  -f luanchao-21310=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-21310/luanchao-21310.fire_peaks.bed \
  -f luanchao-11236=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/common/4_elements/luanchao-11236/luanchao-11236.fire_peaks.bed \
  -b 100 -e hierarchical -S luanchao-21310,luanchao-11236 -C fire -k 2 \
  -y /data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/.extract_cache \
  -o /data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/fire_k2 \
  -l /data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/fire_k2/fiber_cluster.log
```

合成数据冒烟测试（断言：植入模式完全分开、dry-run 无副作用、失败路径、
skip-existing、FIRE FDR 过滤正确）：

```bash
apptainer exec --bind /data/16T:/data/16T --bind /home/luosg/.hermes/cache/scratch:/home/luosg/.hermes/cache/scratch \
  /data/16T/Database/env/fiberseq/fiberseq.sif \
  python /home/luosg/.hermes/cache/scratch/fibercluster_smoke_test.py
```

报告包汇编（定稿图 + 汇总表 + 跨位点总览图 + 敏感性附录，canonical =
fire_k2，即 FIRE 通道 k=2）：

```bash
apptainer exec --bind /data/16T:/data/16T /data/16T/Database/env/fiberseq/fiberseq.sif \
  python /data/16T/luosg/workflow/Omics/modules/fibercluster/bin/fiber_report.py \
  -r fire_k2=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/fire_k2 \
  -r fire_k3=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/fire_k3 \
  -r k2=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/k2 \
  -r k3=/data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot/k3 \
  -c fire_k2 -S luanchao-21310,luanchao-11236 \
  -o /data/16T/luosg/output/luanchao_mulatta/Fiberseq/fiberseq_analysis/fiber_cluster_pilot
```

### 1.5 流程接线（Snakemake）

- **模块文件**：`fibercluster.smk` / `fibercluster.json` / `fibercluster.schema.json`。
  环境复用 `../fiberseq/fiberseq.yaml`（依赖完全一致：numpy/pandas/scipy/
  matplotlib/pysam；按 modules.md「最小环境」原则不另建 yaml/SIF，容器经
  `sif("../fiberseq/fiberseq.yaml")` 解析 env 映射键 `fiberseq`）。
- **规则**:`fiber_cluster`(`{setting}` 通配符,一个通道/k 设置一次引擎调用;
  逐区域产物在 `{outdir}/{setting}/<region>/`,DAG 以标记
  `{outdir}/{setting}/fiber_cluster.done` 跟踪——区域名来自 regions BED,
  按 modules.md「禁用 directory()」用标记文件)。规则按引擎退出码契约放行
  rc=1(= 部分区域失败,清单 `regions_failed.txt`,记 warning 不阻断 DAG);
  仅全区域无产出或 rc∉{0,1} 才判失败(2026-11-22 补丁,否则低覆盖区域
  会卡死 fiber_report)。同批补丁:规则执行前清空本 setting 输出目录(防脏产物
  掩蔽失败检测);回执执行改用 `subprocess.call` 直跑(Snakemake 的 `shell()` 会把
  命令再包一层容器执行且无法捕获退出码),因此命令显式
  `apptainer exec --bind /data/16T <fiberseq.sif> python …`(引擎需 pysam,
  在 fiberseq.sif 内;两条规则同此)。`fiber_report`(汇编报告包,输出 `summary/` 两文件 +
  `sensitivity/sensitivity_table.tsv` + `figures.done`)、`fibercluster_result`
  (收尾聚合,node.py 注册点)。
- **config 键 ↔ 引擎 CLI**:`regions_bed`→`-r`;`samples` + `indir` +
  `Params.fibercluster.fiber_cluster.{msp,nuc,fire}_suffix`→`-m/-u/-f`(形态 A
  拼接 `{indir}/{sample}/{sample}{suffix}`;suffix 置 null 关闭对应通道);
  `settings`(标签列表,即 `{setting}` 通配符取值)每标签一次引擎调用。
  **零配置默认**:`settings` 缺省为 `["auto"]`,此时 `channel` 缺省 msp、
  `k` 缺省 0(引擎 auto-k,见 §3),不需要任何 by_setting/canonical 键。
  敏感性预设仍可用:列出多个标签并为需要偏离全局默认的标签在
  `Params.fibercluster.fiber_cluster.{channel,k}_by_setting` 写覆盖条目
  (无条目的标签回退全局 `channel`/`k`)→`-C/-k`;
  `Params.fibercluster.fiber_cluster.<key>` 与引擎旗标一一对应,仅转发
  config 实际设置的键(engine 默认值单一来源):
  bin_size→-b、bin_mode→-B、min_msp_bp→-P、method→-e、min_fibers→-n、
  min_span_frac→-s、bin_threshold→-t、max_fibers→-a、right_tracks→-T、
  fire_fdr→-q、seed→-x、cache_dir→-y;
  `samples` 列表的顺序即样本显示顺序 → 两条规则的 `-S`(单一来源,无独立
  sample_order 键);`Params.fibercluster.fiber_report.canonical`→fiber_report
  的 `-c`(缺省 = 第一个 setting 标签)。命令回执双写、`OMICS_SCRIPT_ONLY`
  测试模式按 modules.md 规范。

## 2. 图像结构与元件含义

`{region}.clusters.png` 采用 fiberbrowser / 论文 Fig 2H 的单分子纤维画法，
**单面板**、无标题（区域身份在文件名与表格中）：

| 元件 | 呈现 | 含义 |
|---|---|---|
| 分子线（每行黑线） | 沿纤维跨度的黑色细线 | 核小体/中性 DNA（分子主体） |
| 紫色矩形 | linker DNA（非 FIRE 的 MSP） | nuclease 可及的核小体间连接区 |
| 橙/红/深棕红矩形 | FIRE MSP，按显著性三档 | 落在显著 FIRE 峰内且该纤维真实开放的 MSP：橙 0.05<q≤0.10、红 0.01<q≤0.05、深棕红 q≤0.01（RE/P） |
| 行分组 | 不画分隔线；仅靠左侧分色括号与聚类标签区分 | 每段 = 一个聚类；cluster 1 在最上 |
| 子堆叠 | 带内按样本分块 | 一块 = 该类中一个样本的全部（或抽样）纤维 |
| 括号线 + 浅绿圆角标签框（框外左侧） | fiberbrowser 同款：括号线段用该 cluster 的颜色；标签框 `样本 (n=该类内该样本纤维数[, show M])` | 样本身份与该类内的纤维计数；括号颜色 = cluster 颜色 |
| 聚类标签（最左，框外） | `cluster c (n=全类纤维数)`，文字用该 cluster 颜色 | 类别与总数；颜色与括号线、总览图一致（tab10 第 c 档） |
| 行序 | (cluster, 样本, 聚类树叶序) | 相似模式相邻 |
| 右子图（每 cluster 一条覆盖度轨道） | fiberbrowser 覆盖度轨道同款堆叠：红=FIRE 覆盖（底）、紫=linker(MSP) 覆盖（中）、灰=核小体覆盖（顶）；逐碱基计数，纵向恰好占满该 cluster 的行段（与左侧同一行坐标系，严格对齐） | 该类纤维中在该碱基携带该元件的条数（覆盖度）；层可选（`--right-tracks`），堆叠按类内最大总覆盖归一（最大覆盖处顶满行段） |
| 外框 | 无边框，仅保留底部 x 坐标轴 | — |
| x 轴 | 坐标刻度（kb） | 区域实际坐标 |
| 图例（底部，两组分开） | 左子图下方：分子线 / 紫 linker / 三档 FIRE 色；右子图正下方：覆盖度层 | 各自解释本面板元件；右图例随 `--right-tracks` 启用层增减 |

看图要点：某类的结构直接读带内线段排布（FIRE 段多少、位置是否固定、
核小体相位是否一致）；类大小看行数与标签 n；某类几乎单一样本时（标签框
n 一边倒）类别与样本混杂，组间解读需谨慎。纤维数超过 `--max-fibers` 时按
(类, 样本) 比例抽样（标签注 show M），所有计数仍是全量。单类别（无异质性）
位点不出此图，进 `summary/loci_summary.tsv`。

## 3. 算法原理

### 3.1 区域与纤维提取

- 区域合并：所有区域按染色体排序并合并重叠区间，作为重叠判定的查询区间
  （不做坐标改动）。
- 提取访问顺序（快者优先，日志标注实际走的路径）：① **tabix 索引区域跳读**
  ——常规路径；`<file>.tbi` 索引由产出规则 `ft_extract`/`ft_extract_hap`
  随元素 BED 生成（也可用 `tabix -p bed` 预建在源文件旁或 `--cache-dir` 内），
  任意区域集都是秒级；② 区域缓存 `<cache>/<文件名>.subset.bed.gz`（按区域集
  md5 校验有效性）；③ 全文件流式扫描**仅作无索引、无缓存时的兜底**，并回写
  缓存。提取开销与区域数无关，也**不走全文件流式**（流式只是兜底）。
- 记录解析：`blockSizes/blockStarts` 解析为绝对区间列表，保留 span 与任一
  查询区间相交的记录。
- **哨兵 block**：fibertools-rs v0.13.1 的 BED12 每条记录首、末两个 block 是
  边界哨兵（首块 size 0、末块 size 1），不是真实调用；一律丢弃，只保留内部
  block（与 `modules/fibertools/bin/distributed_plot.py` 参考实现一致）。
- 跨度过滤：`span_frac = (min(f_end, end) − max(f_start, start)) / (end − start)`，
  保留 `span_frac ≥ --min-span-frac` 的纤维。未通过的纤维不参与该区域任何
  统计（不是补零）。

### 3.2 bin 二值化

- bin 网格：`n_bins = ceil((end − start) / bin_size)`；bin i 覆盖
  `[start + i·bin_size, min(end, start + (i+1)·bin_size))`，末 bin 允许不足额。
- 每 bin 覆盖 bp = 该纤维所有（哨兵剔除后的）block 与 bin 区间交集长度之和
  （MSP/核小体 block 在记录内互不重叠）。
- 判定模式：
  - `any`（默认）：`覆盖 bp ≥ --min-msp-bp` 记开放。对 100bp bin 即「bin 内
    至少 10bp MSP」，灵敏、保留小 MSP 信号；
  - `frac`：`覆盖 bp / bin 实际宽度 ≥ --bin-threshold` 记开放。更保守，要求
    开放片段铺满 bin 的一半以上。
- 纤维跨度未覆盖的 bin 覆盖 bp = 0 → 记关闭（插补）。对结论的影响由
  `--min-span-frac` 兜底：入选纤维至少铺满区域的 80%（默认）。

### 3.3 聚类通道与 FIRE 语义

- `--cluster-channel msp`：二值化对象 = 全部 MSP block（背景可及性 + 调控
  元件一起进入聚类）。
- `--cluster-channel fire`：二值化对象 = MSP block ∩ 显著 FIRE 峰。峰来源
  `fire_peaks.bed`，保留 `FDR ≤ --fire-fdr` 的峰并按染色体合并；每纤维的
  FIRE block = 该纤维 MSP block 与合并峰集的区间交集（逐对推进先结束方）。
  只有纤维自身带开放 MSP 的峰内区间才算该纤维的 FIRE 信号——FIRE 既要求
  元件显著，又要求这条纤维在该处真实开放。
- FIRE 显著性只能在峰层过滤的原因：`fire.bed.gz`（`ft fire --extract`
  BED9+）承载的是全部 MSP 区间（本项目数据实测与 msp block 逐一相同），
  不含可过滤的逐条 q 值；`fire_peaks.bed` 的 col5 FDR 才是显著性轴
  （如 FDR≤0.1 ⇔ score≥15.32，查 fire_fdr_table.tsv 的阈值→FDR 映射）。
- 覆盖轨道的 nuc 通道固定用 `frac` 模式判定（核小体 block 约 169bp，铺满
  bin 即占位）；fire 轨道与聚类通道同规则。`--cluster-channel fire` 时
  fire 轨道与 open 轨道一致（同一组 block）。

### 3.4 聚类

- 距离：二值向量间的 Hamming 距离 = 两纤维状态不同的 bin 占比
  （`scipy.spatial.distance.pdist(metric="hamming")`）。
- `hierarchical`（默认）：average linkage（UPGMA）合并，
  `scipy.cluster.hierarchy.fcluster(criterion="maxclust", t=k)` 切为 ≤k 类。
- `kmeans`：`scipy.cluster.vq.kmeans2(minit="++", iter=50)`，随机源
  `numpy.random.seed(--seed)`。
- **k 的指定与自动选择**：`-k <int>` 固定 k;`-k 0`（默认，config 键
  `Params.fibercluster.fiber_cluster.k = 0`）为 **auto-k**——对每个区域
  独立在 k=2..6 中选 k，依据同一 Hamming 距离下的 silhouette：
  1. 候选 k 必须切出 k 个非空类，且每类纤维数 ≥ max(2, 5%·n)
     （1–2 条纤维的「类」是噪声，不是染色质状态）；
  2. 全部候选被 guard 排除 = 纤维模式同质 → k=1（单态即真相）；
  3. silhouette 最高者胜；最高分 0.02 以内视为并列,取更小的 k
     （Occam：状态数就简）。
  每个区域落盘 `k_selection.tsv`（k, silhouette 两列）供审阅,日志同步
  打印。合成数据实测：真 2 态→k=2（silhouette 1.000）、真 3 态→k=3、
  同质→k=1。
- 类别重标号：为保证跨运行标签稳定，聚类完成后按「类内平均开放度降序」
  重排标签——**cluster 1 永远是平均最开放的类**，其余依次。
- 退化情形：当二值模式完全相同/大量并列时（如 FIRE 通道下全部纤维为全零
  向量），`maxclust` 可能切不出 k 个非空类；此时输出实际类别数并在日志
  warning，不再报错。纤维总数 < k 仍为错误（该区域记入失败清单）。

### 3.5 派生统计量

- `cluster_patterns.tsv` 的 `open_frac` = 该类所有纤维在该 bin 二值的均值
  （= 该类中该 bin 开放的纤维比例）；`nuc_frac`、`fire_frac` 同定义，分别
  基于核小体/FIRE 二值矩阵。
- `cluster_summary.tsv` 的 `frac_of_sample` = `n(类 ∩ 样本) / n(样本)`；
  同一样本各类之和为 1。`n_fibers` 为原始计数。
- `cluster_assignments.tsv` 的 `n_open_bins` = 纤维二值向量之和；
  `span_frac` 见 §3.1。

### 3.6 图形布局规则（数据侧）

- 显示行序 = (cluster 升序, 样本在 `--sample-order` 中的次序, 聚类树叶序)：
  层次聚类取 `scipy.cluster.hierarchy.leaves_list` 的全局叶序、k-means 取到
  质心距离升序；**不按开放度排序**（会产生人工渐变楔形）。带内按样本分为
  子堆叠，逐一按 fiberbrowser 样式标注。
- 显示抽样：总行数 > `--max-fibers` 时按 (类, 样本) 子堆叠比例等距抽样
  （仅影响显示；计数、比例、表格统计一律用全量纤维）。
- 颜色与元件画法逐项复用 `modules/fiberbrowser/bin/fiberbrowser.py` 的
  调色板与几何（分子线 lw=0.9 黑、linker 矩形高 0.50 紫、FIRE 矩形高
  0.56 按 q 档着色；FDR→档位映射 `fire_tier()`）。
- 所有文字（子堆叠标签框、聚类标签、坐标、图例）均在框外；聚类标签做
  最小间距 0.075 的错位并钳入 [0.03, 0.97]。
- 面板无边框（仅保留底部 x 坐标轴）；括号线段与聚类标签按 cluster 染色
  （tab10 第 c 档，与 `summary/cluster_proportions.png` 的 cluster 配色
  一一对应，两图可互相核对）。
- 右子图覆盖度轨道 = 该类**全量**纤维的逐碱基覆盖计数（该碱基上携带该类
  元件的纤维条数）；三层按 fiberbrowser 覆盖度轨道同款堆叠与配色（自下而上
  fire 红 / linker 紫 / 核小体灰，先画最宽累计层），linker = MSP 去除 FIRE
  段后的剩余区间；纵向按类内最大总覆盖归一（最大覆盖处顶满行段）、恰好占满该类行段（行坐标
  系共享，严格 y 对齐）；层选择 `--right-tracks`（fire/lnk/nuc 任意组合或
  none），图例随启用的层自动增减。

### 3.7 统计语义与限制

- 一条纤维 = 一个细胞的一条等位基因分子。类别是「染色质状态亚群」：可能
  来自不同细胞类型、同一细胞类型的不同状态、或单倍型差异；bulk 数据无细胞
  条码，**无法**把纤维归属到具体细胞，类别≠细胞类型。
- 每组 n=1 生物学重复时，类别比例的组间差异只作描述性比较（region-level），
  不构成样本层显著性；与技术批次的区分需要重复或技术协变量。
- m6A/MSP 是 Hia5 甲基转移酶在开放染色质上沉积的可及性标记，不是天然
  DNA 修饰；`open_frac` 一类指标一律解读为可及性。
- 二值化丢失 MSP 长度与强度信息；`any`/`frac` 两种模式的稀疏度不同，跨
  运行比较必须固定 `--bin-mode`、`--min-msp-bp`/`--bin-threshold`、`--bin-size`。
- 跨度过滤（§3.1）与 bin 插补（§3.2）共同决定了参与统计的纤维集合；低覆盖
  区域（n < `--min-fibers`）不产出结果。

## 4. 参考文献

- Peter CJ, Agarwal A, ..., Stergachis AB, Hasson D, Akbarian S. Single
  chromatin fiber profiling and nucleosome position mapping in the human
  brain. *Cell Reports Methods* 4(12):100911 (2024).
  DOI: [10.1016/j.crmeth.2024.100911](https://doi.org/10.1016/j.crmeth.2024.100911)
  （PMC11704683；Fiber-seq 脑组织适配与单纤维异质性框架、FIRE 流程）
- fibertools-rs v0.13.1（`ft extract` BED12 元件层、`ft fire` FIRE 峰表与
  FDR 表）：[github.com/songlab-cal/fibertools-rs](https://github.com/songlab-cal/fibertools-rs)；
  BED12 哨兵 block 约定与 fire 输出语义已在本项目数据上逐条实测核对
- SciPy 层次聚类与 k-means：
  [docs.scipy.org/doc/scipy/reference/cluster.html](https://docs.scipy.org/doc/scipy/reference/cluster.html)
