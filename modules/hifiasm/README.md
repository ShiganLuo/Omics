# hifiasm — HiFi 基因组组装模块

对 PacBio HiFi BAM 从头组装,产出主 contig FASTA;`hifiasm_qc` 规则对组装
结果做质量评估(统计 + 端粒端捕获)并出图。独立成模块供多个下游消费
(centromere 着丝粒分析、telomere `assembly_scan`/`tidk` 等)。输出布局契约:
`{outdir}/{sample_id}/assembly/asm.bp.p_ctg.fa`,下游按 `assembly_dir` 配置键
或共享 outdir 常量接线。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `inputs.hifiasm.bam` / `inputs.hifiasm.pbi` | dict-mode(优先):样本 ID → 完整 BAM / `.bam.pbi` 路径字典(Fiberseq 的 `{raw}/{sample}/{sample}.bam` staging 布局) |
| `indir` + `bam_substring` | 形态 A(缺省):`{indir}/{sample_id}[.{bam_substring}].bam` + 同名 `.bam.pbi` 索引(PacBio HiFi BAM) |
| `samples` | 样本列表 |

`hifiasm_qc` 的输入是同模块 `hifiasm_assemble` 的输出 FASTA(布局契约路径),
DAG 自动依赖,无需额外配置。

### 1.2 输出

```
{outdir}/{sample_id}/assembly/asm.bp.p_ctg.fa        # 组装主 contig(下游契约路径)
{outdir}/{sample_id}/assembly/asm.bp.p_ctg.gfa       # hifiasm 主 contig 图(工具产物)
{outdir}/{sample_id}/assembly/asm.hifiasm.stamp      # GFA 复用指纹(BAM 路径+大小+mtime)

{outdir}/{sample_id}/assembly/qc/{sample_id}_assembly_stats.tsv   # 每条 contig 组成表
{outdir}/{sample_id}/assembly/qc/{sample_id}_assembly_stats.txt   # 组装统计摘要
{outdir}/{sample_id}/assembly/qc/{sample_id}_telomere_ends.tsv    # 每条 contig 两端端粒调用
{outdir}/{sample_id}/assembly/qc/{sample_id}_telomere_qc.txt      # 端粒捕获质量摘要
{outdir}/{sample_id}/assembly/qc/{sample_id}_assembly_overview.<fmt>  # 组装总览图
{outdir}/{sample_id}/assembly/qc/{sample_id}_telomere_qc.<fmt>       # 端粒 QC 图
{outdir}/{sample_id}/assembly/qc/{sample_id}_telomere_map.<fmt>      # 端粒捕获图谱
```

日志:`{logdir}/{sample_id}/hifiasm_assemble.log`、`{logdir}/{sample_id}/hifiasm_qc.log`。
工具自身输出(hifiasm `[M::...]` 阶段统计、samtools 提示)直接追加在**同一规则
日志**里;注意 hifiasm 在读段加载/纠错期几乎不输出,统计集中在运行结束时成批
打印。引擎另在日志里补时间线:工具版本、Stage 1/2(流式组装)与 Stage 2/2
(GFA→FASTA)的边界与耗时、contig 数/总长/N50,便于中途判断进度。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `logs` | 目录 |
| `bam_substring` | 文件名插入段 | `""` | 见输入(形态 A) |
| `inputs.hifiasm.bam` / `.pbi` | — | 无 | 见输入(dict-mode,优先于形态 A) |
| `Procedure.hifiasm` / `samtools` / `python` | — | 同名命令 | 可执行路径 |
| `--force-reassemble` | 引擎旗标 | 关 | 忽略已有 GFA 与 hifiasm 断点缓存,强制重新组装;缺省时 GFA 带指纹 stamp 可跨重跑复用(阶段 2 失败/中断后重跑不重装) |
| 线程/内存 | — | 32 / `mem_mb=64000` | 真实入口是 `resources.hifiasm_assemble.threads/.mem_mb`(`rule_res()` 解析并按宿主预算钳制;hifiasm 线程扩展到 ~32 收益明显,再往上递减) |
| `Params.hifiasm.qc.enabled` | — | `true` | `hifiasm_qc` 开关(node.py 按此注册输出;false 时整个 QC 不进 DAG) |
| `Params.hifiasm.qc.formats` | `-F` (可重复) | `["png"]` | 出图格式列表,多格式同时输出(合法值见 modules.md §脚本撰写规范 #9) |
| `Params.hifiasm.qc.scan_length` | `--scan_length` | `100000` | 每条 contig 两端扫描窗口(bp) |
| `Params.hifiasm.qc.min_tel_len` | `--min_tel_len` | `1000` | 判定端粒捕获的最小端粒重复长度(bp) |
| `Params.hifiasm.qc.end_margin` | `--end_margin` | `1000` | 端粒簇距 contig 末端的最大间隙(bp),超过判为内部端粒序列(ITS) |
| `Params.hifiasm.qc.top_contigs` | `--top-contigs` | `30` | 端粒捕获图谱绘制的 contig 条数 |
| `Params.hifiasm.qc.expected_chrom_ends` | `--expected-chrom-ends` | `0`(关闭) | 完整单倍体组装应含的染色体端数;node.py 按 `genome.default` 自动填(= 二倍体染色体数,如恒河猴 42、人 46),显式配置优先 |

环境:`conda: "hifiasm.yaml"`(hifiasm、samtools、python、numpy、matplotlib)。
**SIF 复用**:容器模式下显式映射 `"hifiasm": "<env_dir>/centromere/centromere.sif"`
(探针验证该镜像完整覆盖本模块 yaml 依赖,含 numpy/matplotlib 导入冒烟,
2026-11 重验通过;centromere.sif 由 centromere.yaml 治理,升级或重建后须对
本模块重跑探针:import 覆盖 + `--help` 冒烟)。

直接调用示例(模块引擎命令,与回执脚本一致):

```bash
# samtools fastq 直接管道进 hifiasm,无中间 FASTQ;GFA→FASTA 由引擎完成
python modules/hifiasm/bin/hifiasm_assemble.py \
  --bam s1.bam --out-prefix out/s1/assembly/asm \
  --out-fasta out/s1/assembly/asm.bp.p_ctg.fa --threads 32

# 组装质量评估:统计 + 端粒端捕获,单次流式扫描 FASTA
python modules/hifiasm/bin/hifiasm_qc.py \
  --fasta out/s1/assembly/asm.bp.p_ctg.fa --sample_name s1 \
  --output_dir out/s1/assembly/qc --expected-chrom-ends 42
```

## 2. 结果解读

### 2.1 表格输出

`*_assembly_stats.tsv`:每条 contig 一行 —— `contig`、`length`、`gc_bases`/
`gc_frac`(按非 N 碱基计)、`n_bases`/`n_frac`,按长度降序。

`*_assembly_stats.txt`:组装级摘要 —— contig 数、总长、最长/均值/中位数、
N50/L50、N90/L90、auN、GC、N 占比,以及 ≥1kb/10kb/100kb/1Mb 各档 contig
数与合计长度。

`*_telomere_ends.tsv`:每条 contig 一行 —— 两端各报 `tel_5p_length`/
`tel_3p_length`(端粒重复簇跨度)、`tel_*_end_gap`(簇距 contig 末端的间隙)、
`tel_*_motif`(主导链向 CCCTAA/TTAGGG)、`tel_*_supporting_windows`(簇内
达到端粒密度阈值的窗口数)与 `status`(`both`/`5prime_only`/`3prime_only`/
`none`,见 §3.2 判定)。

`*_telomere_qc.txt`:端粒捕获质量摘要 —— ≥1Mb contig 的末端状态分布、
全部 contig 的 5'/3' 捕获数与总数、相对 `expected_chrom_ends` 的完整度
百分比、捕获端粒长度分布(均值/中位/极值)。零捕获时给出可能原因清单。

### 2.2 图像结构与元件含义

**`*_assembly_overview.<fmt>` — 组装总览(2×2 面板,自上而下、从左到右):**

1. **Contig 长度分布** —— log10 长度直方图;绿色虚线 = N50,红色虚线 =
   中位数(图例标注数值)。看图要点:双峰(长染色体臂 contig + 短碎 contig)
   是长读长组装的典型形态,短碎 contig 越少越好。
2. **Nx 曲线** —— X 轴累计组装跨度(%),Y 轴 Nx 长度(log);绿点标注
   N50/N90。曲线平台越长、尾部下坠越晚,连续性越好。
3. **GC vs contig 长度** —— log10 长度(X)× GC%(Y)散点。长 contig 的点
   应聚在窄 GC 带;散布异常宽或出现孤立 GC 群提示污染/嵌合。
4. **Rank-size 曲线** —— contig 按长度排名(log X)× 累计跨度(Y)。
   前几名 contig 占比越高、越早逼近 100%,组装越连续。

**`*_telomere_qc.<fmt>` — 端粒区域组装质量(2×2 面板):**

1. **Contig 末端状态** —— 4 类条形计数(两端都有端粒/仅 5'/仅 3'/无,
   统计对象为 ≥1Mb contig,标题注明 n)。绿色条(两端)即"疑似完整染色体"。
2. **端粒长度分布** —— 捕获端(≥ min_tel_len)的长度直方图,蓝 = 5' 端、
   红 = 3' 端堆叠;绿色虚线 = 中位数。
3. **端粒 landscape** —— contig 按长度排名(X,≥100kb 的前 200 条)×
   端粒长度(Y),蓝圆 = 5' 端、红三角 = 3' 端。点只出现在排名靠前的行 =
   大 contig 才捕获到端粒,符合预期。
4. **端粒完整度** —— 水平堆叠条:绿 = 已捕获端数,灰 = 缺失(相对
   `expected_chrom_ends`,未配置时相对全部 contig 端);标题给完整度百分比。

**`*_telomere_map.<fmt>` — 端粒捕获图谱(按长度取前 N 条 contig):**

- 每条 contig 一根横条(浅蓝 = contig 主体),Y 轴标注 contig 名与长度;
- 两端小色块:绿 = 端粒捕获(旁注长度,如 `12.3kb`),红 = 未捕获(旁注
  `none`);注释一律排在条形外侧,不遮挡数据;
- 看图要点:理想情况下染色体级 contig 两端都是绿块;大片红块 = 组装在
  端粒/亚端粒区断裂。

## 3. 算法原理

### 3.1 组装(hifiasm_assemble)

hifiasm 用串图(BOG)组装 HiFi 读段,纠错后输出单倍型解析
(`bp.hap1/hap2.p_ctg`)与主 contig(`bp.p_ctg`)图;本模块取主 contig 的
GFA `S` 记录转为 FASTA(即 `asm.bp.p_ctg.fa`)。引擎把 `samtools fastq`
的输出流直接接 hifiasm 的 `/dev/stdin`(实测 0.25.0 支持管道输入),
不落中间 FASTQ:省一次全量写读、免单线程 gzip 压缩成本。

### 3.2 组装质量评估(hifiasm_qc)

单次流式扫描 FASTA(不整读内存):每条 contig 累计长度/GC/N 计数,并保留
两端各 `scan_length` 窗口序列。

**统计量**:N50 = 按长度降序累加至总长 50% 处的 contig 长度(L50 为所用
contig 数);N90/L90 同理;auN = Σlength²/Σlength(长度加权连续性指标);
GC 按非 N 碱基计。

**端粒端扫描**:在每端窗口内正则匹配端粒基序(CCCTAA/TTAGGG),相邻命中
间隔 <30bp 聚为一簇;**末端近端优先**——先在距 contig 物理末端 ≤
`end_margin` 的簇中取跨度最大者(真末端端粒),无近端簇时才回退取最大簇
(此时 `end_gap` > end_margin,判为内部端粒序列 ITS,不计入捕获)。窗口
支持数 = 簇内 200bp 滑窗(步进 50bp)中端粒碱基占比 ≥50% 的窗口数。

**捕获判定**:`tel_length >= min_tel_len` 且 `end_gap <= end_margin` 记为
一次捕获;contig 状态按两端捕获组合归为 `both`/`5prime_only`/`3prime_only`/
`none`。**完整度** = 捕获端数 / `expected_chrom_ends`(单倍体组装的染色体
端数 = 2 × 单倍体染色体数 = 二倍体染色体数,如恒河猴 21×2 = 42)。

**判读注意**:捕获数远低于期望端数提示组装未跨越多数端粒;两端捕获长度
完全相同的 contig 需警惕组装伪影;`tel_length > 0` 但 `supporting_windows`
为 0 的信号多为散在重复而非真端粒(已在 `*_telomere_ends.tsv` 中透明呈现,
便于人工复核)。

## 4. 参考文献

1. hifiasm — https://github.com/chhylp123/hifiasm(Cheng H et al.,
   Nature Methods 18:608–610, 2021. DOI: 10.1038/s41592-020-01056-5)
2. 端粒基序 (TTAGGG)n — Moyzis RK et al., A highly conserved repetitive
   DNA sequence, (TTAGGG)n, present at the telomeres of human chromosomes,
   PNAS 85:6622–6626, 1988。
