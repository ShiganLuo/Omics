# 资源平衡机制设计

[← 返回 README](../README.md)

**一句话**：Omics 的资源分配是"主机感知 + 数据量感知 + I/O 感知 + 失败兜底 + 实测校准"的自适应体系，目标是在真实生产负载上最小化整体运行时间（makespan）。

**状态**：全部已实现（2026-11-28 落地，含 I/O 感知部分 §7–§10）；试点评测与参数校准待进行，见行动方案（本机 `~/.hermes/docs/io-aware-scheduling-plan.md`）。

---

## 1. 设计目标（总纲）

**目标：在真实生产负载上最小化整体运行时间（makespan）。**

单一目标函数，所有机制环节（§2）都为它服务。优先级（字典序）：

1. **目标**：makespan 最短；
2. **约束**（必须满足，不是优化对象）：无重试耗尽导致的最终失败、无人工改资源重跑、不改变任何规则的输入输出语义、向后兼容（未超预算的值与静态声明一致；未声明 I/O 资源的规则行为逐比特不变，见 §8.7）、无人为配置每个部署机器（自动探测 + 兜底常量，在任何主机上成立）、无"单 job 需求超预算 → 永不被调度"的假死；
3. **代理指标**（仅当目标难以直接观测时参考）：内存浪费率、核利用率、磁盘 %util、I/O 队列深度（aqu-sz）、await、iowait——它们改善不必然改善 makespan（把并发压到 1 也能让 await 暴跌但 makespan 恶化），冲突时以实测 makespan 为准。

**作用通道声明**（防误读，I/O 维度）：`io_bandwidth_mb_s` 是一个**调度打包维度**——它只决定多少个 job 同时运行，不给任何进程限速（cgroup io.max 之类的进程级限速是另一条通道，明确列为非目标）。对该旋钮的任何效果评估都必须在"打包/并发度"这条通道上做。

**非目标**（明确不做）：多用户公平调度；集群/SLURM 级调度（本机制只管单机 snakemake 进程）；无护栏的自动改写 config（合入只走 `calibrate_resources.py --apply` 的护栏：n≥3、±30% 限幅，见 §6）；追求声明值"精确"本身（只在影响 makespan 的范围内校准）；进程级 I/O 限速（cgroup blkio/io.max）；对未声明规则自动推断带宽需求（I/O 带宽预算只对显式声明的规则生效）；改变 benchmark 采集机制。

**自查条款**（一切魔法数须对目标给出理由或显式标注为待校准近似）：

- clamp 的 0.8 余量、cores 的 90% 份额、份额下限 8 GB、p90×1.3——均为经验近似，由校准闭环（§6）按实测数据重定；
- 多工作流按**输入字节**加权是"工作量大头"的粗代理，不等于关键路径；存在"小数据但串行深的工作流被饿死 → 总时长反升"的理论风险。按待跑 job 数/实测耗时加权（同一 benchmark 数据源）是备选口径，字节加权作为 v1 近似保留并监控。
- I/O 侧：磁盘带宽常量表（§8.3）与 `io_safety_factor = 0.8` 均为**冷启动暂定值**，标注待校准，由 I/O 校准扩展（§10）的实测数据重定；规则权重取 benchmark p90 是"需求上界"的稳健口径，p90 对输入规模分布敏感，锚定配方见 §8.3 参数表；I/O 机制压缩的是"排队等待造成的空转"，不是 I/O 总量——I/O 总量的压缩走临时写入分流（§9）。

---

## 2. 机制总览（一个整体方案）

以下环节共同构成 Omics 的资源平衡方案——它们是同一套闭环的组成部分，不是互相替代的备选项：

| 环节 | 作用 | 背景问题 | 章节 | 状态 |
|---|---|---|---|---|
| 主机探测与全局预算 | cores auto、多工作流预算切分、clamp 防死锁 | 并行度不足、声明超预算假死、多工作流超卖 | §3 | 已实现 |
| 数据量感知的资源缩放 | 资源 = 输入规模的函数，逐规则配 scale | 静态声明与数据量失配 | §4 | 已实现 |
| attempt 升级与自动重试 | 失败自动重试 + 内存翻倍升级 | OOM 一次即失败 | §5 | 已实现 |
| benchmark 校准闭环 | 实测峰值 → 建议表 → 护栏合入，持续校准 | 声明值失真、无法收敛 | §6 | 已实现 |
| I/O 带宽预算 | 规则声明带宽需求，全局上限，调度打包不超预算 | I/O 无预算的并发挤压、多工作流抢盘 | §7–§8 | 已实现 |
| 临时写入分流 | tmp 写入与最终输出分盘 | 同盘写放大、寻道交错 | §9 | 已实现 |
| I/O 校准扩展 | calibrate 增加 io 建议列；可选磁盘带宽探测 | I/O 参数无实测回填通道 | §10 | 已实现 |

分工与闭环：**声明侧**（§4 缩放 + §6 校准）让资源声明贴近真实需求；**调度侧**（§3 主机预算 + §8 I/O 预算）保证并发排得下；**执行侧**（§5）兜住估偏导致的失败；I/O 侧另有降低需求的写入分流（§9）与参数回填（§10）。整个体系形成"声明 → 调度打包 → 执行 → 实测 → 校准 → 更准的声明"的闭环，其中 §4 的缩放参数唯一可靠来源是 §6 的实测数据，§5 兜住 §4 估偏的风险，§10 让 I/O 侧参数随实测收敛。

---

## 3. 主机探测与全局预算（cores auto、多工作流切分、clamp 防死锁）

**背景**：snakemake 的 `--cores` / `--resources` 是静态数字，而部署主机的核数、内存各异——静态声明在不同机器上含义不同，并行度常给不足；多个工作流并发时又会超卖内存；最隐蔽的是单 job 声明超过全局预算时会被 snakemake 永久暂缓（超预算假死）。并行度不足、假死、多工作流超卖这三类问题催生了主机感知的全局预算自动化。

1. **探测主机**：run.py `detect_host_profile()` 输出 `{hostname, cpus, mem_mb}`，cpus 取 `sched_getaffinity` → cgroup `cpu.max` 配额 → `cpu_count` 的最小值，内存沿用 `_detect_total_mem_mb()`（cgroup `memory.max` → meminfo，取 90%）；
2. **cores auto**：`-t` 不传时 `min(floor(cpus × 0.9), cpus − 1)`（90% 份额，至少给系统留 1 核，见 §1 自查条款）；
3. **多工作流预算切分**：每个 Snakemake 进程的 `--resources mem_mb` 按各工作流输入总量的字节比例分配 `总内存 × 0.90`；字节取不到退回均分；份额下限 `min(总/N, 8 GB)`；
4. **clamp 防死锁**：生成的 raw.json 注入 `resource_caps`，`rule_res` 对 `mem_mb ≤ 预算×0.8`、`threads ≤ cpus` 截断并告警；截断做在返回的 dict 子类上，调用点 fallback（`.get("mem_mb", 32000)`）一并封堵。模块 config（手工拼装，无 caps 键）经 `OMICS_RESOURCE_CAPS` 环境变量获得同一份 caps（CmdUtil 按各进程 `--configfile` 导出）。

实现分布在 run.py、common.smk、CmdUtil.py，不涉及任何 .smk 规则声明。

---

## 4. 数据量感知的资源缩放

**背景**：同一规则面对的数据规模跨数量级（历史 raw.json 实测单文件：DR ≈1.4–1.5 GB，Mutation ≈4.4–5.4 GB，scRNAseq ≈7–12 GB，Fiberseq(Totipotency) ≈58–72 GB，Fiberseq(luanchao) ≈430–460 GB），静态声明二选一——按大取值调度密度低，按小取值 OOM。因此把资源从常量改成**输入规模的函数**：输入大就多给、输入小就少给。Snakemake `resources:` 支持 callable，可访问 `input`/`wildcards`/`attempt`，机制现成。

### 4.1 rule_res 扩展（向后兼容，存量调用点零改动）

```python
def rule_res(name, input_size=None, attempt=1):
    """四步，顺序固定：
    1) 静态合并 resources.default < resources[name]；
    2) 配置了 scale 且传入 input_size 时按公式缩放 mem_mb（§4）；
    3) attempt>1 时 base × 2**(attempt-1)（§5）；
    4) 主机 clamp（§3 第 4 点）——必须在最后，否则放大值绕过预算。"""
```

不传 `input_size` → 跳过第 2 步，行为与静态声明完全一致（存量调用点不用一次改完）。

### 4.2 输入规模怎么量（按规则类别）

| 规则类别 | size_metric | 取法 | 例子 |
|---|---|---|---|
| 比对/预处理（per-sample） | `fastq_mb` | `sum(getsize)/1e6` | star_align、cellranger |
| BAM 后处理 | `bam_mb` | 同上，对 bam 求和 | gatk 系、dedup |
| 聚合/差异分析 | `sample_count` | 样本对/文库数 | merge、DESeq2、combine |
| 图表/报表 | — | 不缩放 | plot/report 规则 |

- 字节数而非 reads 数（reads 需扫 fastq 或维护 meta，成本高），即取即用；
- 文件缺失 → `input_size_mb()` 返回 None → 走静态值（**不能返回 0**，0 会把内存压到下限）。

### 4.3 缩放公式

模型按工具内存机理取**固定开销 + 数据相关**两段（单参数幂律已否决——指数无机理含义，不同工具内存形状不同）：

```
mem_mb = clamp(
    base_mem_mb × (fixed_frac + (1 − fixed_frac) × input_size / ref_input_size),
    min_mem_mb,
    max_mem_mb)
```

| 参数 | 含义 | star_align 示例 |
|---|---|---|
| `base_mem_mb` | 基准内存 = 手调值（基准处不推翻已有调优） | 32000 |
| `ref_input_size` | 基准输入量 = 手调值被验证够用时的规模（定法见下） | 按规则定 |
| `fixed_frac` | 固定开销占比 = "索引/模型占多少内存" | 0.9 |
| `min_mem_mb` / `max_mem_mb` | 下限/上限 | 16000 / 200000 |

`fixed_frac` 按**规则（= 具体子命令）的内存行为**给值，不按工具名——同一工具不同子命令行为完全不同（samtools sort 累积 vs view 流式；GATK MarkDuplicates 随输入涨 vs SplitNCigarReads 基本恒定）：

| 内存行为（按规则判定） | 典型例子 | fixed_frac |
|---|---|---|
| 索引/参考驻留，与输入大小基本无关 | STAR align、hisat2、cellranger count | ≈0.9–1 |
| 内存被用户参数限死 | samtools sort（`-m`） | 1 |
| 流式单遍不积累 | samtools view、cutadapt | ≈0–0.3 |
| 中间结果随输入增长 | GATK MarkDuplicates、MSFragger | ≈0–0.3 |
| 矩阵随样本数增长 | DESeq2、merge 类 | 0（规模=样本数） |
| 图表/报表 | plot/report | —（不缩放） |

拿不准的规则宁可不配 scale（保持静态），只给内存行为讲得清的规则配。

**ref_input_size 定法**（历史 raw.json 实测输入量，单文件）：DR ≈1.4–1.5 GB；Mutation ≈4.4–5.4 GB；scRNAseq ≈7–12 GB；Fiberseq(Totipotency) ≈58–72 GB；Fiberseq(luanchao) ≈430–460 GB——差两个数量级，不存在全局 ref。逐规则：查得到调参数据 → 用它；查不到 → 该 workflow 历史中位数；新规则 → 首个真实数据集起步。取偏只影响效率不影响安全；校准闭环（§6）数据充分后可按实测重新锚定（基准点等值重算）。

代入验证（star_align，fixed_frac=0.9）：

```
输入 0.5 GB → 32000 × (0.9+0.1×0.1) =  29 GB
输入   5 GB → 32000 × (0.9+0.1×1)   =  32 GB ← 基准处 = 手调值
输入  50 GB → 32000 × (0.9+0.1×10)  =  61 GB
输入 500 GB → 32000 × (0.9+0.1×100) = 349 GB → 触上限 200 GB
```

threads 一律不缩放：IO-bound 规则加核收益低，核心预算比内存紧张。**公式形态是可替换的实现细节**：集中在 `rule_res` 一处，实测数据充分后可整体换 `a + b×size` 仿射拟合（两段式是其退化形式）。

**公式的定位要说破**：

1. 声明值不限制进程内存，只影响调度器排布——公式的全部价值是**调度密度**（让并发排得下更多 job），不是"省内存"；
2. 因此它对 fixed_frac≈1 的索引驻留型规则几乎是 no-op（输入涨 100 倍，声明值只动 10%）——这类规则**不必配 scale**，维持静态即可；
3. 实测点攒够（≥3）后，对配了 scale 的规则**直接查表/线性插值**优于拿 fixed_frac 参数化外推——公式的长期定位只是**冷启动插值**（新规则、历史覆盖外的输入规模）。

### 4.4 config schema（`resources.<rule>.scale`，可选）

```json
"star_align": {
  "threads": 12, "mem_mb": 32000,
  "scale": {
    "size_metric": "fastq_mb",
    "ref_input_size": 5000,
    "fixed_frac": 0.9,
    "min_mem_mb": 16000,
    "max_mem_mb": 200000
  }
}
```

五字段全部必填（不设隐式默认）；schema 校验 `0 ≤ fixed_frac ≤ 1`、`min_mem_mb ≤ mem_mb ≤ max_mem_mb`、`size_metric ∈ {fastq_mb, bam_mb, sample_count}`；未配 scale 的规则保持静态，默认零变化。

### 4.5 规则侧接入

```python
rule star_align:
    input: fq = ...
    threads: rule_res("star_align").get("threads", 12)   # 不缩放
    resources:
        mem_mb = lambda wildcards, input, attempt: \
            rule_res("star_align",
                     input_size=input_size_mb(input),
                     attempt=attempt).get("mem_mb", 32000)
```

lambda 在每次作业执行前（含重试）求值，拿到运行时真实 input。

### 4.6 边界行为（定死）

| 情况 | 行为 |
|---|---|
| 输入文件缺失或大小为 0 | `input_size_mb()` 返回 None → 静态值 |
| 配了 scale 但调用没传 input_size | 静态值（与现状一致） |
| scale 参数非法 | schema 在 dry-run 期直接拒绝 |
| 缩放结果超主机预算 | clamp 截断 + warning |
| attempt ≥ 2 | 先 attempt 放大，再 clamp |

---

## 5. attempt 升级与自动重试

**背景**：缩放公式估偏不可避免——若一次 OOM 就要人工改资源重跑，自动化不成立。因此失败自动重试、attempt 内存翻倍升级兜底。

- run.py 默认透传 `--restart-times 2`：失败 job 自动重试，v1 不区分失败原因；
- `attempt>1` 时 `mem_mb = base × 2**(attempt-1)`（封顶 cap），只放大内存，不改规则语义；
- 价值：缩放公式估偏也不会一次失败，代价是重跑该规则；OOM 不会演变成人工改资源后全量重跑。

---

## 6. benchmark 校准闭环

**背景**：声明值手调会随工具升级、数据类型变化逐渐失真，且没有任何回填通道；而 benchmark 已随每次运行记录实测峰值，数据现成。因此建立"实测 → 建议表 → 护栏合入"的持续校准闭环——它同时是 §4 缩放公式参数的唯一可靠来源。

### 6.1 数据积累

- 每条执行型规则带 `benchmark:` 指令，记录 `seconds / max_rss / max_vms`，每 job 每 attempt 一行；
- `calibrate_resources.py` 把各运行的 benchmark 追加进 `resource_history.tsv`（append-only，跨运行跨主机）；
- **OOM 事件单独标记**（`.OOM`）：失败 run 的峰值是删失数据（真实需求 ≥ 观测值），拟合按右删失处理，绝不能拿它压低声明值。

### 6.2 "最佳值"定义（目标函数的落地）

makespan 由两条路径受声明值影响：过低 → OOM 重试浪费时间；过高 → 并发排布变少、流水线拉长。故单规则操作化目标 = **"一次跑过"前提下尽量紧**：历史峰值 **p90 × 安全系数**（默认 1.3，数据充分后可降 1.2）。p90 与安全系数是仅有的两个全局超参。

### 6.3 建议与合入

1. ≥3 个数据点才动；样本不足只出建议、不写回；
2. 输出 `resources_suggest.tsv`（rule、n、实测 p90、建议值、当前声明值、变化幅度、删失数、备注）；
3. `calibrate_resources.py --apply` 按护栏写回 config 的 `resources.<rule>.mem_mb`：n≥3、单次 ±30% 限幅、取整、无变化不动、OOM 删失下界（建议值不得低于任何删失观测 × 安全系数）；config 里没有对应资源键的规则不合入；
4. 每次写回是一次可回滚的 config 变更，依据都在 tsv 里可审计。

### 6.4 节奏与防震荡

季度一次或每 N 次生产运行后；工具升级/新数据类型接入后立即。单次调整限幅 ±30%，连续两次同向才允许超限；样本不足不动。

### 6.5 冷启动与优先级

新规则静态起步；优先收编历史 OOM/重试过的规则与大内存 top 规则。

---

## 7. I/O 维度：背景与问题

前几节把 CPU 与内存维度管住了，I/O 维度一度不在预算体系内。2026-11-27 的现场挤压暴露了缺口：16 个比对进程并发把 sda 压到 `%util=100`，而 CPU 还有 36% idle——CPU 有余量、磁盘饱和，正是 I/O 维度缺失的形态。以下先梳理代码与现场事实（7.1），再定位问题（7.2）；后续三节分别回应：I/O 带宽预算（§8）与临时写入分流（§9）直接作用于 makespan，I/O 校准扩展（§10）保证参数长期收敛。

### 7.1 现状梳理（代码与现场事实）

以下均为实际读取的代码/运行数据，不是记忆：

1. **调度只预算 CPU 与内存，无 I/O 维度**。`run.py:425-480` `build_snakemake_cmd()` 组装的 snakemake 命令只有 `--cores <threads>`（run.py:474-475）和 `--resources mem_mb=<mem_mb>`（run.py:479-480）。Snakemake 9.22（本机 `smk` 环境实测 `--help`）支持任意命名资源 `--resources NAME=INT`，机制现成，只是没有被用于 I/O。
2. **per-rule 资源体系已存在且可直接承载新维度**。`modules/common/common.smk:285` `rule_res(name, input_size, attempt)` 四步：静态合并 → scale 缩放（common.smk:248）→ attempt 翻倍 → 主机 clamp（`_RuleRes`，common.smk:266-282，mem_mb/threads 在 get/__getitem__ 上截断）。`resources.<rule>` config 段已覆盖 16 个 workflow（如 `config/RNAseq.json:600-610` `star_align` 的 threads/mem_mb/scale）。`io_bandwidth_mb_s` 作为新命名键走同一条合并/截断链，`rule_res` 主体零改动。
3. **多工作流预算切分已有先例**。`run.py:98` `detect_host_profile()`、`run.py:120` `_workflow_input_bytes()`、`run.py:155` `_split_mem_budget()`（输入字节加权、取不到退均分、份额下限），调用点在 run.py:1263-1296（`_inject_resource_caps` 双通道注入）。I/O 预算的多工作流切分复用同一权重口径。
4. **benchmark 已记录逐规则 I/O 量，校准数据源现成**。实测 `results/46/RNAseq/log/sample/ciTotiSC7/GRCm39/star_align.benchmark.tsv`：表头含 `io_in, io_out` 列；该次 star_align 墙钟 554 s、io_out 16725 MB ≈ **30 MB/s 写需求**（不含输入读）。`src/common/util/calibrate_resources.py` 已消费同一批 tsv 做 mem/threads 校准（只出建议表，人工合入）。
5. **STAR 临时排序写入与最终输出同盘**。`modules/star/star.smk:201`（align 规则）`outTmpDir = outdir + "/{sample_id}/tmp_star"`、star.smk:41（index 规则）同款——tmp 与最终 BAM 都在 `/data/16T`（sda，机械盘）。而 `run.py:501-511` 已有 `tmp_dir` config 键并已加入 singularity bind，`/tmp` 本身在 NVMe（`/dev/nvme0n1p2`），分流通道部分现成。
6. **现场挤压证据**（2026-11-27 实测）：16 个比对进程（STAR ×8 + HISAT2 ×8，各 12 线程）并发时 sda `%util=100`、aqu-sz 10~11、w_await 75~125 ms、聚合写 76~89 MB/s；主机 160 核、load 139、CPU idle 36%——CPU 有余量而磁盘饱和，正是 I/O 维度缺失的形态。
7. **命名规则约束**（既有约定）：`rule_res` 按 snakemake 规则名解析，一个规则服务多个流程阶段时（`star_align` 被 TEtranscripts/featureCounts/SNP/fusion 复用）只允许单键 `star_align`，禁止 per-stage 键。

### 7.2 问题分析

| # | 问题 | 对目标（makespan）的影响 | 对应环节 |
|---|------|--------------------------|----------|
| P1 | I/O 无预算：并发比对的带宽需求之和（8×30 MB/s 量级）远超 HDD 可持续带宽（~100 MB/s），全员排队互拖，单 job 墙钟被拉长 | 直接反作用于目标：并发≠吞吐，越过磁盘饱和点后加 job 只加等待 | I/O 带宽预算（§8） |
| P2 | STAR 临时排序写入与最终 BAM 同盘（star.smk:201）：写放大 + 寻道交错，加剧 P1 | 每 job 的 I/O 需求被抬高，饱和点更早到来 | 临时写入分流（§9） |
| P3 | 校准闭环（§6）只有 mem/threads 维度，无 I/O 建议列；权重/预算无实测数据回填通道 | 参数只能拍脑袋，无法收敛 | I/O 校准扩展（§10） |
| P4 | 多工作流并发时内存有切分（run.py:155）而 I/O 无切分，两个 snakemake 进程各自无上限地压同一块盘 | 多工作流场景 P1 复现且更难排查 | I/O 带宽预算（§8） |

问题按对目标的贡献排序：P1 是主因（现场 makespan 损失的直接来源），P2 降低每 job 需求（与 P1 互补），P3/P4 保证机制长期正确。

---

## 8. I/O 带宽预算（per-rule 声明 + 全局上限，回应 P1/P4）

规则按实测带宽需求声明 `io_bandwidth_mb_s`，run.py 按目标磁盘的可持续带宽给出全局上限 `disk_bandwidth_budget_mb_s`，调度器打包时保证并发 job 的带宽需求之和不超预算，从源头消除"16 个比对同时压一块 HDD"式的 I/O 挤压。

### 8.1 机制与作用通道

Snakemake 的调度打包规则：运行中 job 的各命名资源之和 ≤ `--resources` 给出的全局上限，超出则暂缓调度。把"带宽需求"建模为命名资源，即得磁盘限流：

- 规则侧：`resources: io_bandwidth_mb_s = <该 job 的预期带宽需求，MB/s>`
- 全局侧：`--resources disk_bandwidth_budget_mb_s = <目标磁盘的可持续带宽预算>`

调度器保证任意时刻 Σ(io_bandwidth_mb_s) ≤ disk_bandwidth_budget_mb_s。**作用通道 = 调度打包**，不给进程限速（§1）。一个 job 的实际带宽超出声明值时不会被惩罚——该旋钮只管"多少 job 同时跑"，声明值的准确性影响打包质量而非正确性（低估 → 仍可能挤，高估 → 保守空转，两侧都安全）。

### 8.2 资源命名与语义

| 名称 | 语义 | 声明位置 |
|------|------|----------|
| `io_bandwidth_mb_s` | 该 job 的预期 I/O 带宽需求（io_in+io_out 口径），整数 MB/s | 规则 `resources:` + `config["resources"][<rule>]` |
| `disk_bandwidth_budget_mb_s` | 全局上限：目标磁盘的带宽预算，整数 MB/s | run.py → `--resources`（不进 config） |

语义选"带宽需求 MB/s"而非"无量纲权重"：机制自明（预算与需求同量纲，Σ≤预算即打包约束），且 benchmark 的 (io_in+io_out)/墙钟 可直接换算成声明值（§8.3），校准闭环无需换算系数。规则按**规则名**单键声明（§7.1 第 7 条）。

### 8.3 上限与声明值的确定（公式 + 参数表 + 数值算例）

**全局上限**（run.py 自动探测，一个 helper，公式形状为可替换实现细节）：

```
disk_bandwidth_budget_mb_s = io_safety_factor × sustainable_disk_bandwidth_mb_s
```

| 参数 | 含义 | 取值/来源 | 示例值 |
|------|------|-----------|--------|
| `io_safety_factor` | 给系统/他人进程留的余量 | 暂定 0.8（待校准，由校准扩展（§10）实测重定） | 0.8 |
| `sustainable_disk_bandwidth_mb_s` | 目标盘可持续混合读写带宽 | v1 按盘类型暂定常量（下表）；v2 可选 `--io-probe` 实测；v3 校准回填 | 150 |
| 目标盘判定 | 输出目录所在设备 | `os.stat(outdir).st_dev` → `/sys/dev/block/<maj:min>/queue/rotational` + 设备名含 `nvme` | sda → HDD |

盘类型暂定常量（**全部为待校准暂定值**，数据来源：现场 sda 饱和时聚合吞吐 76~100 MB/s，留余量取 150；SSD/NVMe 为数量级估值，首次真实运行后由校准扩展（§10）重定）：

| 盘类型 | 判定 | sustainable_disk_bandwidth_mb_s（暂定） |
|--------|------|------------------|
| HDD（rotational=1） | queue/rotational = 1 | 150 |
| SATA SSD（rotational=0，非 nvme） | rotational = 0 且设备名无 nvme | 450 |
| NVMe | 设备名含 nvme | 1200 |

**规则声明值**（确定配方，锚定实测而非拍脑袋）：

```
io_bandwidth_mb_s(<rule>) = p90 over benchmarks of (io_in + io_out) / wall_seconds
```

数据源：该规则历史 `*.benchmark.tsv` 的 `io_in/io_out/s` 列（§7.1 第 4 条，校准闭环的 resource_history 同一批数据）。锚定规则：取 p90（稳健上界），≥3 个数据点才给正式值、否则用调用点 fallback 常量；输入规模跨数量级的规则（如 gatk 系）暂用静态值并标注 v1 近似，升级路径 = 复用 scale 机制（与 mem 同款两段式）或按输入字节线性外推，待校准扩展（§10）的数据决定。

**数值算例**（跨 2-3 个数量级，展示打包行为与钳位点；权重为 2026-11-28 落地值，锚定历史 benchmark p90）：

规则权重（实现取值与依据）：

| 规则 | io_bandwidth_mb_s | 依据 |
|------|-------------------|------|
| DESeq2 / 纯计算 | 0 | 无磁盘读写（读的是内存中的小矩阵） |
| bowtie2_align_paired / _single | 3 | benchmark p90 = 3.0（n=39） |
| hisat2_align | 6 | benchmark p90 = 5.6（n=30） |
| MarkDuplicates | 11 | benchmark p90 = 10.4（n=50） |
| star_align | 56 | benchmark p90 = 55.6（n=48，全仓 results 汇总） |
| bam_sort（排序类） | 96 | 无实测；排序类保守 fallback（写≈2×输入），待校准回填 |

打包行为（预算 = 0.8 × 盘带宽）：

| 场景 | 预算 | 打包结果 |
|------|------|----------|
| HDD，2×star_align + hisat2_align | 120 | Σ=2×56+6=118 ≤ 120，3 个 job 并发成立 |
| HDD，第 3 个 star_align 候选 | 120 | Σ=174 > 120，暂缓——这正是 P1 的钳位点（对比现状：16 个全上） |
| HDD，1×bam_sort + 1×star_align | 120 | Σ=152 > 120，bam_sort 独占运行 |
| HDD，声明值 200 的超大 job | 120 | 触发 §8.5 clamp → 压到 120，单独运行（防假死） |
| NVMe 输出 | 960 | 17×star_align 理论可并发，钳位永不触发（限流自动退化为 no-op，正确行为） |
| SATA SSD 输出 | 360 | 6×star_align，通常不绑定，但挡住病态并发 |

### 8.4 run.py 接入

实现项（增量，不动既有函数名——`_detect_total_mem_mb、_collect_bind_paths、_detect_singularity、_merge_singularity_args` 四个被 run_module.py 导入的私有名冻结）：

1. 新增 `_detect_disk_budget_mb_s(outdir)`：按 §8.3 表探测盘类型 → 预算；探测失败回退 HDD 常量并 log 告警（不抛错）。
2. `build_snakemake_cmd()`（run.py:425）追加参数 `disk_budget_mb_s: Optional[int] = None`，非 None 时：

```python
    if mem_mb or disk_budget_mb_s:
        resources = []
        if mem_mb:
            resources.append(f"mem_mb={mem_mb}")
        if disk_budget_mb_s:
            resources.append(f"disk_bandwidth_budget_mb_s={disk_budget_mb_s}")
        cmd += ["--resources", *resources]
```

3. 多工作流切分：run.py:1263 的 `weights = [p["input_bytes"] ...]` 之后，对 disk_budget 用同一权重经 `_split_mem_budget` 同款逻辑分配（新增轻量 `_split_disk_budget(total_mb_s, weights)`，复用份额下限与退均分规则；不改名 `_split_mem_budget`）。每 snakemake 进程拿到各自份额，P4 解决。
4. `--io-budget-mb-s` CLI 旋钮（覆盖自动探测；`0` = 显式关闭限流）与 `run_summary.json` 落盘实际预算值。
5. 默认策略：自动探测开启。兼容性由 §8.7 的权重地板保证——未声明规则只占 1 MB/s，在任何现实并发数下不绑定。

### 8.5 clamp 防死锁（复用 clamp 双通道）

既有故事（§3）：mem_mb 超预算的 job 会被 snakemake 永久暂缓（超预算假死），靠 `_RuleRes` 截断解决。I/O 维度同构：**声明值 ≥ 全局预算的 job 永远凑不齐打包条件**。处理：

- `common.smk` `_clamp_res_value`（common.smk:266-282 调用链）扩展一项：`io_bandwidth_mb_s` 截断到 `disk_bandwidth_budget_mb_s`（caps 值），截断即告警日志；
- caps 双通道注入沿用：raw.json `resource_caps` + `OMICS_RESOURCE_CAPS` 环境变量（run.py:1266 `_inject_resource_caps`），模块 config 侧自动生效；
- run_module.py 独立运行无 caps → 无 clamp 也无预算 → 不限流（与 mem clamp 的现状语义一致）。

### 8.6 边界行为（定死，实现期不得即兴）

| 情形 | 行为 |
|------|------|
| 规则未声明 `io_bandwidth_mb_s` | `--default-resources io_bandwidth_mb_s=1` 地板权重；16 并发 ×1 = 16 ≪ 120，永不绑定——**兼容锚点，未声明规则行为不变** |
| 声明值 = 0 | 允许，纯计算规则不占预算 |
| 声明值 < 0 或非整数 | config load / scale 校验（`_validate_scale` 同款）报错，命名键指向规则 |
| 声明值 ≥ 全局预算 | §8.5 clamp 到预算值，job 单独运行，log 告警 |
| 全局预算探测失败 | 回退 HDD 暂定常量 + 告警，绝不拒绝启动 |
| run_module.py 独立运行 | 无上限（snakemake 对无 `--resources` 上限的命名资源不限流），与 mem 无 clamp 语义对齐 |
| 多工作流并发 | 预算按输入字节加权切分（§8.4 第 3 点），各自进程内打包 |
| attempt 重试 | `io_bandwidth_mb_s` **不**随 attempt 翻倍（带宽需求与尝试次数无关，区别于 mem 的 attempt 翻倍（§5）） |
| scale 配置存在 | v1 不缩放 io 声明值（静态），升级路径见 §8.3（与 mem 的 scale 机制共存不冲突） |
| dry-run | 声明值出现在资源显示中，DAG 与解析行为不变 |

### 8.7 向后兼容（机制，不是口号）

锚点有三：(1) 未声明 `io_bandwidth_mb_s` 的规则权重地板 1，在 ≤120 并发内永不绑定，行为逐比特不变；(2) `--io-budget-mb-s 0` 可整体关闭，回到今天的确切行为；(3) `rule_res` 现有 417 处调用点零改动（新键走既有合并链，`.get("io_bandwidth_mb_s", fallback)` 模式与 mem 完全同款，fallback = JSON 优化值的既有约定照搬）。

---

## 9. 临时写入分流（STAR outTmpDir → 临时目录，默认 /tmp，回应 P2）

STAR 排序产生 ~16.7 GB 临时写再合并成最终 BAM（§7.1 第 4/5 条），目前 tmp 与最终输出同在 sda——写放大且寻道交错，抬高了每个 job 的 I/O 需求（P2）。`tmp_dir` config 键与 singularity bind 已存在（run.py:501-511），`/tmp` 在 NVMe：

实现项：

1. `modules/star/star.smk:201`（align）与 star.smk:41（index）的 `outTmpDir` 改为指向临时目录：`tmp_dir` config 有值时用 `tmp_dir + "/star_align/{sample_id}"`，未配置时**默认 `/tmp/star_align/{sample_id}`，自动分流、默认生效**（2026-11-27 已决：默认生效优先于零行为变化；`tmp_dir` 显式配置可完全覆盖，包括配回 outdir 下旧路径以复现旧行为）；
2. 同款检查扩展到其他大临时写规则（gatk SortSam 的 `--TMP_DIR`、samtools sort `-T`），逐规则排查，不按工具名一把抓（每个子命令的临时行为不同）；
3. tmp 带宽占用仍计入带宽预算（§8）的权重（io_out 实测本来就是全量口径，分流后 benchmark 数字自然回落，权重随之校准下调——带宽预算与写入分流由校准扩展（§10）闭环）。

边界：tmp 盘满 → job 失败走既有重试（§5）；tmp_dir 与 outdir 跨设备时 snakemake 的 `os.replace` 不跨盘——STAR 合并输出直接写 outdir（STAR 本身即 outTmpDir 合并到 outFileNamePrefix），不触发跨盘 move，实现时逐规则确认无 `os.replace` 跨盘调用。**默认 /tmp 的容量风险**：大并发 STAR 的 tmp 可达百 GB 量级，/tmp 在小根分区主机上有写满风险——job 级失败由重试兜底，run.py 落 tmp 前应检查可用空间并告警（阶段 2 实现项）。

---

## 10. I/O 校准扩展（回应 P3）

校准闭环（§6）目前只有 mem/threads 维度；I/O 侧的权重与预算同样是"待校准暂定值"，需要同一套回填通道（P3）：

1. `src/common/util/calibrate_resources.py`：建议表 `resources_suggest.tsv` 增加 `io_bandwidth_mb_s` 列 = p90((io_in+io_out)/s)，沿用既有护栏（§6.3：≥3 点、±30% 限幅、只出建议人工合入）；
2. 可选 `--io-probe`（默认关，2026-11-27 已决保留）：实测 `sustainable_disk_bandwidth_mb_s` 覆盖暂定常量。实现时**优先复用现成 Linux 测速工具（fio、hdparm 等）**直接取磁盘顺序读写带宽，替代手写 256 MB 测速文件（3 次取中位）；现成工具够用则 probe 退化为薄封装。默认关的理由：测速本身写盘、并发运行时测速失真、对用户磁盘有副作用——只在空闲时人工开启；
3. `io_safety_factor` 与盘类型暂定常量随试点实测数据重定（试点方案见行动方案）。

---

## 11. 设计取舍

所有环节共用一套取舍记录（含否决的替代——否决理由也是设计的一部分）：

| 议题 | 决议 | 理由 |
|---|---|---|
| 规模度量 | 字节数（reads 数否决） | 扫 fastq / 维护 meta 成本高 |
| 多工作流内存分割 | 按输入字节加权（均分否决） | 字节是"工作量大头"的粗代理；风险见 §1 自查条款 |
| threads 缩放 | 不缩放，维持静态 | IO-bound 规则加核收益低，核心预算比内存紧张 |
| per-host 覆盖层 | 否决 | 无法枚举所有服务器；异构差异由 clamp 兜底 |
| `--default-resources` | 否决 | config 已有 `resources.default` 两层口径，重复兜底掩盖漏配 |
| 手调较佳值 | 保留为基准点 | 由校准闭环（§6）按实测数据重标定 |
| 缩放公式形态 | 两段式（固定开销+数据相关）；单参数幂律否决 | 幂律指数无机理含义，不同工具内存形状不同 |
| I/O 维度的机制载体 | Snakemake 命名资源 + `--resources` 上限 | 调度器原生打包语义，零侵入规则语义（§7.1 第 1 条实测确认 9.22 支持） |
| I/O 量纲 | 带宽 MB/s，不做无量纲权重 | benchmark 直接换算，机制自明（§8.2） |
| I/O 权重键粒度 | 每规则单键（star_align 不分阶段） | 既有命名规则（§7.1 第 7 条） |
| I/O 默认策略 | 预算自动探测开启 + 权重地板 1 | 兼容锚点（§8.7），且修复默认生效 |
| tmp 分流默认策略 | 默认 `/tmp` 自动分流，未配置即生效（`tmp_dir` 可覆盖，含配回旧路径） | 机制默认生效优先于零行为变化；/tmp 容量风险由重试 + 落盘前空间检查兜底（2026-11-27 已决） |
| 动态反馈调速（外部控制器按 await 实时调并发） | 否决 | Snakemake 无运行时资源重配通道，需常驻外部控制器，效应通道复杂且收益不明 |
| cgroup io.max 进程级限速 | 否决（非目标） | 作用通道是限速而非打包，会让 job 变慢而不是变少，与目标反向 |
| per-host 磁盘配置表 | 否决 | 违背"不枚举部署目标"；自动探测 + 兜底常量在任何主机成立 |
| 启动时强制磁盘测速 | 部分否决（保留可选项） | 测速本身写盘、并发时失真；改为 `--io-probe` 默认关 |

---

## 12. 实现文件

| 文件 | 职责/改动 | 状态 |
|---|---|---|
| `run.py` | detect_host_profile、cores auto、预算加权切分、resource_caps 注入、run_summary.json | 已实现 |
| `modules/common/common.smk` | rule_res 四步合并、clamp（dict 子类）、`input_size_mb`、`_validate_scale` | 已实现 |
| `src/common/util/CmdUtil.py` | `OMICS_RESOURCE_CAPS` 环境通道（模块级 clamp） | 已实现 |
| `src/common/util/calibrate_resources.py` | benchmark 收集、历史累积、建议表、`--apply` 护栏合入 | 已实现 |
| `config/*.json` + `*.schema.json` | `resources` 段（静态声明 + 可选 `scale` 段） | 已实现 |
| `run.py` | `_detect_disk_budget_mb_s`、`build_snakemake_cmd` 追加 `disk_budget_mb_s`、`_split_disk_budget`、`--io-budget-mb-s`、run_summary 落盘 | 已实现（§8） |
| `modules/common/common.smk` | `_clamp_res_value` 扩展 io 键截断（含 `_validate_scale` 同款取值校验） | 已实现（§8） |
| `modules/star/star.smk`、`modules/hisat2/*.smk`、`modules/bowtie2/*.smk`、gatk/samtools 排序规则 | `resources: io_bandwidth_mb_s = rule_res("<rule>").get("io_bandwidth_mb_s", <fallback>)` 声明 | 已实现（§8） |
| `modules/star/star.smk`（:41、:201）及 gatk/samtools 排序规则 | outTmpDir / TMP_DIR 分流到 tmp_dir | 已实现（§9） |
| `config/<Workflow>.json` + `.schema.json` | `resources.<rule>.io_bandwidth_mb_s` 字段（schema 按 per-rule-resources-pattern 的 rule_schema 形态，indent=2 写回） | 已实现（§8） |
| `src/common/util/calibrate_resources.py` | io 建议列 + `--io-probe` | 已实现（§10） |
