# 资源平衡机制设计

[← 返回 README](../README.md)

**一句话**：资源分配是"主机感知 + 数据量感知 + 失败兜底 + 实测校准"的自适应体系，目标是在真实生产负载上最小化整体运行时间（makespan）。

---

## 1. 设计目标（总纲）

**目标：在真实生产负载上最小化整体运行时间（makespan）。**

单一目标函数，所有子设计（cores auto / 预算切分 / clamp / 缩放 / 校准）都为它服务。优先级（字典序）：

1. **目标**：makespan 最短；
2. **约束**（必须满足，不是优化对象）：无重试耗尽导致的最终失败、无人工改资源重跑、不破坏规则语义与向后兼容（未超预算的值与静态声明一致）；
3. **代理指标**（仅当目标难以直接观测时参考）：内存浪费率、核利用率——二者改善不必然改善 makespan，冲突时以实测 makespan 为准。

**非目标**（明确不做）：多用户公平调度；集群/SLURM 级调度（本机制只管单机 snakemake 进程）；无护栏的自动改写 config（合入只走 `calibrate_resources.py --apply` 的护栏：n≥3、±30% 限幅，见 §3.4）；追求声明值"精确"本身（只在影响 makespan 的范围内校准）。

**自查条款**（一切魔法数须对目标给出理由或显式标注为待校准近似）：

- clamp 的 0.8 余量、cores 的 −2/上限 256、份额下限 8 GB、p90×1.3——均为经验近似，由校准闭环（§3.4）按实测数据重定；
- 多工作流按**输入字节**加权是"工作量大头"的粗代理，不等于关键路径；存在"小数据但串行深的工作流被饿死 → 总时长反升"的理论风险。按待跑 job 数/实测耗时加权（同一 benchmark 数据源）是备选口径，字节加权作为 v1 近似保留并监控。

---

## 2. 总体结构

四个机制叠加，共同服务 makespan 目标：

| 方案 | 一句话 | 解决 |
|---|---|---|
| A 主机感知 | 探测主机，cores auto、预算切分、clamp 防死锁 | 并行度不足、声明超预算假死、多工作流超卖 |
| B 数据感知 | 资源 = 输入规模的函数，逐规则配 scale | 静态声明与数据量失配 |
| C attempt 兜底 | 自动重试 + 内存翻倍升级 | OOM 一次即失败 |
| D benchmark 闭环 | 实测峰值 → 建议表 → 护栏合入，持续校准 | 声明值失真，并为 B 锚定参数 |

四者可组合：B 的公式参数唯一可靠来源是 D 的实测数据，C 兜住 B 估偏的风险。

---

## 3. 方案细节

### 3.1 方案 A：主机感知 + 全局预算自动化

1. **探测主机**：run.py `detect_host_profile()` 输出 `{hostname, cpus, mem_mb}`，cpus 取 `sched_getaffinity` → cgroup `cpu.max` 配额 → `cpu_count` 的最小值，内存沿用 `_detect_total_mem_mb()`（cgroup `memory.max` → meminfo，取 94%）；
2. **cores auto**：`-t` 不传时 `min(cpus − 2, 256)`（经验近似，见 §1 自查条款）；
3. **多工作流预算切分**：每个 Snakemake 进程的 `--resources mem_mb` 按各工作流输入总量的字节比例分配 `总内存 × 0.94`；字节取不到退回均分；份额下限 `min(总/N, 8 GB)`；
4. **clamp 防死锁**：生成的 raw.json 注入 `resource_caps`，`rule_res` 对 `mem_mb ≤ 预算×0.8`、`threads ≤ cpus` 截断并告警；截断做在返回的 dict 子类上，调用点 fallback（`.get("mem_mb", 32000)`）一并封堵。模块 config（手工拼装，无 caps 键）经 `OMICS_RESOURCE_CAPS` 环境变量获得同一份 caps（CmdUtil 按各进程 `--configfile` 导出）。

实现分布在 run.py、common.smk、CmdUtil.py，不涉及任何 .smk 规则声明。

### 3.2 方案 B：数据量感知的资源缩放

把资源从常量改成**输入规模的函数**：输入大就多给、输入小就少给。Snakemake `resources:` 支持 callable，可访问 `input`/`wildcards`/`attempt`。

#### 3.2.1 rule_res 扩展（向后兼容，存量调用点零改动）

```python
def rule_res(name, input_size=None, attempt=1):
    """四步，顺序固定：
    1) 静态合并 resources.default < resources[name]；
    2) 配置了 scale 且传入 input_size 时按公式缩放 mem_mb（方案 B）；
    3) attempt>1 时 base × 2**(attempt-1)（方案 C）；
    4) 主机 clamp（方案 A 第 4 点）——必须在最后，否则放大值绕过预算。"""
```

不传 `input_size` → 跳过第 2 步，行为与静态声明完全一致（存量调用点不用一次改完）。

#### 3.2.2 输入规模怎么量（按规则类别）

| 规则类别 | size_metric | 取法 | 例子 |
|---|---|---|---|
| 比对/预处理（per-sample） | `fastq_mb` | `sum(getsize)/1e6` | star_align、cellranger |
| BAM 后处理 | `bam_mb` | 同上，对 bam 求和 | gatk 系、dedup |
| 聚合/差异分析 | `sample_count` | 样本对/文库数 | merge、DESeq2、combine |
| 图表/报表 | — | 不缩放 | plot/report 规则 |

- 字节数而非 reads 数（reads 需扫 fastq 或维护 meta，成本高），即取即用；
- 文件缺失 → `input_size_mb()` 返回 None → 走静态值（**不能返回 0**，0 会把内存压到下限）。

#### 3.2.3 缩放公式

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

**ref_input_size 定法**（历史 raw.json 实测输入量，单文件）：DR ≈1.4–1.5 GB；Mutation ≈4.4–5.4 GB；scRNAseq ≈7–12 GB；Fiberseq(Totipotency) ≈58–72 GB；Fiberseq(luanchao) ≈430–460 GB——差两个数量级，不存在全局 ref。逐规则：查得到调参数据 → 用它；查不到 → 该 workflow 历史中位数；新规则 → 首个真实数据集起步。取偏只影响效率不影响安全；校准闭环（§3.4）数据充分后可按实测重新锚定（基准点等值重算）。

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

#### 3.2.4 config schema（`resources.<rule>.scale`，可选）

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

#### 3.2.5 规则侧接入

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

#### 3.2.6 边界行为（定死）

| 情况 | 行为 |
|---|---|
| 输入文件缺失或大小为 0 | `input_size_mb()` 返回 None → 静态值 |
| 配了 scale 但调用没传 input_size | 静态值（与现状一致） |
| scale 参数非法 | schema 在 dry-run 期直接拒绝 |
| 缩放结果超主机预算 | clamp 截断 + warning |
| attempt ≥ 2 | 先 attempt 放大，再 clamp |

### 3.3 方案 C：attempt 升级 + 自动重试

- run.py 默认透传 `--restart-times 2`：失败 job 自动重试，v1 不区分失败原因；
- `attempt>1` 时 `mem_mb = base × 2**(attempt-1)`（封顶 cap），只放大内存，不改规则语义；
- 价值：方案 B 公式估偏也不会一次失败，代价是重跑该规则；OOM 不会演变成人工改资源后全量重跑。

### 3.4 方案 D：benchmark 校准闭环

#### D.1 数据积累

- 每条执行型规则带 `benchmark:` 指令，记录 `seconds / max_rss / max_vms`，每 job 每 attempt 一行；
- `calibrate_resources.py` 把各运行的 benchmark 追加进 `resource_history.tsv`（append-only，跨运行跨主机）；
- **OOM 事件单独标记**（`.OOM`）：失败 run 的峰值是删失数据（真实需求 ≥ 观测值），拟合按右删失处理，绝不能拿它压低声明值。

#### D.2 "最佳值"定义（目标函数的落地）

makespan 由两条路径受声明值影响：过低 → OOM 重试浪费时间；过高 → 并发排布变少、流水线拉长。故单规则操作化目标 = **"一次跑过"前提下尽量紧**：历史峰值 **p90 × 安全系数**（默认 1.3，数据充分后可降 1.2）。p90 与安全系数是仅有的两个全局超参。

#### D.3 建议与合入

1. ≥3 个数据点才动；样本不足只出建议、不写回；
2. 输出 `resources_suggest.tsv`（rule、n、实测 p90、建议值、当前声明值、变化幅度、删失数、备注）；
3. `calibrate_resources.py --apply` 按护栏写回 config 的 `resources.<rule>.mem_mb`：n≥3、单次 ±30% 限幅、取整、无变化不动、OOM 删失下界（建议值不得低于任何删失观测 × 安全系数）；config 里没有对应资源键的规则不合入；
4. 每次写回是一次可回滚的 config 变更，依据都在 tsv 里可审计。

#### D.4 节奏与防震荡

季度一次或每 N 次生产运行后；工具升级/新数据类型接入后立即。单次调整限幅 ±30%，连续两次同向才允许超限；样本不足不动。

#### D.5 冷启动与优先级

新规则静态起步；优先收编历史 OOM/重试过的规则与大内存 top 规则。

---

## 4. 设计取舍

| 决议 | 内容 |
|---|---|
| 规模度量 | 字节数（reads 数否决：扫 fastq/维护 meta 成本高） |
| 多工作流内存分割 | 按输入字节加权（均分否决）；粗代理，风险见 §1 自查条款 |
| threads 缩放 | 不缩放，维持静态 |
| per-host 覆盖层 | 否决（无法枚举所有服务器）；异构差异由 clamp 兜底 |
| `--default-resources` | 否决（config 已有 `resources.default` 两层口径，重复兜底掩盖漏配） |
| 手调较佳值 | 保留为基准点；由校准闭环（§3.4）按实测数据重标定 |
| 缩放公式形态 | 两段式（固定开销+数据相关）；单参数幂律否决 |

---

## 5. 实现文件

| 文件 | 职责 |
|---|---|
| `run.py` | detect_host_profile、cores auto、预算加权切分、resource_caps 注入、run_summary.json |
| `modules/common/common.smk` | rule_res 四步合并、clamp（dict 子类）、`input_size_mb`、`_validate_scale` |
| `src/common/util/CmdUtil.py` | `OMICS_RESOURCE_CAPS` 环境通道（模块级 clamp） |
| `src/common/util/calibrate_resources.py` | benchmark 收集、历史累积、建议表、`--apply` 护栏合入 |
| `config/*.json` + `*.schema.json` | `resources` 段（静态声明 + 可选 `scale` 段） |
