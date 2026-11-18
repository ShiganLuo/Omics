# 资源消耗平衡机制改进方案

[← 返回 README](../README.md)

**一句话**:把静态资源体系(每规则一组手调常量,不感知数据量与主机)改造成
"主机感知 + 数据量感知 + 失败兜底 + 实测校准"的自适应体系,
目标是在真实生产负载上最小化整体运行时间(makespan)。

**当前状态**:阶段 1(方案 A)已实施未验收(基线未跑,见 §6);
阶段 2-4 未实施。唯一待定项见 §8。

---

## 1. 设计目标(总纲)

**目标:在真实生产负载上最小化整体运行时间(makespan)。**

单一目标函数,所有子设计(cores auto / 预算切分 / clamp / 缩放 / 校准)
都为它服务。优先级(字典序):

1. **目标**:makespan 最短;
2. **约束**(必须满足,不是优化对象):无重试耗尽导致的最终失败、
   无人工改资源重跑、不破坏规则语义与向后兼容(未超预算的值与现状一致);
3. **代理指标**(仅当目标难以直接观测时参考):内存浪费率、核利用率——
   二者改善不必然改善 makespan,冲突时以实测 makespan 为准。

**非目标**(明确不做):多用户公平调度;集群/SLURM 级调度(本方案只管
单机 snakemake 进程);无护栏的自动改写 config(合入只走
calibrate_resources.py --apply 的护栏:n≥3、±30% 限幅,见 §5.4);
追求声明值"精确"本身(只在影响 makespan 的范围内校准)。

**由目标推出的自查条款**:

- 一切魔法数须对目标给出理由或显式标注为待校准近似:clamp 的 0.8 余量、
  cores 的 −2/上限 256、份额下限 8 GB、p90×1.3(§5.4 D.2)——目前均为经验
  近似,方案 D 数据到位后按 makespan 实测重定;
- 多工作流按**输入字节**加权是"工作量大头"的粗代理,不等于关键路径;
  存在"小数据但串行深的工作流被饿死 → 总时长反升"的理论风险。数据到位后
  应改为按待跑 job 数/实测耗时加权(同一 benchmark 数据源),字节加权作为
  v1 近似保留并监控。

---

## 2. 现状(代码事实,2026-11-16)

### 2.1 已完成的前置修复

- 16 个 workflow 的 config/schema 均已补齐 `resources` 段
  (值=各规则 .smk 中的 fallback 默认值),15 个 subworkflow 共 152 处
  模块 config 补传 `resources` 管道;
- `KARRseq.smk` 的本地 `rule_res` 拷贝已删除,统一 include `common.smk`;
- 6 处未走 `rule_res` 的 threads 声明已统一。

### 2.2 仍存在的静态体系

- `common.smk` 的 `rule_res(name)`:`config["resources"]["default"] <
  config["resources"][name]` 两层静态合并,规则里写死用法:
  ```python
  threads: rule_res("trgt_genotype").get("threads", 4)
  resources:
      mem_mb = rule_res("trgt_genotype").get("mem_mb", 32000)
  ```
  所有值是**按"一份典型数据给到刚好不慢"手调的常量**,与输入大小无关(P1)。
- run.py:`-t/--threads` 默认 10,160 核主机只用 6%(P2);内存预算已自动
  探测(cgroup → meminfo,取 94%)但 per-rule 声明值不感知;没有
  `--restart-times`、没有 `benchmark:`、没有任何运行历史(P4/P5)。
- 多 workflow 并行(`-w` 传多个)各起独立 Snakemake 进程,内存预算按整机
  各发一份,N 个并发即超卖(P7)。
- 流水线跑在至少 4 台异构主机上(zhang_c2 / math_zhou / lv / lv_zheda),
  一套 config 静态值共用,单作业声明值可超小主机预算 → 永久 pending(P3)。

---

## 3. 问题清单(按对 makespan 的预期收益排序)

| # | 问题 | 对 makespan 的影响 | 状态 |
|---|------|------|------|
| P2 | `--cores` 固定 10 | 并行度锁死,160 核闲置 94%,量级最大 | 阶段 1 已修(cores auto) |
| P3 | 单作业声明值可超全局预算 | 永久 pending,流水线假死 | 阶段 1 已修(clamp) |
| P7 | 多 workflow 各领整机内存预算 | 超卖 → swap/OOM,全流程拖慢 | 阶段 1 已修(预算切分) |
| P1 | 规则资源与数据量无关 | 小数据挤占调度密度;大数据 OOM/变慢 | 阶段 2(方案 B) |
| P4 | 无 attempt 升级/自动重试 | OOM 一次即失败,人工全量重跑 | 阶段 3(方案 C) |
| P5 | 无 benchmark | 资源值无实测依据,无法校准 | 阶段 4(方案 D) |
| P6 | 资源口径分散 | 工程整洁问题,不直接影响 makespan | 已修完(§2.1) |

---

## 4. 方案总览

| 方案 | 一句话 | 实现成本 | 解决 | 状态 |
|---|---|---|---|---|
| A 主机感知 | 探测主机,cores auto、预算切分、clamp 防死锁 | 低(1-2 天) | P2 P3 P7 | **已实施** |
| B 数据感知 | 资源 = 输入规模的函数,逐规则配 scale | 中 | P1(核心) | 试点未开始 |
| C attempt 兜底 | OOM 识别 + 自动加倍重试 | 低-中 | P4 | 未开始 |
| D benchmark 闭环 | 实测峰值 → 建议表 → 人工合入,持续校准 | 中-高 | P5,校准 B | 未开始 |

四者可组合;B 的公式参数唯一可靠来源是 D 的实测数据,C 兜住 B 估偏的风险。

---

## 5. 方案细节

### 5.1 方案 A:主机感知 + 全局预算自动化(已实施)

1. **探测主机**:run.py `detect_host_profile()` 输出 `{hostname, cpus, mem_mb}`,
   cpus 取 `sched_getaffinity` → cgroup `cpu.max` 配额 → `cpu_count` 的最小值,
   内存沿用 `_detect_total_mem_mb()`(cgroup `memory.max` → meminfo,94%);
2. **cores auto**:`-t` 不传时 `min(cpus − 2, 256)`(经验近似,见 §1 自查条款);
3. **多工作流预算切分**:每个 Snakemake 进程的 `--resources mem_mb` 按各工作流
   输入总量的字节比例分配 `总内存 × 0.94`;字节取不到退回均分;
   份额下限 `min(总/N, 8 GB)`(v1 近似,见 §1 自查条款);
4. **clamp 防死锁(P3)**:生成的 raw.json 注入 `resource_caps`,
   `rule_res` 对 `mem_mb ≤ 预算×0.8`、`threads ≤ cpus` 截断并告警;
   截断做在返回的 dict 子类上,调用点 fallback(`.get("mem_mb", 32000)`)
   一并封堵。模块 config(手工拼装,无 caps 键)经
   `OMICS_RESOURCE_CAPS` 环境变量获得同一份 caps(CmdUtil 按各进程
   `--configfile` 导出)。

验证记录见附录进展 5。改动点:run.py、common.smk、CmdUtil.py;
不改任何 .smk 规则声明。

### 5.2 方案 B:数据量感知的资源缩放(核心方案,试点未开始)

把资源从常量改成**输入规模的函数**:输入大就多给、输入小就少给(P1)。
Snakemake `resources:` 支持 callable,可访问 `input`/`wildcards`/`attempt`。

#### 5.2.1 rule_res 扩展(向后兼容,存量调用点零改动)

```python
def rule_res(name, input_size=None, attempt=1):
    """四步,顺序固定:
    1) 静态合并 resources.default < resources[name](现状,不变);
    2) 配置了 scale 且传入 input_size 时按公式缩放 mem_mb;
    3) attempt>1 时 base × 2**(attempt-1)(方案 C);
    4) 主机 clamp(方案 A 第 4 点)——必须在最后,否则放大值绕过预算。"""
```

不传 `input_size` → 跳过第 2 步,行为与现状完全一致(279 个引用点不用一次改完)。

#### 5.2.2 输入规模怎么量(按规则类别)

| 规则类别 | size_metric | 取法 | 例子 |
|---|---|---|---|
| 比对/预处理(per-sample) | `fastq_mb` | `sum(getsize)/1e6` | star_align、cellranger |
| BAM 后处理 | `bam_mb` | 同上,对 bam 求和 | gatk 系、dedup |
| 聚合/差异分析 | `sample_count` | 样本对/文库数 | merge、DESeq2、combine |
| 图表/报表 | — | 不缩放 | plot/report 规则 |

- 字节数而非 reads 数(reads 需扫 fastq 或维护 meta,成本高),即取即用;
- 文件缺失 → `input_size_mb()` 返回 None → 走静态值(**不能返回 0**,
  0 会把内存压到下限)。

#### 5.2.3 缩放公式

模型按工具内存机理取**固定开销 + 数据相关**两段(单参数幂律已否决——
指数无机理含义,不同工具内存形状不同):

```
mem_mb = clamp(
    base_mem_mb × (fixed_frac + (1 − fixed_frac) × input_size / ref_input_size),
    min_mem_mb,
    max_mem_mb)
```

| 参数 | 含义 | star_align 示例 |
|---|---|---|
| `base_mem_mb` | 基准内存 = 现有手调值(基准处不推翻已有调优) | 32000 |
| `ref_input_size` | 基准输入量 = 手调值被验证够用时的规模(定法见下) | 按规则定 |
| `fixed_frac` | 固定开销占比 = "索引/模型占多少内存" | 0.9 |
| `min_mem_mb` / `max_mem_mb` | 下限/上限 | 16000 / 200000 |

`fixed_frac` 按**规则(= 具体子命令)的内存行为**给值,不按工具名——
同一工具不同子命令行为完全不同(samtools sort 累积 vs view 流式;
GATK MarkDuplicates 随输入涨 vs SplitNCigarReads 基本恒定):

| 内存行为(按规则判定) | 典型例子 | fixed_frac |
|---|---|---|
| 索引/参考驻留,与输入大小基本无关 | STAR align、hisat2、cellranger count | ≈0.9–1 |
| 内存被用户参数限死 | samtools sort(`-m`) | 1 |
| 流式单遍不积累 | samtools view、cutadapt | ≈0–0.3 |
| 中间结果随输入增长 | GATK MarkDuplicates、MSFragger | ≈0–0.3 |
| 矩阵随样本数增长 | DESeq2、merge 类 | 0(规模=样本数) |
| 图表/报表 | plot/report | —(不缩放) |

拿不准的规则宁可不配 scale(保持静态),试点只收编行为讲得清的规则。

**ref_input_size 定法**(历史 raw.json 实测输入量,单文件):
DR ≈1.4–1.5 GB;Mutation ≈4.4–5.4 GB;scRNAseq ≈7–12 GB;
Fiberseq(Totipotency)≈58–72 GB;Fiberseq(luanchao)≈430–460 GB——
差两个数量级,不存在全局 ref。逐规则:查得到调参数据 → 用它;
查不到 → 该 workflow 历史中位数;新规则 → 首个真实数据集起步。
取偏只影响效率不影响安全;方案 D 落地后可重新锚定(基准点等值重算)。

代入验证(star_align,fixed_frac=0.9):

```
输入 0.5 GB → 32000 × (0.9+0.1×0.1) =  29 GB
输入   5 GB → 32000 × (0.9+0.1×1)   =  32 GB ← 基准处 = 现状手调值
输入  50 GB → 32000 × (0.9+0.1×10)  =  61 GB
输入 500 GB → 32000 × (0.9+0.1×100) = 349 GB → 触上限 200 GB
```

threads 一律不缩放(已定,2026-11-16:IO-bound 规则加核收益低,
核心预算比内存紧张)。**公式形态是可替换的实现细节**:集中在 `rule_res`
一处,方案 D 数据到位后可整体换 `a + b×size` 仿射拟合(两段式是其退化形式)。

**公式的定位要说破(评审补充,2026-11-16)**:

1. 声明值不限制进程内存,只影响调度器排布——公式的全部价值是
   **调度密度**(让并发排得下更多 job),不是"省内存";
2. 因此它对 fixed_frac≈1 的索引驻留型规则几乎是 no-op(输入涨 100 倍,
   声明值只动 10%)——这类规则**不必配 scale**,维持静态即可;
3. 方案 D 攒到 ≥3 个实测点后,对试点规则**直接查表/线性插值**优于拿手猜
   fixed_frac 做参数化外推——公式的长期定位只是**冷启动插值**(新规则、
   历史覆盖外的输入规模)。

#### 5.2.4 config schema(`resources.<rule>.scale`,可选)

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

五字段全部必填(不设隐式默认);schema 校验 `0 ≤ fixed_frac ≤ 1`、
`min_mem_mb ≤ mem_mb ≤ max_mem_mb`、
`size_metric ∈ {fastq_mb, bam_mb, sample_count}`;
未配 scale 的规则保持静态,默认零变化。

#### 5.2.5 规则侧接入(试点)

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

lambda 在每次作业执行前(含重试)求值,拿到运行时真实 input。

**试点对象应选线性增长型**(评审修正,2026-11-16:原清单全是索引驻留型,
fixed_frac≈0.9 时公式近似 no-op,验证通过也证明不了价值):
- 线性增长型(主力):gatk MarkDuplicates、大 BAM 的 sort/merge、
  TEcount/TElocal、DESeq2(sample_count 度量)、fragpipe/deepvariant;
- 索引驻留型只留 1 个作对照组(star_align),预期结论是"scale 无增益,
  维持静态"——阴性对照同样是有效结论。

#### 5.2.6 边界行为(先定死)

| 情况 | 行为 |
|---|---|
| 输入文件缺失或大小为 0 | `input_size_mb()` 返回 None → 静态值 |
| 配了 scale 但调用没传 input_size | 静态值(与现状一致) |
| scale 参数非法 | schema 在 dry-run 期直接拒绝 |
| 缩放结果超主机预算 | clamp 截断 + warning |
| attempt=2/3(方案 C 上线后) | 先 attempt 放大,再 clamp |

#### 5.2.7 试点验证标准

2-3 个大小悬殊的数据集(0.5/5/50 GB fastq):
(a) 全程无 OOM;(b) 小数据声明值显著低于静态值;
(c) 资源值随输入单调不减、基准处回归现状值。三条全过再铺开。

### 5.3 方案 C:attempt 升级 + 自动重试(未开始)

- `mem_mb` callable:`base × 2**(attempt-1)`(封顶 cap),run.py 加
  `--restart-times 2`;
- 只对识别为 OOM 的失败升级重试(退出码 137 / 日志含 OOM 字样),
  其他错误照常报错(需小 wrapper 或 Snakemake 侧失败分类);
- 价值:B 公式估偏也不会一次失败,代价是重跑该规则。

### 5.4 方案 D:benchmark 自学习闭环(已落地,每轮运行后校准)

#### D.1 数据积累

- 每条规则加 `benchmark:` 指令,记录 `seconds / max_rss / max_vms`,
  每 job 每 attempt 一行;
- 校准脚本追加进 `resource_history.tsv`(append-only,跨运行跨主机),
  键 = (rule, input_size, attempt),input_size 与 B 同口径;
- **OOM 事件单独标记**:失败 run 的峰值是删失数据(真实需求 ≥ 观测值),
  拟合按右删失处理,绝不能拿它压低声明值。

#### D.2 "最佳值"定义(目标函数的落地)

makespan 由两条路径受声明值影响:过低 → OOM 重试浪费时间;
过高 → 并发排布变少、流水线拉长。故单规则操作化目标 =
**"一次跑过"前提下尽量紧**:历史峰值 **p90 × 安全系数**
(初期 1.3,数据攒多后 1.2)。p90 与安全系数是仅有的两个全局超参。

#### D.3 拟合与建议表

1. ≥3 个数据点才拟合形状(fixed_frac/min/max);1-2 个点只校 base;
   0 个点维持静态;
2. 两段式最小二乘 → fixed_frac 建议值;残差 >±30% → 改仿射
   `a + b×size` 或标注"行为不规则,保持静态";
3. 输出 `resources_suggest.tsv`(rule、当前值、建议值、样本数、输入范围、
   残差),由 `calibrate_resources.py --apply` 按护栏合入 config
   (n≥3、±30% 限幅、取整、无变化不动),git 留痕可回滚;
4. 合入后按 §7 指标验收(makespan 不劣化、浪费率降、OOM=0)。

#### D.4 节奏与防震荡

季度一次或每 N 次生产运行后;工具升级/新数据类型接入后立即。
单次调整限幅 ±30%,连续两次同向才允许超限;样本不足不动;
每次合入是可回滚的 config 变更,依据都在 tsv 里可审计。

#### D.5 冷启动与优先级

新规则静态起步;优先收编历史 OOM/重试过的规则与大内存 top 规则。

---

## 6. 落地路线与当前状态

按 §1 的 makespan 预期收益排序:P2(并行度,量级最大)≫ P3/P7(卡死/超卖)
> P1(配额失配)> P4/P5(兜底与校准)。

| 阶段 | 内容 | 状态 | 验收标准 |
|---|---|---|---|
| 1 = 方案 A | cores auto、预算切分、clamp、统一 rule_res | **代码已落地,验收未完成** | dry-run 无超预算声明;160 核机 auto 用 ≥150 核;并发预算和 = 94% 不超卖;P3 反例被 clamp 且告警 |
| 2 = 方案 B | 试点 8-10 规则配 scale,悬殊数据集实测 | **全量接线完成(417 处声明→lambda;RNAseq 试点 6 规则配 scale),实测验证待真实数据** | §5.2.7 三条 + 小数据声明值显著低于静态、无 OOM |
| 3 = 方案 C | --restart-times + OOM 分类 + attempt 翻倍 | **代码已落地(--restart-times 2 + attempt 翻倍);OOM-only 分类偏差见进展 11** | 注入低估场景自动重试成功,无人工干预 |
| 4 = 方案 D | benchmark 常态化,季度校准 | **代码已落地(376 规则 benchmark + src/common/util/calibrate_resources.py),2026-11-17 RNAseq 首轮校准完成(16 规则合入)** | 建议值 vs 实际峰值偏差中位 ≤20%,建议表可审 |

**⚠ 阶段 1 流程债:改前基线未跑**(§7.1 要求基线先行,代码却已落地)。
补救:改动未 commit,HEAD 仍是旧代码——`git stash` 后跑基线再
`git stash pop`,基线归档后阶段 1 方可视为验收完成。

---

## 7. 效果评判

**改前先跑基线**:用当前(旧)代码把跑分集各跑一遍并归档。没有基线,
"改善"只是感觉。

### 7.1 评测负载(固定跑分集)

- 3 档大小悬殊的数据集 × 2-3 个代表 workflow(建议 0.5 / 5 / 50 GB fastq
  × RNAseq、Fiberseq、scRNAseq),同机同数据;
- **每个配置重复 ≥3 次**(评审补充:makespan 受 conda 建环境、IO 争用、
  缓存状态污染,阶段 2-4 效应量小,单次差异大概率是噪声;报告均值±范围,
  环境预热跑 1 次不计入);
- 每阶段落地后重跑同一跑分集对比,跑分集与数字归档。

### 7.2 指标(§1 总纲的落地)——按阶段声明期望效应量

主指标 = makespan;其余为约束/诊断。**期望效应量按机制分层**(评审补充:
阶段 1 解锁并行宽度,makespan 大幅缩短;阶段 2-4 只提高调度密度与稳健性,
典型 2 样本运行中多数规则仅 2 个 job,调度密度非瓶颈——对其期望 makespan
收益是不现实的):

| 层级 | 指标 | 口径 | 判定 |
|---|---|---|---|
| **主指标** | makespan 比 | 同负载改后/改前总墙钟(≥3 次均值) | **阶段 1**:多工作流并发显著缩短;**阶段 2-4**:不劣化即可,价值以约束项证明 |
| 约束 | OOM/卡死 | OOM 次数 + 永久 pending 数 | 0(违者判失败)——阶段 3/4 的主要价值所在 |
| 约束 | 人工干预 | 手改资源重跑次数 | 0 |
| 诊断 | 内存浪费率 | 1 − Σmax_rss / Σ声明 mem_mb | 方向参考;与 makespan 冲突时以 makespan 为准 |

每阶段收尾:主指标按该阶段期望判定 + 两条约束必须全过,诊断量只解释不挡路。
高并发场景(多工作流、多样本)是阶段 2-4 唯一可能体现 makespan 增益的场合,
跑分集应至少含一个多工作流并发配置。

### 7.3 分阶段验收标准

见 §6 表格第 4 列。

### 7.4 数据来源

- `benchmark:` tsv(方案 D 常态化后全规则覆盖,试点期至少覆盖试点规则);
- Snakemake 日志:job 成败/重试;
- run.py 每次运行输出 run_summary.json:hostname、探测值、各工作流预算、
  makespan、job 统计——已随阶段 1 落地。

---

## 8. 待审问题(仅剩一项)

1. `fixed_frac` 起步值:按规则(子命令)内存行为给(索引驻留型 ≈0.9–1、
   限内存/流式型视情 0–1、增长型 ≈0–0.3)是否可接受?
   注:参数在 config `scale` 段,随时可调、下次运行即生效(dry-run 可核对);
   方案 D 数据到位后按实测重定,起步值只影响效率不影响安全。

---

## 附录 A:已决事项(2026-11-16 审核决议)

| 决议 | 内容 |
|---|---|
| 规模度量 | 字节数(reads 数否决:扫 fastq/维护 meta 成本高) |
| 多工作流内存分割 | 按输入字节加权(均分否决);v1 近似,见 §1 风险条款 |
| threads 缩放 | 不缩放,维持静态 |
| per-host 覆盖层 | 否决(无法枚举所有服务器);异构差异由 clamp 兜底 |
| `--default-resources` | 否决(config 已有 `resources.default` 两层口径,重复兜底掩盖漏配) |
| 手调较佳值 | 先保留为基准点;重标定随方案 D 实测数据在阶段 4 进行 |
| 缩放公式形态 | 两段式(固定开销+数据相关);单参数幂律否决 |

## 附录 B:变更记录

- 2026-11-15 进展 1:KARRseq.smk 本地 rule_res 删除,改 include common.smk;
  其 config/schema 补 resources 段。
- 2026-11-15 进展 2:P6 全部处理完:16 workflow config/schema 补齐
  resources,152 处模块 config 补传管道,6 处散落 threads 统一;
  16 workflow dry-run 全过(301 job 中 288 值不变、13 显式化)。
  顺带发现 9 处 `use rule X` 指向不存在规则(休眠死接线),待议。
- 2026-11-16 进展 3:方案 A 原"per-host 覆盖层"否决删除。
- 2026-11-16 进展 4:决议——字节数、加权切分、threads 不缩放(附录 A)。
- 2026-11-16 进展 5:**阶段 1(方案 A)落地**:run.py 加
  detect_host_profile/cores auto/加权切分/resource_caps 注入/run_summary.json;
  common.smk rule_res clamp(dict 子类封堵调用点 fallback);
  CmdUtil 按进程导出 OMICS_RESOURCE_CAPS(模块级 clamp 通道)。
  验证:18 config × 1904 资源值新旧零差异;16 workflow dry-run 全过
  (Population_genomics/spatial 无测试夹具,与本次无关);
  `--mem-mb 500` 反例 155 处截断告警且通过;双工作流切分与
  caps/env/命令行三处一致;run_module.py 独立运行行为不变。
- 2026-11-16 进展 6:`--default-resources` 否决(附录 A)。
- 2026-11-16 进展 7:手调值保留为基准点(附录 A)。
- 2026-11-16 进展 8:补立 §1 总纲(目标=makespan 最短);缩放公式幂律→
  两段式;ref_input_size 定法=逐规则;D.2 目标函数落地;文档整体重构
  (目标前置、现状更新为当前事实、变更记录移附录)。
- 2026-11-16 进展 9(评审):①公式定位说破——声明值只影响调度密度,
  对索引驻留型近似 no-op,这类规则不配 scale;D 数据足够后以查表/插值
  替代公式,公式退为冷启动插值。②试点对象改为线性增长型(MarkDuplicates、
  sort/merge、TEcount/DESeq2 等),star_align 仅作阴性对照。③§7.1 跑分
  加 ≥3 次重复;§7.2 按阶段声明期望效应量(阶段 1 看 makespan 大幅缩短,
  阶段 2-4 看约束项,makespan 只求不劣化)。
- 2026-11-16 进展 10(**阶段 2 代码落地,试点验证待真实数据**):
  common.smk:`rule_res(name, input_size=None, attempt=1)` 四步扩展
  (静态合并→scale 两段式→attempt 翻倍→clamp 最后),`input_size_mb()`
  (字节数;缺文件/0 字节→None→静态值),`_validate_scale`(dry-run 期拒绝
  非法参数);接线 10 处规则声明(star_align×2 阴性对照、MarkDuplicates×2、
  TEcount×3、TElocal×3,mem_mb 改 lambda,threads 不动);
  config/RNAseq.json 试点 4 规则配 scale 段(fixed_frac:star 0.9 / 其余
  0.2)。验证(未跑真实数据):16 项单元检查全过(锚点回归 32000、下限/
  上限/clamp/attempt/非法参数拒绝);RNAseq dry-run PASS,试点规则 mem_mb
  按输入实际缩放(0 字节输入→静态 32000;dry-run 占位 bam ~78 KB→下限
  8000,属占位文件伪尺寸,真实运行不受影响);Fiberseq/QuantMS dry-run
  回归通过。等价性:不传 input_size 的调用与旧行为逐值一致。
- 2026-11-16 进展 11(**全量改造 + 方案 C/D 代码落地**):
  ①B 全量:全仓库 417 处 `mem_mb = rule_res(...)` 静态声明全部改 lambda
  (`input_size_named(input,'bam'/'fastq'/'vcf')` 命名访问 + 位置输入回退
  整表求和;DESeq2 系用 `sample_count()`,跨 genome_paired/single_samples
  去重);RNAseq.json 增配 DESeq2_TEcount/TElocal(sample_count 度量,
  ref=2 样本,fixed_frac=0)→ 试点共 6 规则。
  ②D 前置:376 个执行型规则加 `benchmark: <log 路径>.benchmark.tsv`
  (与 log 同 wildcards,天然防多 job 写同一文件;128 个聚合/dummy 规则不配);
  `src/common/util/calibrate_resources.py`(D.1 历史 append + D.3 建议:p90×安全系数、
  ≥3 点才动、±30% 限幅、.OOM 标记删失下界;--apply 按护栏合入 config)。
  ③C:`--restart-times`(默认 2)透传 snakemake;attempt 翻倍已在 rule_res。
  偏差记录:OOM-only 重试分类未做(snakemake 无原生按因重试;.OOM 标记
  机制留给 wrapper,当前 v1 对所有失败重试,attempt 升级仅影响 mem_mb)。
  修复过程:benchmark 首版静态路径与 wildcards 冲突→改 log 派生;regex
  插入曾引入 3 类损伤(空行/多行表达式/name= 前缀),全部修复。
  验证:16 workflow dry-run 全过;helper 单元套件全绿(含 sample_count
  去重);calibrate 脚本合成数据测试通过(含删失下界);真实数据验证仍待跑。

## 附录 C:涉及文件

| 文件 | 改动 | 状态 |
|---|---|---|
| `run.py` | detect_host_profile、cores auto、加权切分、resource_caps、run_summary.json | 已改(未 commit) |
| `modules/common/common.smk` | rule_res clamp(dict 子类) | 已改(未 commit) |
| `src/common/util/CmdUtil.py` | OMICS_RESOURCE_CAPS 环境通道 | 已改(未 commit) |
| `subworkflow/KARRseq.smk` | 本地 rule_res 拷贝删除 | 已完成(进展 1) |
| `config/*.json` + `*.schema.json` | resources 段补齐;scale 段待阶段 2 | resources 已完成;**RNAseq.json 4 规则 scale 已配(阶段 2 试点)** |
| 各模块 `.smk` | 试点规则传 `input_size=` | **已接线 10 处(star×2/MarkDup×2/TEcount×3/TElocal×3,未 commit)** |
