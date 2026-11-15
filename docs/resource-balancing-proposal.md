# 资源消耗平衡机制改进方案(草案,待审核)

[← 返回 README](../README.md)

现状:各 rule 通过 `rule_res(name)` 获得一组**静态**的 `threads` / `mem_mb`
"较佳值",既不感知输入数据量,也不感知当前主机的 CPU/内存规模。本文分析
问题并给出 4 个候选方案与推荐落地路线。**本文档只是方案设计,尚未改任何
代码,请审核后再实施。**

> 进展(2026-11-15):§2 P6 中的 KARRseq 问题已先行修复——`KARRseq.smk`
> 的本地 `rule_res` 拷贝已删除,改 include `common.smk`;其 4 个模块 config
> 补传 `resources`;`config/KARRseq.json` + `KARRseq.schema.json` 增加
> `resources` 段(取值=原 fallback,行为不变)。其余 9 个无 resources 段的
> 配置与方案 A-D 待讨论。

> 进展 2(2026-11-15):P6 已全部处理完——16 个 workflow 的 config/schema
> 均补齐 `resources`(值=各规则 .smk 中显式写的 fallback 默认值);
> 15 个 subworkflow 共 152 处模块 config 补传 `resources` 管道;
> 6 处未用 `rule_res` 的 threads 声明已统一(censatbuild_mcg_pileup、
> gtfbuild_liftoff、gtfbuild_wga、function_pseudobulk_go_index、
> spatial_scanpy_advanced、fumitools_dedup)。16 个 workflow dry-run 全过,
> 301 个 job 中 288 个资源值不变,13 个变为显式默认值(明细见任务报告)。
> 顺带发现:9 处 `use rule X` 指向模块中不存在的规则(Snakemake 静默忽略,
> 属休眠死接线),待讨论是否修复。

---

## 1. 现状梳理(代码事实)

### 1.1 per-rule 静态资源

- `modules/common/common.smk` 的 `rule_res(name)`:
  `config["resources"]["default"] < config["resources"][name]` 两层静态合并,
  规则里写死用法:
  ```python
  threads: rule_res("trgt_genotype").get("threads", 4)
  resources:
      mem_mb = rule_res("trgt_genotype").get("mem_mb", 32000)
  ```
- `subworkflow/KARRseq.smk` 有一份**本地拷贝**的 `rule_res`(该文件不 include
  common.smk)。全仓库只有这两处定义,扩展时必须同步或统一。
- 16 个 workflow 配置(`config/*.json`)中只有 6 个带 `resources` 段
  (DR=9 条、Fiberseq=49、Mutation=29、PacVar=36、RNAseq=44、scRNAseq=16),
  其余 10 个(CoCulture/MERIP/ncRNAseq/CLIP/PeakCalling/QuantMS/RIPseq/
  tRNAseq/KARRseq/LRtranscriptome)完全靠各 .smk 里 `.get("threads", N)`
  的散落 fallback 值兜底。
- 所有值都是**按"一份典型数据给到刚好不慢"手工调的常量**,与样本输入大小
  无关。

### 1.2 全局预算(run.py)

- `--cores`:`-t/--threads` 直传,默认 **10**(run.py:107,359)。
  当前主机 160 核 → 默认只用 6% 的核。
- 内存预算:`--mem-mb` 不传时 `_detect_total_mem_mb()` 自动探测
  (cgroup v2 `memory.max` → `/proc/meminfo`,取 94%),
  经 `--resources mem_mb=N` 传给 Snakemake(run.py:63-86,364,1120)。
  这一层已经感知主机内存,但 per-rule 声明值不感知。
- 没有 `--restart-times`,没有 `benchmark:`,没有任何运行历史数据。

### 1.3 异构主机

`node.py` 的 `SERVER_DATABASE_DIR` 显示流水线会跑在至少 4 台主机上
(zhang_c2 / math_zhou / lv / lv_zheda),CPU、内存各不相同,
一套 config 的静态值在所有主机上共用。

---

## 2. 问题分析

| # | 问题 | 后果 |
|---|------|------|
| P1 | 规则资源与数据量无关 | 小数据按大配额跑(浪费,挤占调度);大数据按小配额跑(变慢甚至 OOM kill) |
| P2 | `--cores` 固定 10 | 160 核主机 90%+ 闲置;小核主机可能超订 |
| P3 | 单作业声明值可超过全局预算 | Snakemake 永远不调度该作业,表现为"卡死 pending" |
| P4 | 没有 attempt 升级/自动重试 | OOM 一次即失败,只能人工加大再全量重跑 |
| P5 | 没有 benchmark | 每条规则的真实峰值内存/耗时无据可依,资源值无法校准 |
| P6 | `rule_res` 有 KARRseq 本地拷贝、10 个配置无 resources 段 | 资源口径分散,调一处漏一片 |
| P7 | 多 workflow 并行(`-w` 传多个)各起独立 Snakemake 进程 | 内存预算按整机各发一份,N 个工作流并发时可超卖内存 |

---

## 3. 候选方案

四个方案分别对应:**主机感知、数据感知、失败兜底、长期闭环**,可组合。

### 方案 A:主机感知 + 全局预算自动化(成本低,收益立竿见影)

1. **探测 CPU 与内存**:扩展 run.py 现有 `_detect_total_mem_mb()` 为
   `detect_host_profile()`,输出 `{hostname, cpus, mem_mb}`,
   考虑 cgroup(`cpu.max`/`memory.max`)与 `os.sched_getaffinity`
   (容器/任务绑定下比 `os.cpu_count` 准)。
2. **`--cores` 默认值改为 auto**:`-t` 不传时用
   `min(affinity_cpus - 2, 某上限)`;传了则尊重用户值。
3. **多工作流预算分割**:`-w` 跑 N 个工作流时,每个 Snakemake 进程的
   `--resources mem_mb` 按 `总内存 × 0.94 / N` 下发(或按各工作流输入数据量
   加权,见待审问题 4)。
4. **clamp 防死锁(P3)**:`rule_res` 返回前统一
   `mem_mb = min(declared, budget × 0.8)`,`threads = min(declared, cpus)`,
   截断时打 warning。任何主机上都不会出现"单作业永远装不下"。
5. **per-host 覆盖层**(可选):`config["resources"]["host_profiles"][hostname]`
   作为第三层合并,个别主机需要特殊值时不用改全局。

改动点:run.py、common.smk(+KARRseq.smk 拷贝)。不改任何 .smk 规则调用。

### 方案 B:数据量感知的资源缩放(核心方案)

把资源从常量改成**输入规模的函数**。Snakemake 的 `threads:` / `resources:`
支持 callable(lambda),可访问 `input` / `wildcards` / `attempt`,技术上可行。

1. **向后兼容地扩展 helper**:
   ```python
   def rule_res(name, size_mb=None, attempt=1):
       # 1) 静态合并(default < rule)
       # 2) 若配置了 scale 且给了 size_mb,按公式缩放
       # 3) attempt 升级(见方案 C)
       # 4) 主机 clamp(见方案 A)
   ```
   不传 `size_mb` 的调用行为与现在完全一致 → 279 个引用文件不需要一次改完。
2. **规模度量按规则类别**:
   | 规则类别 | 度量 | 取法 |
   |---|---|---|
   | 比对/预处理(per-sample) | 输入 fastq(.gz) 总字节 | `sum(os.path.getsize(f) for f in input)` |
   | BAM 后处理 | bam 字节 | 同上 |
   | 聚合/差异分析(merge、DESeq2、combine) | 样本数 | `len(config["paired_samples"])+…` |
   | 图表/报表 | 不缩放 | 维持常量 |
3. **缩放公式**(次线性,内存通常随数据量次线性增长):
   ```
   mem_mb = clamp(base_mem × (size_mb / ref_mb) ** exp, min_mem, cap_mem)
   ```
   默认 `exp = 0.7`,`ref_mb = 5000`(即 5 GB 输入时等于现在的手调值,
   **现有"较佳值"原样成为基准点,不推翻已有调优**)。
   threads 缩放更保守(IO-bound 规则加核收益低):通常不缩放,
   只在 CPU-bound 规则(比对)允许小幅缩放并封顶。
4. **config schema 扩展**(`resources.<rule>` 新增可选 `scale` 段):
   ```json
   "star_align": {
     "threads": 12, "mem_mb": 32000,
     "scale": {"metric": "fastq_mb", "ref": 5000, "exp": 0.7,
               "min": 16000, "cap": 200000}
   }
   ```
   未配置 `scale` 的规则保持静态行为。
5. **试点范围**:先只给 top 大内存/长耗时规则配 scale
   (star_align、hisat2_align、cellranger、gatk 系、TEcount/TElocal、
   fragpipe、deepvariant 等),验证公式后再铺开。

### 方案 C:attempt 升级 + 自动重试(对不可预测数据兜底)

- `mem_mb` 走 callable:`base × 2 ** (attempt - 1)`(封顶 cap),
  run.py 加 `--restart-times 2`。
- 只对识别为 OOM 的失败升级重试(退出码 137 / 日志含 OOM 字样);
  其他错误照常报错。需要一个小的 wrapper 或 Snakemake 侧的失败分类。
- 优点:即使方案 B 的公式估偏,也不会一次失败;代价是重跑该规则。

### 方案 D:benchmark 自学习闭环(长期)

- 每条规则加 `benchmark:` 指令,Snakemake 原生记录
  `s / h:m:s / max_rss / max_vms` 等到 tsv。
- 定期跑一个校准脚本:按 (rule, 输入规模) 拟合实际峰值内存/耗时,
  输出"建议 resources 表"供人工审核后合入 config——而不是自动改。
- 这是方案 B 公式参数(exp/ref)的唯一可靠数据来源,也为 P5 提供事实依据。

### 对比

| 方案 | 实现成本 | 解决 | 风险/缺点 |
|---|---|---|---|
| A 主机感知 | 低(1-2 天) | P2 P3 P7 | 不解决数据量问题 |
| B 数据感知 | 中(需配 scale + 改关键规则) | P1(核心) | 公式参数初期靠估计,需要 D 校准 |
| C attempt 兜底 | 低-中 | P4 | 重跑有代价;需错误分类 |
| D benchmark 闭环 | 中-高 | P5,校准 B | 见效慢,属长期项 |

---

## 4. 推荐落地路线(分阶段)

1. **阶段 1 = 方案 A**:cores auto、多工作流内存分割、clamp 防死锁、
   统一 `rule_res` 到一处(消掉 KARRseq 本地拷贝)。不改任何规则声明,
   立刻解决"160 核只用 10 核"和"换主机卡死"。
2. **阶段 2 = 方案 B 试点**:helper 加 `size_mb`/`scale` 支持,
   给 8-10 个大内存规则配 scale,拿 2-3 个大小悬殊的数据集实测对比。
3. **阶段 3 = 方案 C**:`--restart-times` + OOM 分类 + attempt 翻倍。
4. **阶段 4 = 方案 D**:全规则 `benchmark:` 常态化,季度校准一次参数。

---

## 5. 涉及文件(实施时)

| 文件 | 改动 |
|---|---|
| `run.py` | `detect_host_profile()`、cores auto、多工作流预算分割、`--restart-times` |
| `modules/common/common.smk` | `rule_res` 扩展(scale/attempt/clamp) |
| `subworkflow/KARRseq.smk` | 本地 `rule_res` 拷贝改为复用统一实现 |
| `config/*.json` | `resources.<rule>.scale` 段(增量,非破坏) |
| `config/schema.schema.json` 及各 `*.schema.json` | 同步 scale 字段的 schema |
| 各模块 `.smk` | 试点规则传入 `size_mb=`(渐进,非一次性) |

---

## 6. 待审问题(请在审核时给出意见)

1. 数据规模度量:用输入文件**字节数**(便宜、即取即用)还是
   **reads 数**(准确,但要扫 fastq 或维护 meta)?建议字节数起步。
2. 缩放参数:`exp=0.7`、`ref=5 GB` 是否可接受?是否按规则类别
   (比对/聚合/绘图)分档取不同 exp?
3. 现有手调"较佳值":保留为基准点(本文档的假设)还是借机重新标定?
4. 多工作流并行的内存分割:均分,还是按各工作流输入总量加权?
5. 是否引入 Snakemake `--default-resources` 统一兜底(未声明资源的规则
   自动获得按主机推导的默认值)?
6. threads 是否也按数据量缩放?我的倾向是**基本不缩放**
   (IO-bound 规则加核收益低,且核心预算比内存更紧张)。
