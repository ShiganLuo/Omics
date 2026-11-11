---
name: modules编写规范
description: 撰写Snakemake模块时需要遵循此套规范
---

# 目的

统一模块的输入输出、参数命名和日志/环境配置方式，确保模块可被 subworkflow 复用。

# 适用范围

- 本文规范适用于 `workflow/Omics/modules/` 下的所有模块。
- 模块通过 subworkflow 的 `module` + `use rule` 引用，不直接读取 `run.py`。
- 模块只消费 subworkflow 传入的 `config`。

# 核心原则：最小可复用单元

## 输出决定流程走向

Snakemake 的核心机制是**反向推理**：从请求的输出文件出发，反向推导需要执行哪些规则。这是唯一正确的流程控制方式。

**三条铁律：**

1. **输出文件决定规则走向** — node.py 注册哪些输出文件 → Snakemake 反向推导需要跑哪些 rule。不被请求的输出对应的 rule 不会执行。
2. **输入是被动的** — input 函数只为已确定要执行的 rule 准备数据，不能决定 rule 是否执行。
3. **输入永远不可能决定流程走向** — 不要在 input 函数里写条件分支来"阻止" rule 执行，这是用输入控制流程，违反 Snakemake 设计原则。

**错误示例**（用输入控制流程）：
```python
# ❌ 错误：用 input 函数的返回值来阻止 rule 执行
def _get_gtf_for_heatmap(wildcards):
    if per_gene_mode:
        return []  # 想通过返回空来阻止 rule 执行
    return gtf
```

**正确示例**（用输出控制流程）：
```python
# node.py 中：per_gene 模式不注册 heatmap 输出，只注册 heatmap_gene 输出
if per_gene:
    for gene in gene_names:
        outfiles.append(f"{outdir}/heatmap/{sample}_{gene}_heatmap.png")
else:
    outfiles.append(f"{outdir}/heatmap/{sample}_genes_heatmap.png")
# Snakefile 中：input 函数只管准备数据，不做条件判断
def _get_gtf_for_heatmap(wildcards):
    return gtf  # 总是返回 GTF，流程由 node.py 输出控制
```

**判断标准：** 如果删除某个 input 函数的条件分支后，流程行为发生了变化，说明你正在用输入控制流程——这是错误的。

每个模块封装**一个工具**或**一个原子分析步骤**。规则只负责校验和调用，所有逻辑在 `bin/*.py` 中。

## 逻辑收敛到 bin/ 脚本（强制）

规则 run 块里出现下列任一情况，就是把支撑逻辑写错了地方，必须抽到 `bin/*.py`
引擎脚本、规则只留一条调用命令：

- **多阶段编排**：守卫分支（如空 BAM 占位输出）、前置数据加工（格式转换、
  索引、null 构建）、工具调用、清理，两段以上拼在同一个回执 `.sh` 里；
- **shell 控制流**：`if/then/else`、命令替换 `$(...)`、管道 `|` 出现在
  f.write 生成的脚本里（管道只允许作为两条列表构建命令之间的连接符）；
- **inline awk/sed/printf 魔法**：字段下标、转义符堆积的文本加工——这类
  逻辑在 Python 里是流式逐行 + 具名函数，可测试、可读。

标准形态：回执里只有 `python bin/<engine>.py --a x --b y` 一条命令加收尾
echo；引擎脚本自包含（规范 #0）、入口 `main()` 只编排、阶段函数单一职责、
声明型输出（含占位文件）在每条退出路径上都落盘。参照实现：
`modules/fibertools/bin/call_peaks.py` + `fibertools_analysis.smk::ft_call_peaks`。

## 模块不得引用外部脚本（强制）

模块**不得调用其他模块目录下的脚本**（`ROOT_DIR + "/modules/<B>/bin/y.py"`
这类跨模块路径出现在 cmd 里即违规）。可复用的分析能力是一个独立模块，
消费方通过 subworkflow 的 `module` + `use rule` 接线消费它的**输出文件**，
而不是 import/执行它的代码。理由：模块目录整体拷走即能运行是 modules
可迁移性的底线（同规范 #0 的脚本自包含），跨模块脚本引用在拷走后全部断裂。
上游模块只负责准备好输入文件（长表、矩阵、BED），下游模块以文件为接口。
参照实现：`fiberseq/mods` 准备 `{metric}_long.tsv` → subworkflow
`use rule sigtest_run from sigtest` 消费该表 → `fiberseq/mods` 的后处理规则
拆分显著对象；三个模块互不引用对方代码。

**模块 .smk 中禁止条件分支**（`if`/`else` 包裹 rule 定义）。模块应该定义所有可能的 rule，由 subworkflow 决定 `use rule` 哪些。条件逻辑（如是否启用 batch correction、是否跑 velocity）属于 subworkflow 层，不属于模块层。

判断标准：
- 该模块是否可以被独立复用？如果不能，说明它太大了。
- 该模块是否包含两个不同的可执行工具？如果是，拆成两个模块。
- 该模块的 conda 环境是否混入了不相关的依赖？如果是，说明职责不单一。

示例：
- `telomere` 模块只调用 `telogator2`，conda 环境只有 `telogator2`
- `centromere` 模块包含 `hifiasm` + `repeatmasker` 两步（组装→注释），但它们是同一个分析的上下游步骤，合为一个模块合理
- 将 telomere 和 centromere 合为一个模块是**错误的**：两者工具链完全不同

# 目录布局

## A. 简单 3 文件模块（最常见）

```
modules/<tool>/
  <tool>.smk
  <tool>.json
  <tool>.yaml
```

代表：samtools、fastqc、deepvariant、hiphase、pbsv、tabix、trgt

## B. 3 文件 + bin 目录（自定义脚本）

```
modules/<tool>/
  <tool>.smk
  <tool>.json
  <tool>.yaml
  bin/
    summarize.py
    plot.py
```

代表：StringTie、arriba、spectrum、RmrRNA、track、bowtie2、scanpy

`bin/` 存放该模块需要的 Python/R 辅助脚本，通过 `ROOT_DIR + "/modules/<tool>/bin/<script>.py"` 引用。

## C. 父目录 + 子模块目录

```
modules/<tool>/
  <tool>.smk                # 主规则或 prepare 规则
  <tool>.yaml               # 主模块 conda 环境
  <tool>.json
  <subtool>/
    <subtool>.smk           # 子工具规则
    [<subtool>.yaml]        # 子工具独立 conda 环境（可选；省略则共享父目录 yaml）
    <subtool>.json
```

代表：gatk（gatk_prepare.smk + gatk_bqsr/ + gatk_germline/ + gatk_somatic/ + gatk_RNAseq/）、samtools（samtools.smk + sort/）、openms（decoydatabase/ + searchengine/ + psmfdr/ + ...）

子模块默认引用父目录的 `<tool>.yaml`（`conda: "../<tool>.yaml"`）。当子模块需要不同环境时，在自己的目录下放 `<subtool>.yaml`，文件名 stem 必须唯一（见"环境命名规范"）。

## D. 公共模块

```
modules/common/
  common.smk
```

提供所有模块共享的工具函数（`setup_logger`、`sif` 等）。详见 [Common 模块](#common-模块)。

## 导入去重（强制）

common.smk 已 `import` 的名字（`sys / os / time / shutil / tempfile /
subprocess / shlex / re / snakemake.logging.logger`）随 `include` 一并进入
模块命名空间——模块 .smk **不得重复导入**它们，也不得在 include 之前预先
导入。新写规则前先查 common.smk 的 import 区：只有它没有的包（numpy、
pandas、json、gzip …）才在模块里 import。重复导入不报错、不省时间，
纯属噪音（102 个文件 186 行存量曾因此堆积）。

# 配置入口

模块通过 `config.get()` 读取参数。subworkflow 组装 config dict 传入模块，字段约定：

| 字段 | 用途 | 使用场景 |
|---|---|---|
| `indir` | 输入目录 | 所有模块 |
| `outdir` | 输出目录 | 所有模块 |
| `logdir` | 日志目录 | 所有模块 |
| `samples` | 样本 ID 列表 | 大多数模块 |
| `paired_samples` / `single_samples` | 双端/单端样本列表 | 比对模块 |
| `sample_groups` | 样本分组字典 | StringTie |
| `group_pairs` | 对比组定义 | DESeq2 |
| `bam_dir` / `vcf_dir` | 独立 BAM/VCF 目录 | hiphase |
| `Procedure` | 工具可执行路径字典 | 所有模块 |
| `Params` | 工具参数字典 | 需要额外参数的模块 |
| `genome` | 参考文件路径字典 | 需要参考基因组的模块 |
| `ROOT_DIR` | 项目根目录，用于引用 bin/ 脚本 | 有自定义脚本的模块 |
| `env` | conda/SIF 环境映射 | 所有模块（通过 common.smk） |

## per-sample 输入的两种形态与选择（2026-11 sexcheck 定案）

模块消费 per-sample 文件（BAM/VCF/表）有两种合法形态，按**模块的复用面**选择：

**形态 A：`indir` + 布局参数（单一布局模块）** —— `indir`（上游 config 的 `outdir`）+ `bam_suffix` 等布局参数 + 样本列表，模块 input 函数自己拼 `{indir}/{sample_id}/{sample_id}{suffix}`（参照 `msi.smk::_bam`、`cnvkit.smk`、`fragment_size.smk`）。适用于模块只接一种上游布局的场景。

**形态 B：`sample_bams` 字典（跨布局公用模块）** —— config 传「样本 ID → 完整文件路径」字典（参照 `sexcheck.smk`、`hiphase` 的 `bam_dir` 键）。适用于一个模块被**布局互不相同**的多个流程复用（markdup/.dedup/.sorted/STAR 直出/polygenomes 多层）——布局差异本来就该 subworkflow 知道，字典把模块从布局知识里解耦。

**形态 B 的三条强制约束**（否则退回形态 A）：

1. **字典值必须引用上游 config 的 `outdir` 拼**（f-string 引用 `xxx_config["outdir"]`），禁止手拼裸路径副本——路径知识仍只写一次，上游改布局字典自动跟随（同 subworkflow.md 衔接规范）；
2. **DAG 依赖验证**：字典值路径必须能被 Snakemake 反推出生产者规则。dry-run 中消费规则的 reason 应出现 "Input files updated by another job: <字典指向的 BAM>"；若字典指向**未被任何规则声明为 output** 的文件（副产物/外部文件），依赖会静默缺失——这类输入不能只用字典，需另行声明 input 锚点；
3. **样本 ID 必须是纯 sample_id**：per-genome 布局把 genome 放进值的路径段（`f"{cfg['outdir']}/{g}/{s}/..."`），**禁止**把 `genome/sample` 斜杠串当字典键——键污染 wildcards 和日志目录。

节点侧配套：形态 B 同样要在 node.py 注册模块产出（outfiles 驱动 DAG）——注册缺失时模块静默不跑，与写法无关（2026-11 实录：node.py 注册被误 checkout 回退未察觉，sexcheck 消失于 DAG，排查时一度误判为字典写法破坏依赖；dry-run reason 证实字典路径的依赖追踪完好）。

判断标准：模块 config 出现 sample→路径字典时检查上述三条；违反任何一条 = 改回形态 A。

# 脚本撰写规范（任意语言）

模块 `bin/` 下的脚本、辅助工具、交付给下游的任何代码一律遵守，不限语言：

0. **脚本自包含（强制）** — `modules/` 内的脚本不得引用模块目录之外的任何代码：禁止 `sys.path.insert` 引向仓库根、禁止 `from src.common... import`、禁止跨模块 `import`（`modules/<A>/bin/x.py` 导入 `modules/<B>/bin/y.py` 同样禁止）。模块目录整体拷走即能运行，这是 modules 可迁移性的底线。需要 `src/common/util/` 里的工具（如 `setup_logger`）时，把该函数（连同其依赖）直接复制进本模块 `bin/` 下的脚本或私有 helper 文件，保持行为与格式一致即可；`.smk` 规则层不受此条约束（`include: "../common/common.smk"` 由 subworkflow 的 `ROOT_DIR` config 保证，且 Snakemake 模块机制本身支持相对 include）。参照实现：`modules/sigtest/bin/sigtest.py`（直接复制了 LogUtil 的 `setup_logger`）。存量模块中引用 `src.common` 的脚本（fiberseq、DESeq2、arriba 等）在下次触碰时收敛为本规范，不要求一次性回改。

1. **功能注释全覆盖** — 每个函数/方法必须有说明职责的 docstring/头注释：做什么、参数含义与单位、返回什么、何时返回空/报错；关键步骤（算法选择、外部工具调用、边界处理、为什么这样做）配行内注释。禁止复述代码的无信息注释（`# 循环`）。
2. **能抽函数就抽函数** — 独立步骤、重复逻辑、超过 ~30 行的代码块提取为单一职责函数，函数名自解释（动词开头）；入口（`main` / `__main__`）只做编排：解析参数 → 按序调用阶段函数 → 收尾日志。
3. **支持类型的语言必须标注形参与返回类型** — Python 用 PEP 484 注解（`def f(x: str) -> int`，复杂结构先定义类型别名，如 `Rec = Dict[str, Any]`）；R/Rust/Go/Java/TS 等用各自类型系统。类型是文档的一部分，也让静态检查能抓出真实的不精确。
4. **优雅优先** — 命名自解释；控制流扁平（早返回/卫语句）；同一数据的读取与规范化只写一次（单一来源）；删除恒真恒假的死代码；不做多余包装（如单字符串包列表再 `[0]` 解包）；错误信息带修复指引；同一函数的返回值语义一致（不混返 bool/str/None）。
5. **能向量化不用循环** — 表格/数组运算用 pandas/numpy 向量化（布尔掩码、`groupby().agg()`、`np.where`、`.str` 访问器），不写逐行 Python 循环；纯循环只留给无法向量化的控制流。
6. **能并行不串行** — 相互独立的工作单元（多样本、多文件、分块）用多线程（I/O 密集）或多进程（CPU 密集）；注意库级锁陷阱（如 HDF5 全局锁会把多线程串行化，实测更慢，须换多进程）；并行度从 `--threads` 参数传入，合并结果后排序，保证输出确定性。琐碎工作不要为并行而并行。
7. **大文件流式读取** — VCF/GTF/BED/FASTA 等可能上 GB 的文件逐行或分块读取（生成器、`pandas.read_csv(chunksize=...)`），不整读进内存；确需全量驻留时在注释中说明内存预算。
8. **长耗时脚本必须有缓存机制** — 输入不变时复用已有产物，不重算：输入指纹（内容哈希或路径+大小+mtime、全部影响输出的参数）写入 stamp 文件，命中即跳过；缓存粒度对准最贵的步骤（外部工具调用），输入变了要能精确失效；提供强制重算开关（如 `--force-rerun`），默认使用缓存；开关暴露到流程配置 JSON 中。

> 向量化与流式的取舍：向量化要求数据在内存——中小表格直接向量化；超大文件先流式过滤/聚合出小结果集，再对结果向量化。两者冲突时以"不整读大文件"优先。

## 输出声明:禁用 `directory()`,用标记文件或文件清单(强制)

规则的 `output:` 不要写 `directory()`。Snakemake 在每次 job 执行前会**强制
删除** directory() 输出的整个目录(`IOFile.remove` 对 `is_directory` 走
`remove_non_empty_dir=True`,非空目录保护对它无效),脚本级缓存、手工管理
的中间产物会被连锅端;且 directory() 只跟踪"目录存在",内容变化不触发重跑。
正确做法:

1. **标记文件**(推荐):产物目录不声明为输出,规则成功后 touch 一个标记
   (如 `outdir/{comparison}/mei/mei.done`)作为 DAG 跟踪的输出;
2. 或声明**具体文件清单**:真正要交付的那几个文件逐个列出。

长耗时脚本带缓存时尤其依赖此条——缓存放在 directory() 输出里,每次重跑
都会被清空,永远命不中。

# 命名约定

- 规则名描述动作：`*_index` / `*_align` / `*_run` / `*_phase`
- 输出目录按模块名组织：`{outdir}/<module>/{sample_id}/...`
- 收尾规则：`<tool>_result`，仅声明 `input:` 作为依赖聚合点（无 `shell`/`run`）
- 报告规则：`<tool>_report`，聚合所有样本输出生成跨样本报告

# 日志规范

规则日志（`log:` 指令）按**规则粒度**分两种路径，模块从 config 读取 `logdir`：

- **单样本规则**（一次处理一个样本，log 路径含 `{sample_id}` 通配符）：

```python
log:
    logdir + "/{sample_id}/{规则名}.log"
```

- **多样本规则**（一次任务处理多个样本 / 组级聚合，log 路径无 `{sample_id}`）：

```python
log:
    logdir + "/{module名}/{规则名}.log"
```

subworkflow 会把 `logdir` 指到对应分段目录：单样本模块收到 `{总logdir}/sample`，
多样本模块收到 `{总logdir}/group`（分流规则见 `subworkflow/subworkflow.md` 的日志规范）。

## 同时包含单样本与多样本规则的模块

必须读取**两个变量**，subworkflow 会同时传入：

```python
logdir = config.get("logdir", "log")
logdir_combine = config.get("logdir_combine", logdir)
```

- 单样本规则用 `logdir + "/{sample_id}/{规则名}.log"`
- 多样本规则一律用 `logdir_combine + "/{module名}/{规则名}.log"`

参考实现：`modules/trgt/trgt.smk`（纯单样本）、`modules/fragment_size/fragment_size.smk`
（混合模块：`samtools_stats` 单样本、`getFragmentSize`/`plotFragmentSize` 走 `logdir_combine`）。

# 规则编写风格：`run:` 块 + shell 脚本生成

所有规则统一使用 `run:` 块，不允许使用纯 `shell:` 块。命令通过列表构建，写入 `.sh` 脚本后执行。

## 标准模板

```python
rule <tool>_<action>:
    input:
        ...
    output:
        ...
    log:
        logdir + "/{sample_id}/<tool>_<action>.log"
    threads: N
    conda:
        "<tool>.yaml"
    container:
        sif("<tool>.yaml")
    params:
        <tool> = config.get("Procedure", {}).get("<tool>") or "<tool>",
        ...
    run:
        log_path = str(log)
        open(log_path, "w").close()
        rule_logger = setup_logger("<tool>_<action>", log_file=log_path)
        try:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start <tool> <action> for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.xxx))
            os.makedirs(sample_outdir, exist_ok=True)
            os.makedirs(os.path.join(logdir, "cmd", str(wildcards.sample_id)), exist_ok=True)
            # 带时间戳的备份(log/cmd 历史);outdir 只放无时间戳的当前命令(见"命令回执双写")
            script = os.path.join(logdir, "cmd", str(wildcards.sample_id), f"<tool>_<action>_{wildcards.sample_id}_{current_time}.sh")
            plain_script = os.path.join(sample_outdir, f"<tool>_<action>_{wildcards.sample_id}.sh")
            cmd = [
                params.<tool>, "<action>",
                ...
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "<tool> <action> for {wildcards.sample_id} was completed successfully"\n')
            rule_logger.info("Executing: " + " ".join(cmd))
            shutil.copyfile(script, plain_script)
            shell(f"bash {plain_script} >> {log_path} 2>&1")
        except Exception as e:
            rule_logger.error(f"Error occurred during <tool> <action> for sample {wildcards.sample_id}: {e}\n")
            raise RuntimeError(f"Error occurred during <tool> <action> for sample {wildcards.sample_id}: {e}, you can check the log file {log_path} for more details.\n")
```

每个module snakefile必须含有一个面板规则记录module主要输出，便于node.py确定输出文件，格式如下：

```python
rule <module>_result:
    input: ...
```
## 命令回执双写（强制）

每条规则生成的 `.sh` 回执**必须同时写到两个位置**：

| 位置 | 文件名 | 行为 |
|---|---|---|
| 输出目录（`sample_outdir`） | `<rule>_<sample>.sh`（**无时间戳**） | 每次运行覆盖，始终对应当前命令，与结果同在，可直接复制重跑 |
| `<OUTDIR_ROOT>/log/cmd/`（**扁平，无子目录**） | `<rule>_<sample>_<时间戳>.sh` | 备份池：按运行累积，文件名自带规则+样本+时间戳，天然唯一 |

`OUTDIR_ROOT`（大写常量，common.smk 提供）= 流程总输出目录，由 subworkflow
统一传入（`"OUTDIR_ROOT": outdir`），缺省回退 `ROOT_DIR`。配套常量
`RECEIPT_DIR = <OUTDIR_ROOT>/log/cmd`，common.smk 已 `makedirs`，规则只拼
文件路径：

```python
# ✅ 正确
script = os.path.join(RECEIPT_DIR, f"<rule>_{wildcards.sample_id}_{current_time}.sh")
plain_script = os.path.join(sample_outdir, f"<rule>_{wildcards.sample_id}.sh")
...
shutil.copyfile(script, plain_script)
shell(f"bash {plain_script} >> {log_path} 2>&1")

# ❌ 错误
script = os.path.join(logdir, "cmd", str(wildcards.sample_id), f"...")  # 挂在模块 logdir 下,不收敛
script = os.path.join(RECEIPT_DIR, f"{wildcards.sample_id}/xxx.sh")    # 备份池是扁平的,禁子目录
os.makedirs(os.path.join(logdir, "cmd", ...))                          # 目录已由 common 保证
```

写文件顺序：先写时间戳备份版（`script`），`shutil.copyfile` 复制为无时间戳
版（`plain_script`），然后 `shell()` 执行 `plain_script`。注意 plain 文件名若
含通配符插值（`{wildcards.sample_id}`）必须是 f-string。

## 输出目录卫生（强制）

- **中间大文件不留输出目录**：规则不得在自己的 outdir 里复制/重排上游大文件（BAM/大 BED）。需要排序/索引版上游文件时，让上游规则产它或直接引用上游路径
- **未压缩大 BED 一律 bgzip**：BED 超过 1G 的产出用 bgzip+tabix（容器选含 htslib 的环境），下游用流式读取（`open_maybe_gzip`）
- **`.tmp` 中间物用前清理**：`path.unlink(missing_ok=True)`，失败中断的残留下次运行前清掉
- **转换成功即删中间物**：如 bedGraph→bigWig 成功后删除 bedGraph

## 要点

1. **`open(log_path, "w").close()`** — 清空旧日志，避免追加混淆
2. **`setup_logger`** — 从 common.smk 导入，统一日志格式
3. **`rule_logger.info(...)`** — 开始/结束标记，便于定位问题
4. **`current_time` 时间戳** — 脚本名加时间戳避免并发写冲突
5. **`cmd` 列表构建** — 参数逐项添加，条件参数用 `if` 追加
6. **`" ".join(cmd)` vs `" ".join(shlex.quote(str(x)) for x in cmd)`** — 不用 `shlex.join()`（需 Python 3.8+，scTE/mimseq 等模块 conda env Python < 3.8 不兼容）
   - **`" ".join(cmd)`** — 参数不含空格/特殊字符时使用（工具路径、flag、数字、标准路径）
   - **`" ".join(shlex.quote(str(x)) for x in cmd)`** — 参数可能含空格或特殊字符时使用（用户输入、含空格路径、动态文件名），兼容 Python 3.3+
7. **`#!/bin/bash` + `set -euo pipefail`** — 脚本内任何命令失败立即退出，不会执行后续的 echo
8. **`echo` 完成标记** — 脚本末尾写 `echo "...completed successfully"`，执行成功时日志中有明确结尾
9. **`shell(f"bash {script} >> {log_path} 2>&1")`** — 直接路径，不加 `shlex.quote` 包装
10. **`logger.error(...)` + `raise e`** — snakemake logger 记录错误后 re-raise
11. **`import shutil`** — 模块头部需要（回执双写的 `shutil.copyfile`）

## 线程限流环境变量

Snakemake 的 `threads:` 只约束**调度器记账**，不会约束进程内数值库的线程数：
numpy/BLAS、numba 等默认开满全部物理核。实测一个声明 4 线程的 scanpy_auto
进程开出 645 个 OS 线程，多个 job 并行时互相抢占、白白增加切换开销。因此
生成的 shell 脚本必须在命令**之前**导出以下变量，把线程预算传给数值库：

| 变量 | 约束对象 |
|---|---|
| `OMP_NUM_THREADS` | OpenMP 并行（C/C++ 数值库、scipy 部分例程） |
| `OPENBLAS_NUM_THREADS` | OpenBLAS（numpy 矩阵乘法、线性代数） |
| `MKL_NUM_THREADS` | Intel MKL（numpy 链接 MKL 后端时生效） |
| `NUMBA_NUM_THREADS` | numba JIT（scanpy/pynndescent 的 kNN、UMAP） |
| `VECLIB_MAXIMUM_THREADS` | macOS Accelerate/vecLib（Linux 导出无害，保跨平台一致） |

写法（`set -euo pipefail` 之后、命令之前）：

```python
with open(script, "w") as f:
    f.write("#!/bin/bash\nset -euo pipefail\n")
    f.write(f"export OMP_NUM_THREADS={threads}\n")
    f.write(f"export OPENBLAS_NUM_THREADS={threads}\n")
    f.write(f"export MKL_NUM_THREADS={threads}\n")
    f.write(f"export NUMBA_NUM_THREADS={threads}\n")
    f.write(f"export VECLIB_MAXIMUM_THREADS={threads}\n")
    f.write(" ".join(cmd) + "\n")
```

要点：

1. **值取 run 块的 `threads` 变量，禁止写死。** 链路：workflow config JSON 的
   `resources.<rule>.threads` → `rule_res()` → 规则 `threads:` → run 块的
   `threads` 变量。JSON 是资源数量的唯一入口，改线程数只改 config；且
   `threads` 拿到的是调度器**实际分配**的值（声明超出可用核时被自动缩放），
   限流跟着调度走，两头永远一致。
2. 带 `--threads`/`--localcores` 类参数的 CLI 工具（cellranger、STAR 等）应
   把 `threads` 传给工具自身参数；环境变量是给**不接受线程参数的数值库**
   （BLAS/numba）兜底的，两者可以并存。
3. 线程数取值参考真实并行度：numba kNN/harmony 类负载（scanpy_auto、
   scanpy_advanced）12~16 之后收益趋平，申报更多只增加调度压力。

## 带条件参数

```python
cmd = [params.tool, "--required-flag", input.file]
if params.optional_param:
    cmd += ["--optional", params.optional_param]
```

注意：nargs+ 式参数必须逐项拆分，不能 `["-f", "heatmap volcano pca"]`，而要 `"-f", "heatmap", "volcano", "pca"`。

## 多命令序列

```python
with open(script, "w") as f:
    f.write("#!/bin/bash\nset -euo pipefail\n")
    f.write(" ".join(cmd1) + "\n")
    f.write(" ".join(cmd2) + "\n")
    f.write(f'echo "<tool> completed at {current_time}"\n')
```

## Python 数据处理 + shell 命令混合

当需要先用 Python 准备文件再调用外部工具时，Python 逻辑直接写在 `run:` 块中：

```python
run:
    log_path = str(log)
    # Python 准备
    with open(map_file, "w") as f:
        for k, v in mapping.items():
            f.write(f"{k}\t{v}\n")
    # shell 执行
    with open(script, "w") as f:
        f.write("#!/bin/bash\nset -euo pipefail\n")
        f.write(" ".join(cmd) + "\n")
        f.write(f'echo "<tool> completed at {current_time}"\n')
    shell(f"bash {script} >> {log_path} 2>&1")
```

## 为什么不使用纯 `shell:` 块

- 无法记录结构化日志（时间戳、样本名、开始/结束）
- 无法在失败时追加错误信息到日志
- 无法动态构建参数列表

# 动态输入函数模式

比对类模块需要根据样本类型选择不同输入文件：

```python
def get_alignment_input(wildcards):
    """Dynamically determine paired-end or single-end input."""
    paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
    single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

    if wildcards.sample_id in paired_samples:
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        return [single]
    else:
        raise ValueError(f"Sample {wildcards.sample_id} not in paired_samples or single_samples")
```

代表：fastqc、star、bowtie2、hisat2

索引回退模式（优先用配置中的已有索引，否则用模块生成的）：

```python
def get_hisat2_index(wildcards):
    config_prefix = config.get('genome', {}).get('index_prefix')
    if config_prefix:
        first_file = f"{config_prefix}.1.ht2"
        if os.path.exists(first_file):
            return [f"{config_prefix}.{idx}.ht2" for idx in [1,2,3,4,5,6,7,8]]
    return [outdir + f"/index/genome.{idx}.ht2" for idx in [1,2,3,4,5,6,7,8]]
```

代表：bowtie2、hisat2、star

# 输入文件校验与 run 块规范（强制）

## 消费的文件必须声明为 input

规则消费的一切文件（参考序列、GTF、索引、TE 数据库、上游输出）都必须出现在 `input:` 中。Snakemake 在 **DAG 构建期**校验 input 存在性：缺文件在任何 job 启动前抛 `MissingInputException`。藏在 `params:` 或 run 块里的路径则要到**执行期**才报错——dry-run 假通过，job 跑到一半死于缺文件这类非代码原因。

❌ 错误（文件不进 input，执行期才炸）：
```python
params:
    gtf = config.get("genome", {}).get("gtf"),
run:
    if params.gtf:
        cmd += ["--gtf", params.gtf]   # 路径写错/文件缺失，job 跑到一半才发现
```

✅ 正确（一个输入函数声明全部文件，rule 里 `unpack` 展开）：
```python
def sv_mei_inputs(wildcards):
    """声明全部输入文件；缺文件在建图期即失败。"""
    files = {
        "merged_vcf": outdir + f"/{wildcards.comparison}/{wildcards.comparison}_merged.vcf",
        "classifier": classifier,
    }
    for name, key in (("fasta", "fasta"), ("gtf", "gtf")):
        if genome.get(key):
            files[name] = genome[key]   # 可选文件：配置了才进 input
    return files

rule sv_mei:
    input:
        unpack(sv_mei_inputs),
```

## 可选文件的两种合规写法

1. **推荐**：配置了才加入输入函数返回的 dict，run 里 `getattr(input, "gtf", "")` 安全取值；
2. `[x] if x else []` 列表形式，run 里 `if input.x:` 门控。

必填文件永远返回裸字符串；把单个字符串包成单元素列表再 `[0]` 解包是噪音，禁止。静态参考文件（genome fasta/GTF/数据库）也不例外——"它一直都在"不是免检理由。

## run 块禁止模块级全局变量

run 块只允许引用 `input / params / wildcards / output / log`；配置值经 `params` 传入。同一配置块的规范化读取（str→list 归一、默认值）收敛为一个函数，输入函数与 `params` 共用同一来源，禁止重复深链取值：

```python
def _sv_mei_cfg():
    """sv.mei 配置块规范化（dfam_h5 恒为 list）。"""
    mei = dict(config.get("Params", {}).get("sv", {}).get("mei", {}) or {})
    dfam_h5 = mei.get("dfam_h5") or []
    mei["dfam_h5"] = [dfam_h5] if isinstance(dfam_h5, str) else list(dfam_h5)
    return mei

rule sv_mei:
    params:
        mei = _sv_mei_cfg(),   # run 里用 params.mei.get("min_len", 100)
```

## Params 不得无故提升为模块全局变量（强制）

配置值（`Params.*` 及其派生量）**在需要它的位置读取**：input 函数与
`params:` 表达式内直读 `config.get(...)`，或经上面的规范化函数统一取。
只有一种情况允许定义模块级变量：`input:`/`output:`/`log:`/`threads:`
等 **parse 期指令**自身需要的派生值（输出文件名、轴开关、expand 列表）
——指令里无法调用 params，只能用全局可见的量，且定义处注明是为指令服务。

判断标准：如果一个全局变量只被 `params:` 或 run 块引用，放全局**没有任何
收益**（不省解析时间），只损失可读性——删掉，在使用点读取。

反例（禁止）：
```python
_min_fibers = int(p.get("min_fibers", 3))   # 只被 params 引用，却定义在全局
_tmp_dir = config.get("tmp_dir") or ...     # 同上
rule x:
    params:
        min_fibers = _min_fibers,
```
正例：
```python
def _mods_cfg():
    return config.get("Params", {}).get("fiberseq_mods", {}) or {}

_long_name = ...        # output 指令需要，允许（注明用途）
rule x:
    params:
        min_fibers = int(_mods_cfg().get("min_fibers", 3)),
```

## 必填配置缺失 = 建图期报错

输入函数里直接 `raise ValueError("...set Params.x.y...")`，让配置错误在建图期暴露并给出修复指引。物种、基因组等科学参数禁止静默猜测（never guess a species）。

## bin/*.py 入口配套

入口脚本 `main()` 只做编排（解析参数 → 各阶段函数按序调用），阶段逻辑拆单一职责函数；外部工具调用校验返回码与产物存在性，失败信息带修复指引。

# 收尾规则（_result）

每个模块提供 `<tool>_result` 规则，仅声明 `input:` 作为依赖聚合点，供 subworkflow 的 `use rule ... as ...` 引用：

```python
rule deepvariant_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
```

```python
rule scanpy_result:
    input:
        h5ad = expand(auto_dir + "/{t}_auto.h5ad", t=tissues)
```

# 报告模块（_report）

报告模块消费上游所有模块的输出（TSV/CSV/PNG/BAM/log），生成 PPTX + XLSX 报告。不重复计算上游分析，只读取已有结果文件。

代表：RNAseq_report、ncRNAseq_report、PeakCalling_report

## 目录结构

```
modules/<Workflow>_report/
  <Workflow>_report.smk        # generate_report + report_result 两个 rule
  <Workflow>_report.json
  <Workflow>_report.yaml       # python-pptx, matplotlib, pandas, openpyxl
  bin/
    generate_report.py
```

## 两个 rule

```python
REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "<Workflow>_report", "bin", "generate_report.py")

rule generate_report:
    input:
        # 用 expand() 声明所有上游依赖文件
        ...
    output:
        report = outdir + "/<Workflow>_report.pptx",
        file_inventory = outdir + "/<Workflow>_report_files.xlsx",
    log:
        logdir + "/<Workflow>_report.log"
    threads: 1
    conda:
        "<Workflow>_report.yaml"
    container:
        sif("<Workflow>_report.yaml")
    params:
        samples = samples,
        title = config.get("Params", {}).get("report", {}).get("title") or "Default Title",
        lang = config.get("Params", {}).get("report", {}).get("lang") or "zh",
        img_dir = outdir + "/ppt_results",
        script = REPORT_SCRIPT,
    run:
        # 标准 run 块模板
        ...

rule report_result:
    input:
        report = outdir + "/<Workflow>_report.pptx",
        file_inventory = outdir + "/<Workflow>_report_files.xlsx",
```

## subworkflow 集成

```python
# subworkflow/<Workflow>.smk 末尾
<Workflow>_report_config = {
    "ROOT_DIR": ROOT_DIR,
    "env": config.get("env", {}),
    "outdir": outdir,
    "logdir": logdir,
    "samples": all_samples,
    "Params": {
        "report": config.get("Params", {}).get("report", {}),
    },
}
module <Workflow>_report:
    snakefile: "../modules/<Workflow>_report/<Workflow>_report.smk"
    config: <Workflow>_report_config
use rule generate_report from <Workflow>_report as <Workflow>_generate_report
use rule report_result from <Workflow>_report as <Workflow>_report_result
```

# Subworkflow 跨模块衔接（强制）

subworkflow 中后一个模块消费前一个模块的输出时，**输入目录一律引用上游
config 的 `outdir`**，不得手拼同一路径的 f-string 副本：

```python
# ✅ 正确：路径知识只写一次，wiring 跟随 DAG
fiberseq_footprint_config = {
    "indir": elements_dir,                        # = extract 步的 outdir
    "motif_dir": fiberseq_motif_config["outdir"], # 上游模块 config["outdir"]
    "outdir": f"{analysis_outdir}/footprint",
    ...
}

# ❌ 错误：同一路径出现两个手拼副本，改一处漏一处
"motif_dir": f"{outdir}/fiberseq_analysis/motif",   # 与 motif 模块的 outdir 重复定义
```

规则：
- 上下游目录在 subworkflow 里定义**模块级常量**（`elements_dir`、
  `phased_dir`），由产出方 config 引用同一常量；Snakefile 是顺序执行的，
  **禁止前向引用**后面才定义的 config 字典（NameError）
- 模块 .smk 内部接收这类输入用独立 config 键（`motif_dir`/`quant_dir`/
  `annot_dir`/`hap_dir`/`vcf_dir`），fallback 用 `config.get("xxx") or indir`
- 判断标准：一个物理路径在整个 subworkflow 里只出现一次字面定义

# 模块测试（run_module.py --test）

模块没有固定输入（不走 meta.tsv），但必须可独立加载、可建 DAG。测试框架
`run.py --test all` 会自动追加模块阶段（内部调 `run_module.py --test all`），
对每个 `modules/**/<stem>.smk + <stem>.json` 单元检查：

- config 加载 + schema 校验（"Required field is empty" 视为待填模板，不判失败）
- `default_targets()` 目标展开 + wrapper 生成 + snakemake 干跑
- 测试路径自动注入（复用 `SchemaValidator.generate_test_paths`，与 workflow
  测试同一机制）；输入缺失/InputFunction/未绑定通配符/配置门控的 DAG 停止
  算通过（无固定输入的预期），语法/接线/配置名错误判失败

单独跑：`python run_module.py --test all`（或 `--test <module>`）。

# Conda 环境模板（.yaml）

```yaml
name: <tool>
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - <tool>>=<version>
```

> channel 顺序必须是 `conda-forge` -> `bioconda` -> `defaults`。bioconda 依赖 conda-forge 的包，顺序颠倒会导致依赖解析失败。

## 环境命名规范（强制）

YAML 文件名、`name:` 字段、SIF 文件名必须三者一致，等于模块目录名：

| 文件 | 命名 |
|------|------|
| YAML 文件名 | `<module_dir>.yaml` |
| YAML `name:` 字段 | `<module_dir>` |
| SIF 文件名 | `<module_dir>.sif` |

`common.smk` 中的 `sif()` 函数用 YAML 文件名 stem 查找 SIF 路径，不读取 YAML 内容。命名不一致会导致容器路径解析失败。

> **例外**：如确实需要名称与目录不同，在 config 的 `env` 字典中显式映射：`"yaml_stem": "/path/to/actual.sif"`。

# Common 模块

`modules/common/common.smk` 提供所有模块共享的工具函数。

## 包含内容

```python
import sys, os, time, shutil
from snakemake.logging import logger
import shlex

ROOT_DIR = config.get("ROOT_DIR", ".")

# 自动将 ROOT_DIR/src 加入 sys.path
from common.util.LogUtil import setup_logger

# SIF 容器路径解析
def sif(yaml_filename: str) -> str:
    """用 YAML 文件名 stem 查找 SIF 路径。
    1. config["env"][<stem>] -- 精确映射
    2. config["env"]["env_dir"] + <module_dir>/<stem>.sif -- 约定回退
    """
```

## 使用方式

```python
# modules/<tool>/<tool>.smk
include: "../common/common.smk"

# 子目录模块
include: "../../common/common.smk"

# subworkflow
include: "../modules/common/common.smk"
```

# Params传递规范

对应模块参数传递最好先在头部定义模块参数config,这样以后修改参数位置整体移动时可以方便修改，比如

```python
decoil_params = config.get("Params", {}).get("decoil", {})
```

- **参数链条保持完整**：配置字段从模板到消费方保持同一条完整嵌套路径（如 `Params.sv.mei.te_buckets`），subworkflow 组装模块 config 时保留层级，不缩短、不改根、不拍平。
- **只传用得到的参数**：subworkflow 组装模块 config 的 `Params` 时只放该模块实际读取的子树/字段，用不到的参数不传。

# Pitfalls

## 1. 子目录 conda 路径

引用父目录的共享 yaml 时用 `conda: "../<parent>.yaml"`；引用子目录自己的 yaml 时用 `conda: "<subtool>.yaml"`。不能写 `conda: "../modules/<tool>/<tool>.yaml"`。

## 2. run 块中的 shell 脚本路径

使用 `current_time` 时间戳避免并行冲突：
```python
current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
script = f"{outdir}/{wildcards.sample_id}/tool_{current_time}.sh"
```

## 3. outfiles 路径一致性

run.py 中 `outfiles` 的路径必须与规则的 `output` 完全匹配，否则 Snakemake 静默跳过或全量重建。

## 4. expand() 中的通配符转义

在 `expand()` 内引用外部通配符时用双花括号 `{{}}`：
```python
expand(outdir + "/{case}/{case}.{{genome_version}}.msisensor", case=samples.index)
```

## 5. bin/*.py 使用 logging 而非 print

所有 bin/*.py 脚本必须使用 Python `logging` 模块，不能用 `print()` 输出信息。在 `main()` 中调用 `setup_logging()` 初始化。

## 6. nargs+ 参数拆分

命令列表中 nargs+ 参数必须逐项拆分：`"-f", "heatmap", "volcano", "pca"`，不能 `"-f", "heatmap volcano pca"`（shlex 会认为是一个字符串）。

## 7. output 不能是函数（只有 input 可以）

Snakemake 只允许 input 指定为函数（`Only input files can be specified as
functions`），output 必须是 parse 期可确定的 str 或 str 列表。配置驱动的
可变输出（如按 `Params.*.formats` 出多种图）用 parse 期列表推导展开：

```python
output:
    plot = [outdir + "/{sample_id}/viz/struct." + fmt
            for fmt in (p.get("formats") or ["png"])],
```

用字符串拼接保留 `{sample_id}` 给通配符解析（f-string 要写成 `{{sample_id}}`）。
同一份 config 列表放 `params` 回读给 run 块，保证脚本产出与声明一一对应。
运行时才知道产物的场景用 `directory()` / marker / checkpoint，不要把可能
产物全声明成 output 靠"能产哪个产哪个"（MissingOutput 不允许部分满足）。

## 8. directory() 输出会在规则运行前被清空

声明 `directory(path)` 为 output 的规则，Snakemake 在运行前删除整个目录
（含其中"别人写的"文件）。多个规则向同一目录写文件时，后跑的会清掉先跑
的产出（sv 模块 `enrichment/` 被 `sv_exp_enrichment` 清空、exp_specific
的 RepeatMasker 结果被误删即是实例）。规矩：一个 directory 输出只属于一条
规则；其他规则要写就写自己的目录。

## 9. outdir 里的时间戳回执会无限堆积

带时间戳的 `.sh` 每次重跑多一份，结果目录因此逐渐变脏。回执一律双写：
outdir 放无时间戳版（覆盖），`log/cmd/` 放时间戳版（历史），见
「命令回执双写」。

## 10. 扁平输出目录下 plain_script 不带样本名会并发互踩

多样本规则的输出若都在同一目录（如 `outdir/idxstats/{sample}.tsv`，dirname
相同），`plain_script` 写成固定名 `rule.sh` 时并行 job 互相覆盖，A 样本会执行
到 B 样本的命令（症状：产出串样本 / MissingOutput，日志却显示成功）。规矩：
`plain_script = <rule>_{wildcards.<key>}.sh`，key 必须含区分并行 job 的通配符
（样本/对比/基因组）；out-dir 含通配符目录（`outdir/{sample_id}/`）的规则不受
影响。实例：sexcheck_idxstats 曾因固定名在 6 并行下 5 个 job 串台。

# QC 指标接入规范（新增 QC 产物时，强制）

任何模块产出 QC 类指标（metrics/统计表/日志推断值）并希望进入
`qc_summary.tsv` 时，必须一次完成以下四件事，缺一即视为接入未完成：

1. **解析器**：在 `modules/qc_summary/bin/qc_summary.py` 新增 `parse_<stage>`
   函数（输入 stage_dir + sample，返回扁平 `{metric: value}` dict，键小写、
   值保持字符串；文件缺失返回 `{}`），并注册进 `PARSERS`。选入指标走模块级
   白名单常量（如 `_HS_METRICS_FIELDS`），不要整行收录无关列。
2. **subworkflow 接线**：对应 `subworkflow/<Workflow>.smk` 的
   `qc_summary_config` 加三处——`qc_dirs` 加 `"<stage>": <producer>_config["outdir"]`
   （引用上游 config dict，禁止重拼路径字符串）；产物被 node.py 注册或必然
   生成时同步加 `qc_anchors` 的 `{sample}` 模板锚点（保证排序与重跑联动），
   产物可能缺席（如按 assay 门控）时锚点用条件语句注入、不得无条件声明；
   锚点只对部分样本存在时加 `qc_anchor_samples` 覆盖。
3. **指标登记**：在 `modules/qc_summary/metrics.md` 对应 stage 章节登记每个
   新列：指标名（=tsv 列名）、含义、计算原理（原始文件哪行哪字段、如何换算）、
   反映什么/怎么看（含典型范围与阈值）；stage 概览表同步加一行。三要素不齐
   全不算收录完成（metrics.md 顶部既有约定，此处升格为模块级强制规范）。
4. **列序与出图**：核心指标加入 `build_table()` 的 `core_order`（数值核心项
   排前），需要出图的加 `render_plots()` 分支并扩 `--plot` choices。

node.py 侧按 assay/配置门控该 QC 规则的 outfiles 注册（如 exome 才有
collectHsMetrics），WGS 等不适用场景的列为 NA 属预期，不构成接入缺失。
消费方（报告着色阈值等）同步要求见 metrics.md 的「指标登记规范」。
