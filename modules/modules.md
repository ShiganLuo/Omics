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

# 命名约定

- 规则名描述动作：`*_index` / `*_align` / `*_run` / `*_phase`
- 输出目录按模块名组织：`{outdir}/<module>/{sample_id}/...`
- 收尾规则：`<tool>_result`，仅声明 `input:` 作为依赖聚合点（无 `shell`/`run`）
- 报告规则：`<tool>_report`，聚合所有样本输出生成跨样本报告

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
        try:
            open(log_path, "w").close()
            rule_logger = setup_logger("<tool>_<action>", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            rule_logger.info(f"Start <tool> <action> for sample {wildcards.sample_id} at {current_time}")
            sample_outdir = os.path.dirname(str(output.xxx))
            os.makedirs(sample_outdir, exist_ok=True)
            script = os.path.join(sample_outdir, f"<tool>_<action>_{wildcards.sample_id}_{current_time}.sh")
            cmd = [
                params.<tool>, "<action>",
                ...
            ]
            with open(script, "w") as f:
                f.write("#!/bin/bash\nset -euo pipefail\n")
                f.write(" ".join(cmd) + "\n")
                f.write(f'echo "<tool> <action> for {wildcards.sample_id} at {current_time} completed successfully"\n')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"Error occurred during <tool> <action> for sample {wildcards.sample_id}: {e}\n")
            logger.error(f"Error occurred during <tool> <action> for sample {wildcards.sample_id}: {e}")
            raise e
```

## 要点

1. **`open(log_path, "w").close()`** — 清空旧日志，避免追加混淆
2. **`setup_logger`** — 从 common.smk 导入，统一日志格式
3. **`rule_logger.info(...)`** — 开始/结束标记，便于定位问题
4. **`current_time` 时间戳** — 脚本名加时间戳避免并发写冲突
5. **`cmd` 列表构建** — 参数逐项添加，条件参数用 `if` 追加
6. **`" ".join(cmd)`** — 不用 `shlex.join()`
7. **`#!/bin/bash` + `set -euo pipefail`** — 脚本内任何命令失败立即退出，不会执行后续的 echo
8. **`echo` 完成标记** — 脚本末尾写 `echo "...completed successfully"`，执行成功时日志中有明确结尾
9. **`shell(f"bash {script} >> {log_path} 2>&1")`** — 直接路径，不加 `shlex.quote` 包装
10. **`logger.error(...)` + `raise e`** — snakemake logger 记录错误后 re-raise

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
        h5ad = expand(de_dir + "/{t}_de.h5ad", t=tissues),
        table = expand(de_dir + "/{t}_markers.tsv", t=tissues)
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
