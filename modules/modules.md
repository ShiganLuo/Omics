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

每个模块封装**一个工具**或**一个原子分析步骤**，不允许将多个独立分析合并到一个模块中。

判断标准：
- 该模块是否可以被独立复用（不依赖其他模块的输出）？如果不能，说明它太大了。
- 该模块是否包含两个不同的可执行工具？如果是，拆成两个模块。
- 该模块的 conda 环境是否混入了不相关的依赖？如果是，说明职责不单一。

示例：
- `telomere` 模块只调用 `telogator2`，conda 环境只有 `telogator2`
- `centromere` 模块包含 `hifiasm` + `repeatmasker` 两步（组装→注释），但它们是同一个分析的上下游步骤，不可独立复用，所以合为一个模块是合理的
- 将 telomere 和 centromere 合为一个模块是**错误的**：两者工具链完全不同，可独立使用

# 基础结构

## 三种目录布局

### A. 简单 3 文件模块（最常见）

```
modules/<tool>/
  <tool>.smk
  <tool>.json
  <tool>.yaml
```

代表：samtools、fastqc、deepvariant、hiphase、pbsv、tabix、trgt

### B. 3 文件 + bin 目录（自定义脚本）

```
modules/<tool>/
  <tool>.smk
  <tool>.json
  <tool>.yaml
  bin/
    summarize.py
    plot.py
```

代表：StringTie、arriba、spectrum、RmrRNA、track、bowtie2

`bin/` 存放该模块需要的 Python/R 辅助脚本，通过 `ROOT_DIR + "/modules/<tool>/bin/<script>.py"` 引用。

### C. 父目录 + 子模块目录（共享 conda 环境）

```
modules/<tool>/
  <tool>.yaml               # 共享 conda 环境
  <tool>.smk                # 主规则或 prepare 规则
  <tool>.json
  <subtool>/
    <subtool>.smk           # 子工具规则，conda: "../<tool>.yaml"
    <subtool>.json
```

代表：gatk（gatk_prepare.smk + gatk_bqsr/ + gatk_germline/ + gatk_somatic/ + gatk_RNAseq/）、samtools（samtools.smk + sort/）、msisensor-pro（tumor-normal/ + tumor-only/）

子目录规则的 `conda:` 必须引用父目录的 yaml：`conda: "../<tool>.yaml"`

### D. 公共模块（共享工具函数）

```
modules/common/
  common.smk
```

提供所有模块共享的工具函数（如日志、路径处理等）。详见 [common 模块](#common-模块)。

# 统一配置入口

模块通过 `config.get()` 读取参数，常见字段：

| 字段 | 用途 | 使用场景 |
|---|---|---|
| `indir` | 输入目录 | 所有模块 |
| `outdir` | 输出目录 | 所有模块 |
| `logdir` | 日志目录 | 所有模块 |
| `samples` | 样本 ID 列表 | 大多数模块 |
| `paired_samples` | 双端样本列表 | 比对模块 |
| `single_samples` | 单端样本列表 | 比对模块 |
| `sample_groups` | 样本分组字典 | StringTie |
| `sample_somatic_vcf_dict` | 样本→体细胞 VCF 映射 | spectrum |
| `sample_group_dict` | 样本→组映射 | spectrum |
| `bam_dir` / `vcf_dir` | 独立 BAM/VCF 目录 | hiphase |
| `Procedure` | 工具可执行路径字典 | 所有模块 |
| `Params` | 工具参数字典 | 需要额外参数的模块 |
| `genome` | 参考文件路径字典 | 需要参考基因组的模块 |
| `ROOT_DIR` | 项目根目录，用于引用 bin/ 脚本 | 有自定义脚本的模块 |

# 命名约定

- 规则名描述动作：`*_index` / `*_align` / `*_run` / `*_phase`
- 输出目录按模块名组织：`{outdir}/<module>/{sample_id}/...`
- 收尾规则命名：`<tool>_result`，仅声明 `input:` 作为依赖聚合点（无 `shell`/`run`）
- 报告规则命名：`<tool>_report`，聚合所有样本输出生成跨样本报告

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
    params:
        <tool> = config.get("Procedure", {}).get("<tool>") or "<tool>"
    run:
        try:
            open(log[0], "w").close()
            logger = setup_logger(logger_name="<tool>_<action>", log_file=log[0])
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start <tool> <action> for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/<tool>_<action>_{current_time}.sh")
            cmd = [
                params.<tool>, "<action>",
                ...
            ]
            with open(script, "w") as f:
                f.write(" ".join(cmd) + "\n")
            shell(f"bash {script} >> {log[0]} 2>&1")
        except Exception as e:
            with open(log[0], "a") as f:
                f.write(f"<tool> <action> failed for sample {wildcards.sample_id} with error: {e}\n")
            raise f"<tool> <action> failed for sample {wildcards.sample_id} with error: {e}"
        finally:
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"successfully activated  <rule> for sample {wildcards.sample_id} at {current_time}")
```

## 要点

1. **`open(log[0], "w").close()`** — 清空旧日志，避免追加混淆
2. **`setup_logger`** — 从 common.smk 导入，统一日志格式
3. **`current_time` 时间戳** — 脚本名加时间戳避免并发写冲突
4. **cmd 列表构建** — 参数逐项添加，条件参数用 `if` 追加，不拼接字符串
5. **`shell(f"bash {script} >> {log[0]} 2>&1")`** — stdout 和 stderr 都追加到日志
6. **try/except** — 捕获异常，在日志末尾写入错误信息后 re-raise

## 带条件参数

```python
cmd = [params.tool, "--required-flag", input.file]
if params.optional_param:
    cmd += ["--optional", params.optional_param]
```

## 多命令序列

```python
with open(script, "w") as f:
    f.write(" ".join(cmd1) + "\n")
    f.write(" ".join(cmd2) + "\n")
```

## Python 数据处理 + shell 命令混合

当需要先用 Python 准备文件再调用外部工具时，Python 逻辑直接写在 `run:` 块中，shell 命令仍通过脚本执行：

```python
run:
    # Python 准备
    with open(map_file, "w") as f:
        for k, v in mapping.items():
            f.write(f"{k}\t{v}\n")
    # shell 执行
    with open(script, "w") as f:
        f.write(" ".join(cmd) + "\n")
    shell(f"bash {script} >> {log[0]} 2>&1")
```

## 为什么不使用纯 `shell:` 块

- 纯 `shell:` 块无法记录结构化日志（时间戳、样本名、开始/结束）
- 纯 `shell:` 块无法在失败时追加错误信息到日志
- 纯 `shell:` 块无法动态构建参数列表
- 统一风格降低维护成本，所有规则可预期相同的行为模式

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

# 跨样本聚合规则

当模块产生逐样本输出需要跨样本汇总时：

```python
rule arriba_report:
    input:
        passed_fusions = expand(outdir + "/{sid}/{sid}_passed_fusions.tsv", sid=samples),
        discarded_fusions = expand(outdir + "/{sid}/{sid}_discarded_fusions.tsv", sid=samples)
    output:
        report = outdir + "/../arriba_report/arriba_fusion_report.html"
    log:
        logdir + "/all/arriba_report.log"
    conda:
        "arriba.yaml"
    params:
        summary_script = os.path.join(ROOT_DIR, "modules/arriba/bin/summarize_arriba_fusions.py")
    run:
        current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
        script = f"{outdir}/arriba_report.{current_time}.sh"
        cmd = [
            "python", params.summary_script,
            "-p", ",".join(input.passed_fusions),
            "-d", ",".join(input.discarded_fusions),
            "-o", outdir + "/../arriba_report"
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")
```

关键：
- 聚合脚本放在模块 `bin/` 目录下
- 使用 `-p` / `-d` 传递逗号分隔的文件列表（不用 `--indir` 扫描）
- 输出到子目录（如 `<tool>_report/`），不输出单个扁平文件
- 规则的 `output` 指向一个代表性文件用于依赖追踪

# 收尾规则（_result）

每个模块建议提供 `<tool>_result` 规则，仅声明 `input:` 作为依赖聚合点，供 subworkflow 的 `use rule ... as ...` 引用：

```python
rule deepvariant_result:
    input:
        vcf = outdir + "/{sample_id}/{sample_id}.vcf.gz",
        tbi = outdir + "/{sample_id}/{sample_id}.vcf.gz.tbi"
```

```python
rule gatk_bqsr_result:
    input:
        bam = outdir + "/{genome}/gatk/bqsr/{sample_id}.sorted.markdup.BQSR.bam",
        bai = outdir + "/{genome}/gatk/bqsr/{sample_id}.sorted.markdup.BQSR.bam.bai"
```

# 输出与中间文件

- 结果输出固定落在 `outdir` 目录树
- 中间产物可用 `temp(...)` 标记
- `directory(...)` 用于索引目录的规则输出
- 调试标记文件：`touch {output.flag}` 作为规则完成标记（如 fastqc）

# 配置字典约定（.json）

## 标准模板（单样本列表）

```json
{
    "indir": "input",
    "outdir": "output",
    "logdir": "logs",
    "samples": [],
    "Procedure": {
        "<tool>": null
    },
    "Params": {
        "<tool>": {}
    },
    "genome": {
        "fasta": null
    }
}
```

## 比对模块模板（双/单端样本）

```json
{
    "indir": "input",
    "outdir": "output",
    "logdir": "logs",
    "paired_samples": [],
    "single_samples": [],
    "Procedure": {
        "<tool>": null
    },
    "Params": {
        "<tool>": {}
    },
    "genome": {
        "fasta": null,
        "gtf": null
    }
}
```

## 非标准输入模块模板

```json
{
    "indir": "input",
    "outdir": "output",
    "logdir": "log",
    "sample_somatic_vcf_dict": {},
    "sample_group_dict": {},
    "genome": {
        "fasta": ""
    }
}
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

> **注意**：channel 顺序必须是 `conda-forge` → `bioconda` → `defaults`。bioconda 依赖 conda-forge 的包，顺序颠倒会导致依赖解析失败。

# Common 模块

## 用途

`modules/common/common.smk` 提供所有模块共享的工具函数，避免代码重复。

## 包含内容

```python
# 标准库
import sys
import os
import time
import shutil

# 从 config 获取
ROOT_DIR = config.get("ROOT_DIR", ".")

# 从 src/common 导入
from common.LogUtil import setup_logger
```

## 使用方式

### 在模块中 include

```python
# modules/<tool>/<tool>.smk
include: "../common/common.smk"

# 现在可以使用 setup_logger, time, shutil 等
```

### 在子目录模块中 include

```python
# modules/<tool>/<subtool>/<subtool>.smk
include: "../../common/common.smk"
```

### 在子工作流中 include

```python
# subworkflow/<Workflow>.smk
include: "../modules/common/common.smk"
```

### 在规则中使用

```python
rule your_rule:
    input: "input.txt"
    output: "output.txt"
    log: "logs/your_rule.log"
    run:
        # setup_logger 已从 common.smk 导入
        open(log, "w").close()
        logger = setup_logger(logger_name="your_rule", log_file=log)
        
        try:
            logger.info("Processing...")
            shell("some command > {log} 2>&1")
        except Exception as e:
            logger.error(f"Failed: {e}")
            raise e
```

## 优势

1. **避免代码重复**: 定义一次，到处使用
2. **统一日志**: 所有模块使用相同的 logger 配置
3. **易于维护**: 更新日志逻辑只需改一处
4. **依赖清晰**: 每个模块明确声明需要的工具

## 故障排除

### Import Error: `No module named 'common'`

如果出现此错误：
1. 确保 `ROOT_DIR` 在配置中正确设置
2. 检查 `src/common/LogUtil.py` 是否存在
3. 验证 include 路径相对于 snakefile 是否正确

### 路径问题

common 模块会自动将 `ROOT_DIR/src` 添加到 `sys.path`。如果有问题：
1. 检查 `ROOT_DIR` 是否为绝对路径
2. 验证 `src` 目录结构：
   ```
   workflow/Omics/
   ├── src/
   │   └── common/
   │       └── LogUtil.py
   ├── modules/
   │   └── common/
   │       └── common.smk  ← 此文件
   └── subworkflow/
       └── <Workflow>.smk
   ```

# Pitfalls

## 1. 子目录 conda 路径

子目录规则必须使用 `conda: "../<parent>.yaml"`（相对路径到父目录的 yaml），不能写 `conda: "../modules/<tool>/<tool>.yaml"`。

## 2. run 块中的 shell 脚本路径

使用 `current_time` 时间戳避免并行冲突：
```python
current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
script = f"{outdir}/{wildcards.sample_id}/tool_{current_time}.sh"
```

## 3. logger 来源

两种引入方式都存在，保持一致即可：
```python
from snakemake.logging import logger          # 推荐
# 或
import logging; logger = logging.getLogger()   # star 模块使用
```

## 4. outfiles 路径一致性

run.py 中 `outfiles` 的路径必须与规则的 `output` 完全匹配，否则 Snakemake 静默跳过或全量重建。

## 5. expand() 中的通配符转义

在 `expand()` 内引用外部通配符时用双花括号 `{{}}`：
```python
expand(outdir + "/{case}/{case}.{{genome_version}}.msisensor", case=samples.index)
```

## 6. 旧式模块（annovar 风格）

存在少量旧式模块直接用 `configfile:` 加载 YAML、`config["key"]` 取值。新模块应使用 `config.get()` 方式，不写 `configfile:`。
