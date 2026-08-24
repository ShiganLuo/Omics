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

每个模块封装**一个工具**或**一个原子分析步骤**，不允许将多个独立分析合并到一个模块中，另外规则只负责校验和调用。

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

### C. 父目录 + 子模块目录

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

子模块默认引用父目录的 `<tool>.yaml`（`conda: "../<tool>.yaml"`）。当子模块需要不同环境时，在自己的目录下放 `<subtool>.yaml`，文件名 stem 必须唯一（见下文"环境命名规范"），不能多个子模块共用同名 yaml。

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

- 所有规则统一使用 `run:` 块，不允许使用纯 `shell:` 块。命令通过列表构建，写入 `.sh` 脚本后执行。

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
        <tool> = config.get("Procedure", {}).get("<tool>") or "<tool>"
    run:
        log_path = str(log)
        try:
            open(log_path, "w").close()
            logger = setup_logger(logger_name="<tool>_<action>", log_file=log_path)
            current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
            logger.info(f"Start <tool> <action> for sample {wildcards.sample_id} at {current_time}")
            script = os.path.join(outdir, f"{wildcards.sample_id}/<tool>_<action>_{current_time}.sh")
            cmd = [
                params.<tool>, "<action>",
                ...
            ]
            with open(script, "w") as f:
                f.write(shlex.join(cmd) + "\n")
                f.write(f'echo "<tool>_<action> called for {wildcards.sample_id} at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} was successfully completed"')
            shell(f"bash {script} >> {log_path} 2>&1")
        except Exception as e:
            with open(log_path, "a") as f:
                f.write(f"<tool> <action> failed for sample {wildcards.sample_id} with error: {e}\n")
            raise RuntimeError(f"<tool> <action> failed for sample {wildcards.sample_id} with error: {e}")
```
注意：
1. 对于有些脚本接受nargs+式参数，命令列表不能["-f", "heatmap volcano pca"]的形式；而是要"-f", "heatmap", "volcano", "pca",因为shlex会把它们认为一个字符串进行传递，导致报错

## 要点

1. **`open(log_path, "w").close()`** — 清空旧日志，避免追加混淆
2. **`setup_logger`** — 从 common.smk 导入，统一日志格式
3. **`current_time` 时间戳** — 脚本名加时间戳避免并发写冲突
4. **cmd 列表构建** — 参数逐项添加，条件参数用 `if` 追加，不拼接字符串
5. **`shell(f"bash {script} >> {log_path} 2>&1")`** — stdout 和 stderr 都追加到日志
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
    log_path = str(log)
    # Python 准备
    with open(map_file, "w") as f:
        for k, v in mapping.items():
            f.write(f"{k}\t{v}\n")
    # shell 执行
    with open(script, "w") as f:
        f.write(" ".join(cmd) + "\n")
    shell(f"bash {script} >> {log_path} 2>&1")
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

# 报告模块（_report）

报告模块是特殊的跨样本聚合模块，消费上游所有模块的输出（TSV/CSV/PNG/BAM/log），生成 PPTX + XLSX 报告。不重复计算上游分析，只读取已有结果文件。

代表：RNAseq_report、ncRNAseq_report

## 目录结构

采用 B 类布局（3 文件 + bin/）：

```
modules/<Workflow>_report/
  <Workflow>_report.smk        # generate_report + report_result 两个 rule
  <Workflow>_report.json       # 配置 schema
  <Workflow>_report.yaml       # conda 环境 (python-pptx, matplotlib, pandas, openpyxl)
  <Workflow>_report.def        # Apptainer SIF 构建文件
  <Workflow>_report.dockerfile # Docker 构建文件
  bin/
    generate_report.py         # 报告生成脚本
```

## 两个 rule

```python
REPORT_SCRIPT = os.path.join(ROOT_DIR, "modules", "<Workflow>_report", "bin", "generate_report.py")

rule generate_report:
    input:
        # 用 expand() 声明所有上游依赖文件
        per_sample_bams = expand(outdir + "/common/4_per_gene_bam/{sample}/{sample}.bam", sample=samples),
        per_sample_tails = expand(outdir + "/common/4_per_gene_bam/{sample}/{sample}_tail.csv", sample=samples),
        # ... 其他上游输出
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
        # title/subtitle/pipeline/genome/date/lang 从 config["Params"]["report"] 读取
        title = config.get("Params", {}).get("report", {}).get("title") or "Default Title",
        lang = config.get("Params", {}).get("report", {}).get("lang") or "zh",
        img_dir = outdir + "/ppt_results",
        script = REPORT_SCRIPT,
    run:
        # 标准 run 块模板：日志 -> 构建命令 -> 写 .sh -> 执行
        ...

rule report_result:
    input:
        report = outdir + "/<Workflow>_report.pptx",
        file_inventory = outdir + "/<Workflow>_report_files.xlsx",
```

关键约定：
- **两个输出**：`report`（PPTX）+ `file_inventory`（XLSX），XLSX 汇总所有上游 TSV/CSV 的数据
- **`report_result`** 是纯依赖聚合点，仅声明 `input:`，供 subworkflow `use rule` 引用
- **`img_dir`**：中间图表 PNG 的持久化目录（`outdir/ppt_results/`），与临时文件分离
- **`--lang zh`**：报告支持中英双语，通过 `config["Params"]["report"]["lang"]` 控制

## generate_report.py 结构

```python
#!/usr/bin/env python3
"""Generate PPT report for <Workflow> pipeline results."""

# 1. 中文字体自动检测（WenQuanYi / Noto Sans CJK / SimHei）
# 2. I18N 字典（zh / en 两套标题、标签翻译）
# 3. 数据收集函数：解析上游 TSV/CSV/log 文件
# 4. 图表生成函数：matplotlib 柱状图/堆叠图/分布图
# 5. Slide 构建函数：python-pptx 创建标题页、图表页、数据表页、总结页
# 6. Excel 文件清单：openpyxl 写入所有上游数据 + 交叉表
# 7. main()：组装所有 slide，保存 PPTX + XLSX
```

### 图表生成

图表通过 matplotlib 生成临时 PNG，嵌入 PPTX：

```python
class TempImageStore:
    """生成临时 PNG，可选持久化到 img_dir，结束时自动清理。"""
    def save_fig(self, fig, stem: str) -> str:
        # 保存到 tempfile，同时 copy 到 img_dir
    def cleanup(self):
        # 删除临时文件
```

### Slide 布局常量

```python
SLIDE_W = 10.0       # 16:9 宽度（英寸）
SLIDE_H = 5.625      # 16:9 高度
HEADER_H = 0.65      # 标题栏高度
MARGIN_L = 0.45      # 左边距
CONTENT_W = SLIDE_W - MARGIN_L - MARGIN_R  # 内容区宽度
```

### 助手函数

| 函数 | 用途 |
|------|------|
| `_header(slide, text)` | 深色标题栏 + 白色标题文字 |
| `_textbox(slide, ...)` | 文本框，支持字号/粗体/颜色/对齐 |
| `_bullets(slide, ...)` | 项目符号列表 |
| `_table(slide, ...)` | 数据表，表头深色背景、隔行浅色 |
| `_add_picture(slide, path, ...)` | 等比缩放图片居中放入指定区域 |

### Excel 文件清单

```python
def write_file_inventory(output_path, analysis_dir, samples, ...):
    """每个上游 TSV/CSV 成为一个 sheet，包含完整数据（上限 5000 行）。"""
    # Sheet 1: Overview（分析目录、样本数、对比数等）
    # Sheet 2: Sample_Stats（每样本一行，所有统计指标）
    # Sheet 3+: 每样本的 manifest / tail CSV / 上游结果
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
    "paired_samples": paired_samples,
    "single_samples": single_samples,
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

## config 中的 report 参数

```json
{
    "Params": {
        "report": {
            "title": "ncRNAseq 分析报告",
            "subtitle": "",
            "pipeline": "FASTQ -> ... -> Report",
            "genome": "GRCh38",
            "date": "",
            "lang": "zh"
        }
    }
}
```

## 与普通聚合规则的区别

| 特征 | 普通聚合规则（如 arriba_report） | 报告模块（_report） |
|------|------|------|
| 输出 | 单个 HTML/TSV | PPTX + XLSX |
| 依赖 | 同模块的上游输出 | **跨模块**的所有上游输出 |
| 脚本 | 模块自带 `bin/` 脚本 | 独立 `bin/generate_report.py` |
| 环境 | 共享上游模块的 conda | 独立 conda（python-pptx/matplotlib/openpyxl） |
| subworkflow | `use rule` 单条 | `use rule generate_report` + `use rule report_result` 两条 |
| 图表 | 无或简单 | matplotlib 生成多张图表嵌入 PPTX |

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

> **注意**：channel 顺序必须是 `conda-forge` -> `bioconda` -> `defaults`。bioconda 依赖 conda-forge 的包，顺序颠倒会导致依赖解析失败。

## 环境命名规范（强制）

YAML 文件名、`name:` 字段、SIF 文件名必须三者一致，等于模块目录名：

| 文件 | 命名 |
|------|------|
| YAML 文件名 | `<module_dir>.yaml` |
| YAML `name:` 字段 | `<module_dir>` |
| SIF 文件名 | `<module_dir>.sif` |

示例：
```
modules/star/star.yaml       name: star       -> star.sif
modules/bedtools/bedtools.yaml  name: bedtools  -> bedtools.sif
modules/openms/psmfdr/psmfdr.yaml  name: psmfdr  -> psmfdr.sif
```

**禁止**：
- 多个子模块共享同一 YAML 文件名（如全部叫 `openms.yaml`）
- YAML `name:` 字段与文件名 stem 不一致（如 `gatk.yaml` 写 `name: gatk4`）

`common.smk` 中的 `sif()` 函数用 YAML 文件名 stem 作为 key 查找 SIF 路径，
不读取 YAML 内容。命名不一致会导致容器路径解析失败。

> **例外**：如确实需要名称与目录不同（如上游版本锁定），在 config 的 `env`
字典中显式映射：`"yaml_stem": "/path/to/actual.sif"`。这是最后手段，不应作为常规做法。

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

# 容器路径解析
_ENV_MAP = config.get("env", {})

def sif(yaml_filename: str) -> str:
    """用 YAML 文件名 stem 查找 SIF 路径。
    1. config["env"][<stem>] -- 精确映射
    2. config["env"]["env_dir"] + <module_dir>/<stem>.sif -- 约定回退
    """

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

子目录规则引用 yaml 时使用相对路径。引用父目录的共享 yaml 时用 `conda: "../<parent>.yaml"`（如 gatk_RNAseq 引用 `../gatk.yaml`）；引用子目录自己的 yaml 时用 `conda: "<subtool>.yaml"`。不能写 `conda: "../modules/<tool>/<tool>.yaml"`。

当子模块需要独立环境时，必须有自己的 `<subtool>.yaml`，且文件名 stem 必须唯一（见上文"环境命名规范"），不能多个子模块共用同名 yaml。

## 2. run 块中的 shell 脚本路径

使用 `current_time` 时间戳避免并行冲突：
```python
current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
script = f"{outdir}/{wildcards.sample_id}/tool_{current_time}.sh"
```

## 3. outfiles 路径一致性

run.py 中 `outfiles` 的路径必须与规则的 `output` 完全匹配，否则 Snakemake 静默跳过或全量重建。

## 4. expand() 中的通配符转义

在 `expand()` 内引用外部通配符时用双花括号 `{{}}`：
```python
expand(outdir + "/{case}/{case}.{{genome_version}}.msisensor", case=samples.index)
```

