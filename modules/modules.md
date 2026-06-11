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

# 三种 run/shell 风格

默认采用风格2，除非特殊指定

## 风格 1：纯 `shell:` 块（简单命令）

适用于单条命令、参数固定的工具。

```python
rule bam_flagstat:
    input:
        bam = indir + "/{sample_id}/{sample_id}.bam"
    output:
        flagstat = outdir + "/{sample_id}/{sample_id}_flagstat.txt"
    log:
        logdir + "/{sample_id}/flagstat.log"
    conda:
        "samtools.yaml"
    params:
        samtools = config.get("Procedure", {}).get("samtools") or "samtools"
    shell:
        """
        {params.samtools} flagstat {input.bam} > {output.flagstat}
        """
```

代表：samtools、fastqc、gatk_prepare（addReadsGroup、MarkDuplicates）

## 风格 2：`run:` 块 + shell 脚本生成（复杂命令）

适用于需要动态构建参数列表、条件分支、管道连接的工具。生成的 `.sh` 脚本同时方便调试和复现。

```python
rule pbsv_call:
    input:
        svsig = outdir + "/discover/{sample_id}/{sample_id}.svsig.gz",
        fasta = fasta
    output:
        vcf = outdir + "/call/{sample_id}/{sample_id}.vcf"
    log:
        logdir + "/{sample_id}/pbsv_call.log"
    threads: 8
    conda:
        "pbsv.yaml"
    params:
        pbsv = config.get("Procedure", {}).get("pbsv") or "pbsv"
    run:
        current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        logger.info(f"Start pbsv call for sample {wildcards.sample_id} at {current_time}")
        script = os.path.join(outdir, f"pbsv_call_{current_time}.sh")
        cmd = [
            params.pbsv, "call",
            "--num-threads", str(threads),
            input.fasta,
            input.svsig,
            output.vcf
        ]
        with open(script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(" ".join(cmd) + "\n")
        shell("bash {script} > {log} 2>&1")
```

代表：pbsv、deepvariant、hiphase、star、bowtie2、hisat2

### 带条件参数的变体

```python
run:
    cmd = [params.tool, "--required-flag", input.file]
    if params.optional_param:
        cmd += ["--optional", params.optional_param]
    # ...生成脚本
```

### 多命令序列的变体

```python
run:
    script = f"{outdir}/{wildcards.sample_id}/align.{current_time}.sh"
    cmd1 = [params.aligner, "-x", index, "-1", r1, "-2", r2, "|", "samtools", "view", "-bS", "-", ">", output.bam]
    cmd2 = ["samtools", "index", output.bam]
    with open(script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(" ".join(cmd1) + "\n")
        f.write(" ".join(cmd2) + "\n")
    shell("bash {script} > {log} 2>&1")
```

代表：star（align + mv + index）、bowtie2（align + mv unmapped）、hisat2（align pipe samtools）

## 风格 3：`run:` 块 + 直接 Python 逻辑

适用于需要 Python 数据处理（pandas、json、文件聚合等）的规则。

```python
rule merge_msi_result:
    input:
        expand(outdir + "tumor_normal_output/{case}/{case}.{genome_version}.msisensor",
               case=samples.index)
    output:
        outdir + "tumor_normal_output.{genome_version}.merge.tsv"
    run:
        output_info = pd.DataFrame(columns=[...])
        for case in samples.index:
            value_info = [i.split() for i in open(...)]
            output_info.loc[case] = value_info[1]
        output_info.to_csv(f"{output}", sep="\t")
```

代表：msisensor-pro（merge）、spectrum（生成映射文件 + shell）、track（写 JSON config + shell）

### 混合 Python + shell 的变体

```python
run:
    # 先用 Python 准备输入文件
    with open(map_file, "w") as f:
        for k, v in mapping.items():
            f.write(f"{k}\t{v}\n")
    # 再调用 shell 命令
    shell("python {params.script} -i {map_file} -o {output} > {log} 2>&1")
```

代表：spectrum、arriba_report、track (igv)

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
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - <tool>>=<version>
```

# Pitfalls

## 1. 子目录 conda 路径

子目录规则必须使用 `conda: "../<parent>.yaml"`（相对路径到父目录的 yaml），不能写 `conda: "../modules/<tool>/<tool>.yaml"`。

## 2. shell vs run 选择

- 简单单命令 → 用 `shell:` 块
- 需要条件参数、管道、多步命令 → 用 `run:` 块生成脚本
- 需要 Python 数据处理 → 用 `run:` 块
- 不要在 `run:` 块中直接写多行 shell 字符串，生成脚本文件更易调试

## 3. run 块中的 shell 脚本路径

使用 `current_time` 时间戳避免并行冲突：
```python
current_time = time.strftime("%Y%m%d.%H:%M:%S", time.localtime())
script = f"{outdir}/{wildcards.sample_id}/tool_{current_time}.sh"
```

## 4. logger 来源

两种引入方式都存在，保持一致即可：
```python
from snakemake.logging import logger          # 推荐
# 或
import logging; logger = logging.getLogger()   # star 模块使用
```

## 5. outfiles 路径一致性

run.py 中 `outfiles` 的路径必须与规则的 `output` 完全匹配，否则 Snakemake 静默跳过或全量重建。

## 6. expand() 中的通配符转义

在 `expand()` 内引用外部通配符时用双花括号 `{{}}`：
```python
expand(outdir + "/{case}/{case}.{{genome_version}}.msisensor", case=samples.index)
```

## 7. 旧式模块（annovar 风格）

存在少量旧式模块直接用 `configfile:` 加载 YAML、`config["key"]` 取值。新模块应使用 `config.get()` 方式，不写 `configfile:`。
