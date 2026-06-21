---
name: subworkflow编写规范
description: 撰写subworkflow时需要遵循此套规范
---

# 目的

本规范用于统一 Snakemake 子工作流（subworkflow）的结构、变量命名、模块引用方式和参数传递方式，保证可复用、可维护、可扩展。

# 通用结构

1) 统一基础变量与默认值

- `indir` / `outdir` / `logdir`
- `paired_samples` / `single_samples`
- `outfiles`
- `ROOT_DIR` — 项目根路径，必须传入每个模块配置

说明：当前所有配置由项目根目录的 `run.py` 统一生成并下发给子工作流消费。

2) `rule all` 只依赖 `outfiles`

3) 以模块化方式组织流程

- `module` + `use rule` 复用模块规则
- 每个模块配置以 `*_config` 字典命名
- `Procedure` 与 `Params` 从 `config` 提取

4) 需要条件分支时，用 `aligner` / `trimmer` 等参数控制

5) 每个模块配置记录日志（`logger.info`）

注意：暂不支持在子工作流内部拼接/串联其他 subworkflow；如需组合流程，请在 `run.py` 侧编排并生成 `outfiles`。

# 命名约定

- 子工作流中导入规则统一加前缀，避免重名：
  - 例：`use rule star_align from star as RNAseq_star_align`
- 配置对象统一 `_config` 后缀：`star_config_for_TEtranscripts`
- `Procedure` / `Params` 使用 `config.get(...)` 读取，保证容错

# 配置字典规范

推荐包含以下字段（按需选用）：

- `indir` / `outdir` / `logdir`
- `paired_samples` / `single_samples` / `samples`
- `Procedure`：可执行程序名或工具路径
- `Params`：工具参数
- `genome`：基因组相关路径（`fasta` / `gtf` / `index_dir` 等）
- `ROOT_DIR`：项目根路径，**所有模块配置必须传入**（common.smk 依赖此变量定位 `src/` 和 `bin/` 脚本）

# 典型模式示例

## 1) 条件分支选择 trimmer（来自 RNAseq.smk）

```python
if trimmer == "cutadapt":
	cutadapt_config = {
		"indir": indir,
		"outdir": f"{outdir}/fastq/cutadapt",
		"logdir": logdir,
		"ROOT_DIR": ROOT_DIR,
		"Procedure": {
			"trim_galore": config.get("Procedure", {}).get("trim_galore")
		}
	}
	module cutadapt:
		snakefile: "../modules/cutadapt/cutadapt.smk"
		config: cutadapt_config
	use rule trimming_Paired from cutadapt as RNAseq_trimming_Paireds
	use rule trimming_Single from cutadapt as RNAseq_trimming_Single
elif trimmer == "trimmomatic":
	trimmomatic_config = {
		"indir": indir,
		"outdir": f"{outdir}/trimmomatic",
		"logdir": logdir,
		"ROOT_DIR": ROOT_DIR,
		"Procedure": {
			"trimmomatic": config.get("Procedure", {}).get("trimmomatic")
		},
		"Params": {
			"trimmomatic": {
				"adapter_pe": config.get("Params", {}).get("trimmomatic", {}).get("adapter_pe"),
				"adapter_se": config.get("Params", {}).get("trimmomatic", {}).get("adapter_se")
			}
		}
	}
	module trimmomatic:
		snakefile: "../modules/trimmomatic/trimmomatic.smk"
		config: trimmomatic_config
	use rule trimmomatic_Paired from trimmomatic as RNAseq_trimmomatic_Paireds
	use rule trimmomatic_Single from trimmomatic as RNAseq_trimmomatic_Singles
else:
	raise ValueError(f"Unsupported trimmer: {trimmer}")
```

## 2) 条件分支选择比对器（来自 CLIP.smk / RNAseq.smk）

```python
if aligner == "hisat2":
	hisat2_config = {
		"indir": cutadapt_config["outdir"],
		"outdir": f"{outdir}/hisat2",
		"logdir": logdir,
		"paired_samples": paired_samples,
		"single_samples": single_samples,
		"ROOT_DIR": ROOT_DIR,
		"Procedure": {
			"hisat2": config.get("Procedure", {}).get("hisat2")
		},
		"genome": {
			"fasta": config.get("genome", {}).get("fasta"),
			"index_dir": config.get("genome", {}).get("hisat2_index_dir")
		}
	}
	module hisat2:
		snakefile: "../modules/hisat2/hisat2.smk"
		config: hisat2_config
	use rule hisat2_align from hisat2 as CLIP_hisat2_align
	use rule hisat2_index from hisat2 as CLIP_hisat2_index
elif aligner == "star":
	star_config = {
		"indir": cutadapt_config["outdir"],
		"outdir": f"{outdir}/star",
		"logdir": logdir,
		"paired_samples": paired_samples,
		"single_samples": single_samples,
		"ROOT_DIR": ROOT_DIR,
		"Procedure": {
			"star": config.get("Procedure", {}).get("star")
		},
		"genome": {
			"fasta": config.get("genome", {}).get("fasta"),
			"gtf": config.get("genome", {}).get("gtf"),
			"index_dir": config.get("genome", {}).get("star_index_dir")
		}
	}
	module star:
		snakefile: "../modules/star/star.smk"
		config: star_config
	use rule star_align from star as CLIP_star_align
	use rule star_index from star as CLIP_star_index
else:
	raise ValueError(f"Unsupported aligner: {aligner}")
```

## 3) 多模块串联（来自 Mutation.smk）

```python
module bwa_mem2:
	snakefile: "../modules/bwa-mem2/bwa-mem2.smk"
	config: bwa_mem2_confg
use rule bwaMem2_index from bwa_mem2 as Mutation_bwaMem2_index
use rule bwaMem2_alignment from bwa_mem2 as Mutation_bwaMem2_alignment

module gatk_prepare:
	snakefile: "../modules/gatk/gatk_prepare.smk"
	config: gatk_prepare_config
use rule addReadsGroup from gatk_prepare as Mutation_addReadsGroup
use rule MarkDuplicates from gatk_prepare as Mutation_MarkDuplicates
```

## 4) 外部模块批量导入（来自 RNAseq.smk / CoCulture.smk）

```python
module TEtranscripts:
	snakefile: "../modules/TEtranscripts/TEtranscripts.smk"
	config: TEtranscripts_config
use rule * from TEtranscripts as RNAseq_*
```

# 输出文件规范

- `outfiles` 应在最外层由主 run.py 生成并传入 subworkflow
- 子工作流只消费 `outfiles`，不负责全局拼接
- 如需跨 subworkflow 聚合输出，由 `run.py` 统一编排
- 输出目录统一挂载在 `outdir` 下，保持层级清晰

# 最小模板

```python
shell.prefix("set -x; set -e;")
from snakemake.logging import logger

indir = config.get("indir", "data/fastq")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "logs")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])
outfiles = config.get("outfiles", [])
ROOT_DIR = config.get("ROOT_DIR", ".")

rule all:
	input:
		outfiles

module example:
	snakefile: "../modules/example/example.smk"
	config: {
		"indir": indir,
		"outdir": f"{outdir}/example",
		"logdir": logdir,
		"ROOT_DIR": ROOT_DIR,
		"Procedure": {
			"example": config.get("Procedure", {}).get("example")
		}
	}
use rule run_example from example as MyWorkflow_example
```

# Pitfalls

## 1. 忘记传 ROOT_DIR

每个模块配置必须包含 `"ROOT_DIR": ROOT_DIR`。common.smk 通过 `config.get("ROOT_DIR", ".")` 读取此值，用于定位 `src/common/LogUtil.py` 和模块的 `bin/` 脚本。如果不传，默认值 `"."` 会导致路径错误（找不到脚本或导入失败）。

## 2. 子目录模块的 conda 路径

子目录规则必须使用 `conda: "../<parent>.yaml"`（相对路径到父目录的 yaml），不能写 `conda: "../modules/<tool>/<tool>.yaml"`。