---
name: modules编写规范
description: 撰写Snakemake模块时需要遵循此套规范
---

# 版本记录

- 2026-05-28: 初始化规范，补充模块目录组成与子目录 conda 继承说明

# 目的

统一模块的输入输出、参数命名和日志/环境配置方式，确保模块可被 subworkflow 复用。

# 适用范围

- 本文规范适用于 `workflow/Omics/modules/` 下的所有模块。
- 模块通过 subworkflow 的 `module` + `use rule` 引用，不直接读取 `run.py`。
- 模块只消费 subworkflow 传入的 `config`。

# 基础结构

0) 目录文件组成

- 每个模块目录包含：`<module>.smk` + `<module>.json` + `<module>.yaml`
	- `.smk`：流程规则
	- `.json`：模块配置模板
	- `.yaml`：conda 环境
- 若模块存在子目录：子目录通常只有 `.smk` + `.json`
	- 子目录规则的 `conda:` 引用父目录的 `.yaml`

1) 统一配置入口

- `indir` / `outdir` / `logdir`
- `paired_samples` / `single_samples` / `samples`
- `Procedure` / `Params`
- `genome` / `ROOT_DIR`

2) 必须有可复用的 `rule`

- 规则命名清晰、功能单一（如 `star_index` / `star_align`）
- 输出路径固定在 `outdir` 下

3) 统一日志与环境

- `log:` 指向 `logdir` 下的规则日志
- `conda:` 使用模块目录内的 yaml
- 若是子目录规则，`conda:` 指向父目录的 yaml

# 命名约定

- 规则名描述动作：`*_index` / `*_align` / `*_result`
- 输出目录按模块名组织：`{outdir}/<module>/<sample_id>/...`
- 使用 `*_result` 作为模块收尾标记（可选）

# 配置字典约定

模块内只使用 `config` 读取参数，建议字段如下：

- `indir` / `outdir` / `logdir`
- `Procedure`：工具名或可执行路径
- `Params`：参数字典
- `genome`：基因组文件路径
- `paired_samples` / `single_samples` / `samples`
- `ROOT_DIR`：需要运行项目脚本时使用

# 常见模式示例

## 1) Index + Align 双规则（来自 star.smk / hisat2.smk）

```python
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")
indir = config.get("indir", "input")
paired_samples = config.get("paired_samples", [])
single_samples = config.get("single_samples", [])

rule star_index:
	input:
		fasta = config.get("genome", {}).get("fasta"),
		gtf = config.get("genome", {}).get("gtf")
	output:
		index_file = directory(outdir + "/index")
	log:
		logdir + "/index/star_index.log"
	conda:
		"star.yaml"
	params:
		STAR = config.get("Procedure", {}).get("STAR") or "STAR",
		index_dir = outdir + "/index"
	shell:
		"""
		mkdir -p {params.index_dir}
		{params.STAR} --runMode genomeGenerate \
			--genomeDir {params.index_dir} \
			--genomeFastaFiles {input.fasta} \
			--sjdbGTFfile {input.gtf} \
			> {log} 2>&1
		"""
```

## 2) 动态输入选择（来自 fastqc.smk / star.smk）

```python
def get_alignment_input(wildcards):
	paired_r1 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_1.fq.gz"
	paired_r2 = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}_2.fq.gz"
	single = f"{indir}/{wildcards.sample_id}/{wildcards.sample_id}.single.fq.gz"

	if wildcards.sample_id in paired_samples:
		return [paired_r1, paired_r2]
	if wildcards.sample_id in single_samples:
		return [single]
	raise ValueError(f"Sample {wildcards.sample_id} not defined")
```

## 3) 合并输出（来自 TEtranscripts.smk / featureCounts.smk）

```python
def get_cnt_tables(wildcards):
	return [f"{outdir}/TEcount/{sid}.TEcount.cntTable" for sid in samples]

rule combine_TEcount:
	input:
		fileList = get_cnt_tables
	output:
		outfile = outdir + "/TEcount/all_TEcount.tsv"
	conda:
		"TEtranscripts.yaml"
	params:
		combineTE = ROOT_DIR + "/modules/TEtranscripts/bin/combineTE.py",
		indir = outdir + "/TEcount"
	shell:
		"""
		python {params.combineTE} -p TEcount -i {params.indir} -o {output.outfile}
		"""
```

# 输出与中间文件

- 结果输出固定落在 `outdir` 目录树
- 中间产物可用 `temp(...)` 标记
- `directory(...)` 用于索引目录的规则输出

# 最小模板

```python
from snakemake.logging import logger

indir = config.get("indir", "input")
outdir = config.get("outdir", "output")
logdir = config.get("logdir", "log")

rule example:
	input:
		infile = f"{indir}/{{sample_id}}.txt"
	output:
		outfile = f"{outdir}/{{sample_id}}.out.txt"
	log:
		logdir + "/{sample_id}/example.log"
	conda:
		"example.yaml"
	params:
		tool = config.get("Procedure", {}).get("example") or "example"
	shell:
		"""
		{params.tool} -i {input.infile} -o {output.outfile} > {log} 2>&1
		"""
```
