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

6) 每个 subworkflow 收尾必须包含一个 **QC 汇总模块**和一个 **report 模块**（见「QC 汇总与 report 模块规范」）

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
- `env`：SIF 容器路径映射字典，**所有模块配置必须传入**（common.smk 的 `sif()` 函数依赖此字段解析 `.sif` 路径）

# 日志规范

subworkflow 顶层从 config 读取 `logdir`，下发给模块配置时**按规则粒度分流**到两个固定分段目录：

- 单样本规则 → `logdir/sample`
- 多样本（组级/跨样本）规则 → `logdir/group`

## 传法

| 模块内规则构成 | 模块配置传入 | 模块内日志落点 |
|---|---|---|
| 仅单样本规则 | `"logdir": f"{logdir}/sample"` | `{logdir}/sample/{sample_id}/{规则名}.log` |
| 仅多样本规则 | `"logdir": f"{logdir}/group"` | `{logdir}/group/{module名}/{规则名}.log` |
| 单样本 + 多样本混合 | `"logdir": f"{logdir}/sample"`<br>加 `"logdir_combine": f"{logdir}/group"` | 单样本规则用 `logdir`，多样本规则用 `logdir_combine` |

混合模块必须**同时传两个变量**（`logdir` 与 `logdir_combine`），模块侧据此把两类规则的路由分开（模块内部怎么消费见 `modules/modules.md` 的日志规范）。

## 示例（Mutation.smk，混合模块 fragment_size）

```python
fragment_size_config = {
    ...
    "logdir": f"{logdir}/sample",
    "logdir_combine": f"{logdir}/group",
    ...
}
```

## 约束

- 分段目录名固定为 `sample` / `group`，不要新造名字（不要 `logs`、`combine`、`multi` 等）。
- 不要把顶层 `logdir` 原样传给模块——那会让样本日志与组级日志混在同一级目录。
- 本规范只约束规则日志（`log:` 指令）；命令回执池 `OUTDIR_ROOT/log/cmd/` 由 `OUTDIR_ROOT` 派生，与 logdir 分段无关。

# 模块衔接规范

subworkflow 内模块按执行顺序自上而下编排。**下游模块消费上游输出时，输入目录一律引用上游 config 的 `outdir`**（`xxx_config["outdir"]`），禁止手拼同一路径的 f-string 副本——路径知识只写一次，上游改布局下游自动跟随：

```python
# ✅ 正确：indir 引用上游 config["outdir"]
cutadapt_config = {
    "outdir": f"{outdir}/fastq/cutadapt",
    ...
}
star_config = {
    "indir": cutadapt_config["outdir"],   # 上游 outdir，不重拼
    "outdir": f"{outdir}/star",
    ...
}

# ❌ 错误：同一路径出现两个手拼副本，改一处漏一处
star_config = {
    "indir": f"{outdir}/fastq/cutadapt",  # 与 cutadapt_config["outdir"] 重复定义
    ...
}
```

## 模块级常量

一个输出目录被多个下游消费时，在产出方 config **之后**定义模块级常量，下游统一引用常量：

```python
extract_config = {
    "outdir": f"{analysis_outdir}/extract",
    ...
}
elements_dir = extract_config["outdir"]      # 常量引用产出方 config

footprint_config = {
    "indir": elements_dir,
    ...
}
```

## 条件分支的衔接

上游是二选一（如 aligner）时，用条件表达式取实际选中的那个 config：

```python
# 来自 CLIP.smk
"indir": star_config["outdir"] if aligner == 'star' else hisat2_config["outdir"],
```

## 约束

- Snakefile 顺序执行：config 字典与常量**先定义后使用**，禁止前向引用后面才定义的 config（NameError）。
- 判断标准：一个物理路径在整个 subworkflow 里只出现**一次字面定义**，其余位置全部通过 `xxx_config["outdir"]`（或引用它的常量）获取。
- 衔接只打通路径；是否真的运行由 node.py 注册的 outfiles 驱动，DAG 自动解析依赖（下游 enabled 自动拉起上游）。
- 模块 .smk 内部如何消费这些输入键（`motif_dir` / `vcf_dir` 等独立 config 键）见 `modules/modules.md` 的「Subworkflow 跨模块衔接」。

# QC 汇总与 report 模块规范（强制）

背景：Mutation 流程曾只有零散的 QC 产物（fastqc raw/trimmed、flagstat 等），各自落在独立目录、无人汇总，导致交付时需要人工翻目录拼 QC 结论。为避免每个流程重复出现这个问题，规定**每个 subworkflow 收尾必须编排两个模块**：一个 QC 汇总模块 + 一个 report 模块。

## 1) QC 汇总模块（`qc_summary`）

**职责**：跨样本聚合全流程各环节的 QC 指标，产出机器可读的汇总表 + 人可读的可视化，是流程质量结论的单一入口。QC 汇总是**数据层**：收指标、对齐样本、给阈值判定，供 report 模块和人工排查消费。

**模块落地**：`modules/qc_summary/qc_summary.smk`（一个模块被所有 subworkflow 复用；各流程的指标差异通过 config 的指标目录清单表达，不复制多份模块）。subworkflow 侧命名 `<WF>_qc_summary_config`，规则别名 `<WF>_qc_summary`。

**config 字段**：

- `qc_dirs`：指标目录字典，key 为环节名（`fastqc_raw` / `fastqc_trimmed` / `flagstat` / `frip` / `insert_size` / `msi` 等），value 为该环节产物目录。**一律引用上游 config 的 `outdir`**（模块衔接规范同样适用于 QC 目录），禁止手拼重复路径
- `qc_anchors`：**DAG 锚点**模板字典（`"{sample}/{sample}_flagstat.txt"` 这类含 `{sample}` 占位符的相对路径），只声明**保证会产出**的文件（node.py 注册的输出、或规则必然写的 log/metrics）。锚点会展开成 qc_summary 规则的 `input:`，作用有二：① 排序——qc_summary 一定在上游之后跑；② 失效——上游任一锚点被刷新，qc_summary 自动重跑。**输出未注册的环节（如 fastqc）不要锚定**，否则 DAG 直接 MissingInputException；这类环节留在 `qc_dirs` 里做 best-effort 扫描（缺了降级为 NA 列，不报错）
- `qc_anchor_samples`：可选的按环节覆盖样本列表（如 MSI 锚点只存在于配对设计中的 tumor 样本），缺省用全样本：

```python
qc_summary_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": f"{outdir}/QC/qc_summary",
    "logdir": f"{logdir}/group",          # 跨样本聚合 → group 分段
    "samples": paired_samples + single_samples,
    "qc_dirs": {
        "fastqc_raw": fastqc_raw_config["outdir"],
        "fastqc_trimmed": fastqc_trimmed_config["outdir"],
        "flagstat": samtools_config["outdir"],
    },
    # 只锚定保证产出的环节；fastqc 未注册输出，不锚定
    "qc_anchors": {
        "flagstat": "{sample}/{sample}_flagstat.txt",
    },
}
module qc_summary:
    snakefile: "../modules/qc_summary/qc_summary.smk"
    config: qc_summary_config
use rule qc_summary from qc_summary as Mutation_qc_summary
```

**产出约定**：

- `{outdir}/QC/qc_summary/qc_summary.tsv` — 样本 × 指标汇总表（每行一个样本）
- `{outdir}/QC/qc_summary/qc_summary.html`（或 `.png` 集合）— 可视化
- 允许按流程追加环节专属指标列；缺该环节产物时对应列留空（`NA`），不报错——模块对缺失目录必须容错，因为可选环节（msi/cnvkit/burden 等）是否运行由 node.py 决定。

**日志**：QC 汇总是组级规则（一次消费全部样本），`logdir` 传 `f"{logdir}/group"`。

## 2) report 模块（`<WF>_report`）

**职责**：面向交付的最终汇总报告（PPT/HTML），聚合 QC 结论 + 各分析环节的关键结果图表。report 是**呈现层**：只消费已注册的产出文件（QC 汇总表、spectrum/burden/fusion 等结果图），不重算任何指标。

**模块落地**：每个流程一个专属模块 `modules/<WF>_report/`（内容因流程而异，无法共用），规则名 `generate_report`，收尾提供 `<WF>_report_result` 面板规则（仅 `input:`，见 modules.md 收尾规则约定）。subworkflow 侧命名 `<WF>_report_config`，规则别名 `<WF>_generate_report`。

**config 字段**：

```python
Mutation_report_config = {
    "ROOT_DIR": ROOT_DIR,
    "OUTDIR_ROOT": outdir,
    "env": config.get("env", {}),
    "outdir": outdir,
    "logdir": f"{logdir}/group",
    "samples": paired_samples + single_samples,
    "qc_summary": qc_summary_config["outdir"],   # 引用 QC 汇总模块 outdir
    "spectrum_dir": somatic_spectrum_config["outdir"],
    "burden_dir": burden_config["outdir"],
    "Params": {
        "report": config.get("Params", {}).get("report", {})
    }
}
module Mutation_report:
    snakefile: "../modules/Mutation_report/Mutation_report.smk"
    config: Mutation_report_config
use rule generate_report from Mutation_report as Mutation_generate_report
```

- 输入目录一律引用上游 config 的 `outdir` 常量；report 的输入项允许超集（模块内对缺失目录容错），可选环节未注册时报告对应章节自动缺省。
- **产出必须是 PPT + Excel 双件**（对齐 PeakCalling_report / ncRNAseq_report 现状）：
  - `{outdir}/<WF>_report.pptx` — 汇报用幻灯片（QC 结论 + 各环节关键结果图表）
  - `{outdir}/<WF>_report.xlsx` — 结果收录工作簿，交付后可离线检索全部数据。固定 sheet 结构：
    1. `Overview` — 运行元信息（pipeline/genome/日期/样本/配对）
    2. `QC Summary` — qc_summary.tsv 整表收录
    3. 每个分析模块一个 sheet（变异/峰/表达等，逐行收录该模块的核心表格数据；表大的收录汇总层，行级明细保留在原文件）
    4. `File Inventory`（收尾） — 关键交付文件清单：module / path / exists / bytes，供核对产物完整性
  - 模块未运行时对应 sheet 写一行 note 占位，**不报错、不省略 sheet**（交付物结构稳定）
- Excel 由 `generate_report` 同一条规则产出（`--excel-output` 参数），不另设规则；node.py 把 pptx 和 xlsx 一起注册进 outfiles。

## 3) node.py 侧注册（两模块缺一不可）

- QC 汇总：`outfiles.append(f"{outdir}/QC/qc_summary/qc_summary.tsv")`
- report：受 `Params.report.enabled` 门控（`_enabled(datajson, "report")`），**pptx 和 xlsx 成对注册**（`f"{outdir}/<WF>_report.pptx"` + `.xlsx`）。report enabled 时自动拉起上游（DAG 反推），QC 汇总随之运行。
- 两个模块**无条件接线**（subworkflow 里始终 `module` + `use rule`），是否真的运行由 outfiles 注册驱动——与「输出决定流程走向」铁律一致，不要用 if 包裹模块定义。

## 4) sexcheck 模块（性别核查，产 BAM 流程必备）

`modules/sexcheck/`：通过比对 BAM 在性染色体上的覆盖度推断样本性别（XY/XX/unknown），防样本混淆。所有产全基因组 BAM 的 subworkflow 必须接入：

- **接线**：subworkflow 里 `sexcheck_config` 的 `sample_bams` 逐样本指向**比对/markdup 产出的 BAM**（引用上游 config 的 `outdir`，per-genome 流程用 `genome/sample` 元组展开）；`logdir` 传 sample 分段、`logdir_combine` 传 group 分段（模块含单样本 idxstats 规则 + 组级判型规则）
- **判型**：X/Y 覆盖度按染色体长度归一化（reads/Mb）后取比值；XY= Y 有覆盖且比值低，XX= Y 近零，unknown= 参考无 chrY 或落在模糊带。阈值可经 `Params.sexcheck.{ratio_xy_threshold, y_min_frac}` 覆盖（仅接受数值，null 走脚本默认）
- **输出**：`{outdir}/QC/sexcheck/sexcheck.tsv`（样本 × X/Y 覆盖度指标 + sex 列）；node.py 无条件注册该输出
- **容器**：env 映射复用同时含 samtools + pandas 的 SIF（当前 fragment_size.sif）
- 无全基因组 BAM 的流程（scRNAseq/spatial/QuantMS/tRNAseq）不适用，不接

## 5) 检查清单（新增流程 / 审查存量流程时逐条过）

- [ ] subworkflow 末尾有 `qc_summary_config` + `module qc_summary` + `use rule ... as <WF>_qc_summary`
- [ ] subworkflow 末尾有 `<WF>_report_config` + `module <WF>_report` + `use rule generate_report ... as <WF>_generate_report`
- [ ] `qc_dirs` / report 输入目录全部引用上游 config `["outdir"]`，无手拼重复路径
- [ ] 保证产出的 QC 环节都声明了 `qc_anchors`（qc_summary 的 input 依赖上游，上游刷新自动重跑汇总）；未注册输出的环节只进 `qc_dirs` 不锚定
- [ ] QC 汇总输出注册进 outfiles；report 输出（pptx + xlsx 成对）经 `Params.report.enabled` 门控注册
- [ ] 两模块日志走 `logdir/group`
- [ ] 产全基因组 BAM 的流程接入了 sexcheck（`sample_bams` 引用比对产出，node.py 注册 sexcheck.tsv）
- [ ] 存量流程缺口记录在案：Mutation 已补齐（qc_summary + Mutation_report，2026-11）；CoCulture / KARRseq / QuantMS / tRNAseq / PacVar / Population_genomics / Fiberseq / spatial 等（无 QC 汇总或无 report，逐个补齐）
- [ ] subworkflow/README.md 的流程简介同步列出这两个收尾模块

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
env = config.get("env", {})

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
		"env": env,
		"Procedure": {
			"example": config.get("Procedure", {}).get("example")
		}
	}
use rule run_example from example as MyWorkflow_example

# 收尾：QC 汇总模块（modules/qc_summary，全流程共用）+ report 模块（modules/<WF>_report）
example_config = {}  # 占位：上文各模块 config，供 qc_dirs / report 引用 outdir

qc_summary_config = {
	"ROOT_DIR": ROOT_DIR,
	"OUTDIR_ROOT": outdir,
	"env": env,
	"outdir": f"{outdir}/QC/qc_summary",
	"logdir": f"{logdir}/group",
	"samples": paired_samples + single_samples,
	"qc_dirs": {
		"example": f"{outdir}/example",   # 实际写法：example_config["outdir"]
	},
}
module qc_summary:
	snakefile: "../modules/qc_summary/qc_summary.smk"
	config: qc_summary_config
use rule qc_summary from qc_summary as MyWorkflow_qc_summary

MyWorkflow_report_config = {
	"ROOT_DIR": ROOT_DIR,
	"OUTDIR_ROOT": outdir,
	"env": env,
	"outdir": outdir,
	"logdir": f"{logdir}/group",
	"samples": paired_samples + single_samples,
	"qc_summary": qc_summary_config["outdir"],
	"Params": {
		"report": config.get("Params", {}).get("report", {})
	}
}
module MyWorkflow_report:
	snakefile: "../modules/MyWorkflow_report/MyWorkflow_report.smk"
	config: MyWorkflow_report_config
use rule generate_report from MyWorkflow_report as MyWorkflow_generate_report
```

# Pitfalls

## 1. 忘记传 ROOT_DIR

每个模块配置必须包含 `"ROOT_DIR": ROOT_DIR`。common.smk 通过 `config.get("ROOT_DIR", ".")` 读取此值，用于定位 `src/common/LogUtil.py` 和模块的 `bin/` 脚本。如果不传，默认值 `"."` 会导致路径错误（找不到脚本或导入失败）。

## 2. 忘记传 env

每个模块配置必须包含 `"env": config.get("env", {})`。common.smk 的 `sif()` 函数通过 `config["env"]` 查找 SIF 容器路径。如果不传，使用 `container: sif("xxx.yaml")` 的规则会在运行时报 `KeyError` 或找不到 `.sif` 文件。

## 3. 子目录模块的 conda 路径

子目录规则必须使用 `conda: "../<parent>.yaml"`（相对路径到父目录的 yaml），不能写 `conda: "../modules/<tool>/<tool>.yaml"`。

## 4. 默认参数

尽量不要在subworkflow中使用or去兜底json文件中的null值，从而提供默认值。确保json文件即是实际运行配置

## Fiber-seq 工作流

单分子染色质可及性测序分析流程。

**参考文献:**
- Stergachis et al., 2020, Science (Fiber-seq原始论文, DOI: 10.1126/science.aaz1646)
- Jha, Bohaczuk et al., 2024, Genome Research (fibertools-rs工具, DOI: 10.1101/gr.279095.124)
- https://fiberseq.github.io/

**流程步骤:**
1. `ft predict-m6a` - 从PacBio HiFi CCS数据预测m6A修饰
2. `ft add-nucleosomes` - 添加核小体调用
3. `ft fire` - 调用FIRE调控元件
4. `pbmm2` - 比对到参考基因组
5. `ft extract` - 提取数据为BED格式
6. `ft call-peaks` - FDR peak calling（shuffled null）；`ft fire --extract` - per-fiber FIRE BED
7. 下游分析（fiberseq_analysis 模块）- 注释/motif/footprint/核小体定位/定量矩阵/组间差异/共活化/单倍型/着丝粒专题/轨道
8. enrichment（GO/KEGG）、ldsc（遗传力分区，接口）

**输入:**
- PacBio HiFi CCS BAM (带kinetics标签)
- 参考基因组FASTA

**输出:**
- Fiber-seq BAM (含m6A、核小体、MSP、FIRE注释)
- BED12文件 (m6A.bed.gz, nuc.bed.gz, msp.bed.gz, fire.bed.gz)

**工具:**
- fibertools-rs (ft): https://github.com/fiberseq/fibertools-rs
- pbmm2: PacBio比对器
- hiphase: 单倍型分相 (可选)
