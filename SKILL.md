---
name: Omics
description: Omics 工作流入口技能与子技能组织说明
---

# 目的

说明 `run.py` 如何构建工作流配置、启动 Snakemake，以及相关 skill 的组织方式。

# 管理的技能

- 子工作流规范: [subworkflow/subworkflow.md](subworkflow/subworkflow.md)
- 模块规范: [modules/modules.md](modules/modules.md)
- meta.tsv 填充指南: [assests/meta/meta_tsv_filling.md](assests/meta/meta_tsv_filling.md)
- mimseq 模块文档: [modules/mimseq/README.md](modules/mimseq/README.md)

# run.py 职责

## 1) 生成每个工作流的配置 (raw.json)

- 从 [workflow/Omics/config](Omics/config) 加载模型模板
- 将 CLI 参数合并进模板
	- 扁平参数更新顶层字段
	- 点号参数 (例如 `Params.STAR.alignEndsType=Local`) 更新嵌套字段
- 尽可能自动转换 CLI 参数类型 (bool/int/float)
- 填充运行时字段:
	- `ROOT_DIR`, `indir`, `outdir`, `logdir`
	- `paired_samples`, `single_samples`
	- `outfiles` (按工作流)
- 在输出目录写出 `raw.json`

## 2) 选择工作流与 snakefile

`run.py` 支持的工作流:

- `CoCulture` -> CoCulture.smk
- `MERIP` -> MERIP.smk
- `RNAseq` -> RNAseq.smk
- `CLIP` -> CLIP.smk
- `Mutation` -> Mutation.smk

每个工作流由对应的 `run<Workflow>()` 函数负责生成 `outfiles` 等专用字段。

## 3) 启动 Snakemake

`run.py` 构建的命令包含:

- `--configfile` 指向生成的 `raw.json`
- `--cores` 使用 CLI 参数
- `--use-conda` + `--conda-prefix` + `--conda-frontend`
- `--rerun-triggers` 使用 CLI 参数
- 可选 `--dry-run`
- 通过 `--snakemake-args` 透传额外参数

# 扩展清单

新增工作流时:

1) 在 [workflow/Omics/config](Omics/config) 添加模型模板
2) 在 `run.py` 添加 `run<Workflow>()` 函数
3) 将工作流名称加入 `--workflow_name` 选项
4) 在 `workflow/Omics/subworkflow/` 提供对应 snakefile
5) 更新 `outfiles` 逻辑以匹配预期输出
6) 更新subworkflow/README.md和README.md
