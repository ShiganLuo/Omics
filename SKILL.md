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

| 工作流 | snakefile | 说明 |
|--------|-----------|------|
| CoCulture | CoCulture.smk | 共培养分析 |
| MERIP | MERIP.smk | MeRIP-seq |
| RNAseq | RNAseq.smk | RNA-seq |
| ncRNAseq | ncRNAseq.smk | 非编码/小RNA-seq |
| CLIP | CLIP.smk | CLIP-seq |
| Mutation | Mutation.smk | 体细胞突变 |
| PacVar | PacVar.smk | PacBio 变异 |
| KARRseq | KARRseq.smk | KARR-seq |
| PeakCalling | PeakCalling.smk | ChIP-seq/DIP-seq |
| QuantMS | QuantMS.smk | 定量蛋白质组学 |
| tRNAseq | tRNAseq.smk | tRNA-seq |

每个工作流由对应的 `run<Workflow>()` 函数负责生成 `outfiles` 等专用字段。

## 3) 启动 Snakemake

`run.py` 构建的命令包含:

- `--configfile` 指向生成的 `raw.json`
- `--cores` 使用 CLI 参数
- `--use-conda` + `--conda-prefix` + `--conda-frontend`
- `--rerun-triggers` 使用 CLI 参数
- 可选 `--dry-run`
- 通过 `--snakemake-args` 透传额外参数

## 4) 测试模式 (--test)

```bash
cd /rna_seq_1/luoshg/Chipseq_20260709/workflow/Omics

# 测试所有工作流 (dry-run)
python run.py --test

# 测试单个工作流
python run.py --test ncRNAseq
python run.py --test PeakCalling

# 所有注册工作流均可测试
# PeakCalling, ncRNAseq, RNAseq, Mutation, CLIP, MERIP, KARRseq, CoCulture, tRNAseq, PacVar, QuantMS
```

测试资源位于 `assests/test/`:

```
assests/test/
  data/
    ref/           # 参考基因组 (touch 占位)
      GRCm39.fa
      GRCm39.gtf
      chrom.sizes
    fastq/         # 测试 FASTQ (touch 占位)
      TestIP1/     # PE IP 样本
      TestInput1/  # PE Input 样本
      TestSE1/     # SE 样本
      TestSample1/ # PE 普通样本
  meta_<Workflow>.tsv  # 11 个工作流的 meta 文件
```

`--test` 模式自动:
- 从 `WORKFLOW_DISPATCH` 获取所有注册工作流（不硬编码）
- 输出到临时目录 (`tempfile.mkdtemp`)，不在源码树内
- `atexit` 注册清理函数，测试结束后自动删除输出
- 按工作流名动态查找 `assests/test/meta_<Workflow>.tsv`
- 注入 genome.fasta / genome.gtf 路径
- 使用本地 conda-prefix（避免权限问题）

# 扩展清单

新增工作流时:

1) 在 [workflow/Omics/config](Omics/config) 添加模型模板
2) 在 `run.py` 添加 `run<Workflow>()` 函数
3) 将工作流名称加入 `--workflow_name` 选项
4) 在 `workflow/Omics/subworkflow/` 提供对应 snakefile
5) 更新 `outfiles` 逻辑以匹配预期输出
6) 更新subworkflow/README.md和README.md
7) 在 `test/` 添加对应的 meta_<Workflow>.tsv 和测试数据
8) 在 `TEST_WORKFLOWS` 字典中注册
