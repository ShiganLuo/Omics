# Omics Workflow

这里是 `Omics` 相关的统一工作流入口。当前目录下的 `run.py` 会根据样本元信息和 `workflow_name` 选择对应的 Snakemake 子流程，并自动生成各流程所需的 `raw.json` 配置文件。

- README.md为人类阅读文档
- 目录名.md为agent加载的skill

snakemake version: >= 9.16.3

## 起因

做这个项目的核心原因，是希望包容更多真实分析场景里的复杂度，并且更好地掌控那些过去已经分析过的流程。

在实际科研工作中，很多流程并不是一次性跑完就结束，而是会不断迭代：样本类型会变化，参数会调整，分析分支会增加，历史流程也常常需要回溯、复现、比较和复用。如果只是追求“先跑通一次”，这些复杂度往往会被分散在脚本、手工记录和临时改动里，后续维护成本会越来越高。`Omics` 想做的是把这些复杂度尽可能收拢到一个统一、可追踪、可扩展的工作流框架里，让过往分析过的流程真正沉淀下来，而不是随着项目结束而散失。

在选型过程中，也比较过主流的流程工具：

- `Snakemake`：规则组织直观，Python 生态结合紧密，适合逐步抽象、持续演化，也更方便对已有分析流程进行细粒度改造和接管。
- `Nextflow`：云原生和大规模调度能力强，社区生态成熟，但在当前场景下，对已有流程做细致接入和日常维护时，心智负担相对更高。
- `Cromwell`：在 WDL 生态下标准化程度高，适合强调任务描述规范和平台化执行的场景，但对我这里这种需要频繁调整、快速迭代和兼容历史分析实现的工作方式来说，不够灵活。

综合比较之后，最终选择了 `Snakemake`。原因不是它在所有场景里都最好，而是它最适合这个项目当前的目标：在保持结构化和可维护性的同时，尽可能包容复杂度，并把过去分析过的流程逐步纳入一个自己能够真正掌控的体系里。

## 目录结构

- `run.py`：统一入口脚本，读取元信息、合并配置并调用 Snakemake。
- `run.sh`：当前项目里常用的一键运行示例。
- `config/`：各工作流的模板配置文件。
- `subworkflow/`：各分析流程的 Snakemake 主入口。
- `modules/`：可复用模块定义。
- `example/`：示例配置和示例流程文件。

sample_id不能包含.
## 支持的工作流

`run.py` 通过 `-w/--workflow_name` 选择工作流：

| 工作流 | 说明 | 典型输出 |
| --- | --- | --- |
| `CoCulture` | 共培养样本分析，支持多个物种 | 物种区分后的 BAM、下游统计结果 |
| `MERIP` | MeRIP-seq / m6A-seq 分析 | dedup BAM、peak 结果、IGV 可视化 |
| `RNAseq` | 常规转录组分析 | count 矩阵、TE 表达结果 |
| `CLIP` | iCLIP / CLIP-seq 分析 | 质控、比对、PureCLIP、bedGraph / bigWig、IGV 页面 |
| `Mutation` | 体细胞突变分析（tumor vs normal） | Mutect2 VCF、Spectrum 可视化 |
| `PacVar` | PacBio 长读长变异检测 | 结构变异 VCF、SNP VCF、phasing 结果、端粒长度（4种方法）、着丝粒统计 |
| `KARRseq` | Kethoxal-Assisted RNA-RNA interaction sequencing | RNA-RNA 相互作用 pairs 文件 |
| `ncRNAseq` | 非编码 RNA 分析 | ncRNA 表达量矩阵 |
| `RNA_SNP` | RNA 变异检测 | SNP/INDEL 结果 |
| `PeakCalling` | ChIP-seq / DIP-seq peak calling 分析 | trimming、bowtie2 比对、MACS3 peak 结果 |
| `QuantMS` | 定量蛋白质组学分析（TMT/LFQ/DIA） | mzTab 定量结果、MSstats 统计分析 |
| `tRNAseq` | tRNA 修饰诱导错配测序分析（mim-tRNAseq） | 覆盖度、修饰定量、CCA 分析、DESeq2 差异表达 |

### CLIP

igv模块准备:

```json
  "igv": {
        "js": "/data/pub/zhousha/Reference/igv.min.js",
        "id": "mm39",
        "name": "Mouse (GRCm39/mm39)",
        "publicPathMap": {
            "/data/pub/zhousha/": "/data/",
            "/data/pub/zhousha/Reference/": "/ref/"
        },
        "fastaURL": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa",
        "indexURL": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai",
        "cytobandURL": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/mm39.cytoBand.txt",
        "tracks": [
            {
                "name": "Gencode vM38 genes",
                "format": "gtf",
                "url": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.basic.gene_exon.sorted.gtf.gz",
                "indexUrl": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/gencode.vM38.basic.gene_exon.sorted.gtf.gz.tbi"
            },
            {
                "name": "Gencode rmsk repeats",
                "format": "gtf",
                "url": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.sorted.gtf.gz",
                "indexUrl": "/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.sorted.gtf.gz.tbi"
            }
        ]
    }
```
1. 注释轨道

- 注释轨道需要自己建立索引，以防止浏览器全量加载注释卡死（F12观察返回码是否为200,200则是没有配置索引）
```sh
tabix -g gff [gft.gz]
```
- nginx需要确保sendfile是开启状态(on)

2. publicPathMap

为了让nginx能够获取流程生成html内的内容，特意配置了这个map，以生成相对路径让nginx能够读取
![alt text](assests/img/image.png)

## scRNAseq 和空间转录组

### `scRNAseq`

标准 Scanpy 流程从 `.h5ad` 或多个样本 h5ad 开始，规则位于 `modules/scanpy/`，`subworkflow/scRNAseq.smk` 只负责编排。流程包括：

- QC：基因数、线粒体比例、归一化、log1p、高变基因
- 降维和聚类：PCA、neighbors、UMAP、Leiden
- marker 和差异表达
- 可选高级分析：DPT 拟时序、RNA velocity、LIANA 细胞通讯、infercnvpy CNV

最小配置：

```json
{
  "input_h5ad": "/path/input.h5ad",
  "Params": {
    "scanpy": {
      "advanced": {
        "trajectory": true,
        "velocity": false,
        "communication": false,
        "cnv": false
      }
    }
  }
}
```

### `spatial_transcriptomics`

空间转录组流程位于 `subworkflow/spatial_transcriptomics.smk` 和 `modules/spatial_scanpy/`，以 10x Visium 为主，同时支持带 `obsm["spatial"]` 坐标的自定义 spot×gene h5ad。

Visium 输入：

```json
{
  "visium_h5": "/path/filtered_feature_bc_matrix.h5",
  "spatial_dir": "/path/spatial"
}
```

自定义 h5ad 输入：

```json
{
  "input_h5ad": "/path/spatial.h5ad"
}
```

流程包括空间 QC、归一化、高变基因、PCA、空间聚类、空间 marker、Squidpy 空间邻域和 Moran's I 空间自相关。cell2location 依赖已在环境模板中预留，但只有在明确参考 scRNA 表达矩阵和细胞类型 signature 输入后才应启用。

## 快速开始

1. 准备输入数据。
   - 如果传入的是元信息文件，通常应包含样本、物种、测序布局等信息。
   - 如果传入的是 fastq 目录，`run.py` 会直接从目录中解析样本信息。
2. 准备运行环境：选择 conda 模式或容器（Apptainer/SIF）模式。
   - **conda 模式**：需指定 `--conda-prefix`，`run.py` 会自动添加 `--use-conda`。
   - **容器模式**：指定 `--sdm`，`run.py` 会省略 `--use-conda`，改为让 Snakemake 用 `container:` 指令解析的 SIF 镜像运行。SIF 镜像由 `EnvUtil.py` 从 conda yaml 生成（见 `src/common/util/EnvUtil.py`），`common.smk` 的 `sif()` 函数根据 yaml 文件名 stem 查找路径。
3. 使用 `run.py` 启动对应流程。

推荐示例（conda 模式）：

```bash
python workflow/Omics/run.py \
  -m data/meta/fastq \
  -w CLIP \
  -o output \
  -t 48 \
  --log log/CLIP.log \
  --conda-prefix /data/pub/zhousha/env/mutation_0.1 \
  --Params.trim_galore.quality 10
```

容器模式示例：

```bash
python workflow/Omics/run.py \
  -m data/meta/fastq \
  -w CLIP \
  -o output \
  -t 48 \
  --sdm \
  --rerun-triggers mtime
```

容器模式下 `run.py` 会自动从配置 JSON 中提取路径并生成 `--singularity-args --bind ...`，
无需手动指定。如需覆盖，用 `--singularity-args`。

如果只是想检查流程而不真正执行，可加上 `--dry-run`。

### 测试模式（`--test`）

`run.py` 支持内置测试模式，用于快速检查某个工作流是否能正确解析配置、生成 `raw.json` 并完成 Snakemake dry-run。

```bash
# 测试单个工作流
python workflow/Omics/run.py \
  --test ncRNAseq \
  -o output

# 测试所有工作流
python workflow/Omics/run.py \
  --test all \
  -o output
```

测试模式会自动：
- 使用 `workflow/Omics/assests/test/meta_<workflow>.tsv`
- 从 `workflow/Omics/assests/test/data/` 生成测试路径
- 输出到 `output/test/`（或你指定的输出目录下的 `test/`）
- 开启 `--dry-run`
- 使用本地 `.conda` 前缀，避免污染共享环境

如果某个工作流没有对应的测试 meta 文件，`run.py` 会跳过该工作流并提示警告。

## 运行特定步骤

`run.py` 提供 `--forcerun` 参数,用于精确重跑某个规则的特定 job,且不会触发下游规则。

### 用法

```bash
# 重跑某个规则(所有 wildcards)
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun function_gsea

# 重跑某个规则的特定样本
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun trimming_Paired:sample_id=S1

# 同时重跑多个 target
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun trimming_Paired:sample_id=S1 trimming_Paired:sample_id=S2
```

### 原理

`--forcerun` 内部转换为 snakemake 的 `--until` + `--forcerun`:

- `--until` 将 DAG 截断到目标规则(含上游依赖),不包含下游
- `--forcerun` 强制重跑目标规则,即使输出文件已存在

### 规则名自动补全

subworkflow 中通过 `use rule ... as <WorkflowName>_...` 重命名了规则。`run.py` 会自动补前缀,用户只需写原始规则名:

| 用户输入 | 自动转换为 |
| --- | --- |
| `function_gsea` | `RNAseq_function_gsea` |
| `trimming_Paired:sample_id=S1` | `RNAseq_trimming_Paired:sample_id=S1` |
| `RNAseq_function_gsea` | `RNAseq_function_gsea`(已带前缀,不变) |
| `all` | `all`(特殊规则名,不变) |

### 查看 DAG 和调试

```bash
# 查看将要执行的 job(配合 --dry-run)
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun function_gsea --dry-run

# 透传其他 snakemake 参数
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --forcerun function_gsea \
  --snakemake-args --printshellcmds

# 生成 DAG 图
snakemake -s workflow/Omics/subworkflow/RNAseq.smk \
  --configfile output/RNAseq/raw.json \
  --cores 1 --sdm apptainer \
  --until RNAseq_function_gsea --forcerun RNAseq_function_gsea \
  --dag | dot -Tpng > dag.png
```

### 跳过指定规则

如需跳过某个规则及其下游,可通过 `--snakemake-args` 透传 snakemake 的 `--omit-from`:

```bash
python workflow/Omics/run.py \
  -m data/meta.tsv -w RNAseq -o output --sdm \
  --snakemake-args --omit-from RNAseq_generate_report
```

## `run.py` 参数说明

- `-m, --meta`：元信息文件或 fastq 目录。
- `-w, --workflow_name`：工作流名称，可选 `CoCulture`、`MERIP`、`RNAseq`、`CLIP`、`Mutation`、`PacVar`、`KARRseq`、`PeakCalling`、`tRNAseq`。
- `-o, --output_dir`：输出目录。
- `-t, --threads`：Snakemake 线程数。
- `--dry-run`：只生成计划，不执行。
- `--log`：日志文件路径。
- `--conda-prefix`：conda 包缓存目录（conda 模式必需，容器模式下省略）。
- `--sdm`：使用 Apptainer 容器后端（SIF 镜像）。设置后 `run.py` 不添加 `--use-conda`，改用 `container:` 指令解析的 SIF 运行。自动从配置 JSON 提取路径生成 `--singularity-args --bind ...`。
- `--singularity-args`：手动指定 singularity/apptainer 参数，如 `--singularity-args '--bind /path1,/path2'`（与 `--sdm` 配合使用，覆盖自动生成的 bind 列表）。
- `--rerun-trigger`：Snakemake 的重跑触发条件，默认 `input`。可选值及含义：

  | 参数 | 含义 |
  | --- | --- |
  | `code` | rule 定义代码（.smk 文件）发生变化时重跑 |
  | `input` | 输入文件内容（哈希）发生变化时重跑。**需要 `.snakemake/metadata/` 中的历史哈希记录；metadata 为空时无法对比，所有 rule 都会重跑** |
  | `mtime` | 输入文件修改时间比输出文件新时重跑 |
  | `params` | rule 的 params 发生变化时重跑 |
  | `software-env` | conda 环境发生变化时重跑 |

  不指定 `--rerun-trigger` 时，Snakemake 默认使用全部五个触发器。指定 `--rerun-trigger input` 表示**仅**检查 input 内容变化，不检查 code/mtime/params/software-env，更轻量但依赖 metadata。

- `--conda-frontend`：`conda` 或 `mamba`。
- `--snakemake-args`：透传给 Snakemake 的额外参数，放在这个标志后面，例如 `--snakemake-args --keep-going --rerun-incomplete`。

`run.py` 还支持额外参数透传给配置文件：

- `--key value`
- `--key=value`
- 嵌套字段：`--Params.trim_galore.quality 10`

## 输入约定

- 单端文件：通常识别为单个 `fq.gz` / `fastq.gz` 文件。
- 双端文件：通常识别为成对的 `*_1.fq.gz`、`*_2.fq.gz`，或 `*_R1.fq.gz`、`*_R2.fq.gz`。
- `trim_galore` 只是包装命令，实际运行时仍需要 `cutadapt`。
- 如果样本物种名包含空格，配置文件内部会统一规范化，例如 `Mus musculus` -> `Mus_musculus`。

## 输出约定

每次运行都会在 `output/<workflow_name>/` 下生成对应结果，同时写出 `raw.json` 和 `log/` 目录。

对于 `CLIP` 流程，当前还会额外生成：

- `bedtools/`：用于覆盖度和可视化的中间结果。
- `track/igv_track_iclip.html`：可直接打开的 IGV 浏览页面。

其中 `track` 模块会把 bigWig 和参考基因组资源整理成可在浏览器中访问的路径，因此如果在本机或服务器查看 IGV 页面，需要保证这些资源由 nginx 或其他静态服务正确暴露。

## 当前流程特点

- 支持多个工作流统一入口。
- 支持单端和双端测序。
- 配置通过模板 JSON 合并生成，便于复用和覆盖。
- 日志和输出目录由流程自动创建。

## 运行示例

当前仓库中已有的示例脚本可以直接参考 `run.sh`。它对应的典型执行方式是：

```bash
bash workflow/RNA-SNP/run.sh
```

如果需要手动调用 Snakemake，也可以参考 `run.py` 最终拼接出来的命令形式。

## 备注

- `modules/track/README.md` 说明了 IGV / UCSC track 的生成方式。
- `subworkflow/README.md` 说明了各子流程的职责和输入输出。
- 如果后续新增 workflow，建议同步补充：
  - `config/<workflow>.json`
  - `subworkflow/<workflow>.smk`
  - `subworkflow/<workflow>.json`
  `subworkflow/<workflow>.yaml`
  - `run.py` 中的分发逻辑
- 各软件传递参数的默认值均为软件或者适配流程的默认值
- " ".join(cmd)。cmd不能包含None
## 待做

- [x] 实际执行包装成shell，兼容HPC
- [x] 完善meta设计
- [x] 添加项目skill文档
- [x] 整合所有曾经分析过的流程
- [x] 添加json值校验模块





