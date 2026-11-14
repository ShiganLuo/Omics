# cellranger/polygenomes — 多基因组单细胞 RNA-seq 定量子模块(Cell Ranger)

polygenomes(每基因组独立目录)布局下的 10x Cell Ranger 封装:
每基因组建参考(`cellranger mkref`)→ 每样本每基因组 `cellranger count`
出基因 × 细胞矩阵 → 自研脚本转 h5ad(供 Scanpy)。共享父目录
`../cellranger.yaml` 环境与 `modules/cellranger/bin/` 脚本
(`cellranger_ref.py`、`cellranger_to_h5ad.py`)。Cell Ranger 为 10x
官方二进制,经 `Procedure.cellranger` 注入(非 conda 包)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `cellranger_input_dict` (config) | 样本 → `{fastq_dir, sample_prefix}` 字典(FASTQ 目录 + `--sample` 前缀);`fastq_dir` 不存在直接 `FileNotFoundError` |
| `genome.<g>.fasta` / `.gtf` (config) | 每基因组参考 FASTA/注释(mkref 输入) |
| `genome.<g>.genome_name` (config) | 参考名(缺省用 genome ID) |
| `genome.<g>.cellranger_ref_dir` (config) | 可选预建参考;存在则跳过 mkref,否则要求 fasta+gtf 齐全并使用本模块构建的 `{outdir}/cellranger_ref/{genome}/{genome_name}` |
| `samples` / `outdir` / `logdir` / `logdir_ref` / `h5ad_outdir` (config) | 样本列表与目录 |

### 1.2 输出

```
# cellranger_ref(每基因组)
{outdir}/cellranger_ref/{genome}/                      # mkref 参考目录(directory)

# cellranger_count(每基因组 × 样本;cellranger outs/ 移出后)
{outdir}/{genome}/{sample_id}/{sample_id}.bam                      # possorted BAM
{outdir}/{genome}/{sample_id}/raw_feature_bc_matrix/               # 原始矩阵目录
{outdir}/{genome}/{sample_id}/filtered_feature_bc_matrix/          # 过滤后矩阵目录
{outdir}/{genome}/{sample_id}/filtered_feature_bc_matrix.h5        # 过滤矩阵 h5

# cellranger_to_h5ad(每基因组 × 样本)
{h5ad_outdir}/{genome}/{sample_id}/{sample_id}_cellranger.h5ad     # Scanpy 输入
```

回执:`{outdir}/cellranger_ref/{genome}/cellranger_ref.sh`、
`{outdir}/{genome}/cellranger_count_{sample_id}.sh`、
`{h5ad_outdir}/{genome}/{sample_id}/cellranger_to_h5ad.sh`;日志
`{logdir_ref}/cellranger_ref/{genome}/...`、`{logdir}/{sample_id}/{genome}/...`。

### 1.3 参数

config 位于 `Params.cellranger`、`Procedure`、`genome` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.cellranger.mkref.version` | `cellranger_ref.py --version` | `2024-A` | Cell Ranger 参考版本标签 |
| `Params.cellranger.count.create-bam` | `--create-bam` | true | 是否输出 BAM |
| `Params.cellranger.count.nosecondary` | `--nosecondary` | true | 丢弃次要比对 |
| `genome.<g>.genome_name` | `--genome` / `--transcriptome` 名 | genome ID | 参考命名 |
| `Procedure.cellranger` | 引擎二进制 | `cellranger` | 10x 官方二进制 |
| (threads) | `--localcores` / `--nthreads` | 16(ref)/ 8(count) | `rule_res` 可覆盖 |
| (mem_mb) | — | 64000(ref)/ 32000(count) | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
# 建参考(父目录 bin/ 封装 mkref)
python modules/cellranger/bin/cellranger_ref.py \
    --fasta genome.fa --gtf genome.gtf \
    --output out/cellranger_ref/GRCh38 --genome GRCh38 \
    --version 2024-A --cellranger cellranger --nthreads 16

# 定量(在 out/{genome}/ 下执行,outs/ 产物随后移入样本目录)
cellranger count --id sample1 \
    --transcriptome out/cellranger_ref/GRCh38/GRCh38 \
    --fastqs /path/to/fastqs --sample sample1 \
    --localcores 8 --create-bam true --nosecondary

# 转 h5ad
python modules/cellranger/bin/cellranger_to_h5ad.py \
    --input out/GRCh38/sample1/filtered_feature_bc_matrix \
    --output out/3_h5ad/GRCh38/sample1/sample1_cellranger.h5ad \
    --sample-id sample1
```

注:`cellranger_count` 运行前会**删除**已存在的样本输出目录(cellranger
对非 pipestance 存量目录会失败),重跑即重建。

## 2. 算法原理与生物学意义

Cell Ranger `count` 把 10x 单细胞读段(Barcode + UMI + cDNA)比对到
转录组参考,按细胞 Barcode 聚类划分真实细胞,输出基因 × 细胞计数矩阵
(raw 保留全部 Barcode,filtered 只留判为细胞的列)。polygenomes 变体
的关键差异:同一 FASTQ 对每个基因组各建参考、各跑一次 count——适用于
跨物种(人/鼠共培养、异种移植)或双组装比较,`{genome}` 段进全部输出
路径,同一细胞在不同参考下的归属/污染(如双物种 disambiguate 场景)
可对照评估。h5ad 转换(`cellranger_to_h5ad.py`)把 filtered 矩阵包成
AnnData,供 Scanpy 聚类/注释(下游 scanpy 模块)。已知 scTE/Cell Ranger
细胞数不一致问题见父目录 `KNOWN_ISSUES.md`。

## 3. 参考文献

1. Cell Ranger(10x Genomics)—— https://www.10xgenomics.com/support/software/cell-ranger
2. 父模块 bin/ 脚本(`cellranger_ref.py` / `cellranger_to_h5ad.py`)为
   本仓库自研封装。
