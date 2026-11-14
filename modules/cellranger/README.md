# cellranger — 单细胞 RNA-seq 比对定量与 h5ad 转换模块

Cell Ranger 三规则:`cellranger_ref`(建参考)→ `cellranger_count`(比对
+ 基因-细胞表达矩阵)→ `cellranger_to_h5ad`(过滤矩阵转 AnnData)。
Cell Ranger 以二进制提供(非 conda 包),经 `Procedure.cellranger` 指定。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 规则 | 说明 |
|---|---|---|
| `genome.fasta` / `genome.gtf` | cellranger_ref | Ensembl 风格 FASTA/GTF(脚本可从 URL 下载或用本地路径) |
| `genome.cellranger_ref_dir` | cellranger_count | 预建参考目录(存在则优先,跳过 mkref) |
| `cellranger_input_dict`(config) | cellranger_count | 样本 → `{"fastq_dir": ..., "sample_prefix": ...}`(fastq 目录 + 样本名前缀过滤);目录不存在抛 `FileNotFoundError` |
| `{outdir}/{sample_id}/filtered_feature_bc_matrix` | cellranger_to_h5ad | 过滤后表达矩阵目录;`ctd_input_mode == "cellranger_dir"` 时改取 `{outdir}/{sample_id}/outs/filtered_feature_bc_matrix` |

### 1.2 输出

```
{outdir}/cellranger_ref/{genome_name}/                    # mkref 参考目录(directory())
{outdir}/{sample_id}/{sample_id}.bam                      # 比对 BAM(由 possorted_genome_bam.bam 改名)
{outdir}/{sample_id}/raw_feature_bc_matrix/               # 原始计数矩阵目录
{outdir}/{sample_id}/filtered_feature_bc_matrix/          # 过滤后计数矩阵目录
{outdir}/{sample_id}/filtered_feature_bc_matrix.h5        # HDF5 矩阵
{h5ad_outdir}/{sample_id}/{sample_id}_cellranger.h5ad     # AnnData(默认 h5ad_outdir = {outdir}/3_h5ad)
{outdir}/cellranger_ref/cellranger_ref.sh 等              # 各规则回执
logdir/{sample_id}/cellranger_count.log                   # 运行日志(ref 日志在 logdir_ref)
```

`cellranger_count` 运行前会删除残留的样本目录(cellranger 对非 pipestance
目录会报错),运行后把 `outs/` 里的产物 mv 到声明输出。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.cellranger` | 可执行文件 | null(回退 `cellranger`) | Cell Ranger 二进制 |
| `Params.cellranger.mkref.genome_name` | `--genome` | `GRCh38` | 参考名(拼输出目录 `{outdir}/cellranger_ref/{genome_name}`) |
| `Params.cellranger.mkref.version` | `--version` | `2024-A` | 参考版本标签 |
| `Params.cellranger.count.create-bam` | `--create-bam` | false | 是否产出 BAM(json 默认 true) |
| `Params.cellranger.count.nosecondary` | `--nosecondary` | true | 跳过二对比对;经 `or True` 回退,实际恒为 true |
| `cellranger_input_dict` | `--fastqs` / `--sample` | {} | 样本 → fastq 目录与样本前缀 |
| `ctd_input_mode` | — | `bam` | to_h5ad 取矩阵位置:`bam` 用整理后目录,`cellranger_dir` 用原始 outs/ |
| `h5ad_outdir` | — | `{outdir}/3_h5ad` | h5ad 输出根 |
| `logdir_ref` | — | `logs/ref` | ref 规则日志目录 |
| `indir` / `outdir` / `logdir` / `samples` | — | `input` / `output` / `logs` / [] | 基础键 |
| `genome.fasta` / `gtf` / `cellranger_ref_dir` | — | — | 见输入 |

线程/内存经 `rule_res(...)` 配置:cellranger_ref 16/64000,cellranger_count
8/32000(`--localcores`),cellranger_to_h5ad 1/8000。

直接调用示例(模块 bin/ 引擎脚本):

```bash
python modules/cellranger/bin/cellranger_ref.py \
  --fasta genome.fa --gtf genome.gtf --output out/cellranger_ref \
  --genome GRCh38 --version 2024-A --cellranger /path/to/cellranger --nthreads 16

python modules/cellranger/bin/cellranger_to_h5ad.py \
  --input out/s1/filtered_feature_bc_matrix --output out/3_h5ad/s1/s1_cellranger.h5ad --sample-id s1
```

## 2. 算法原理与生物学意义

`cellranger_ref.py`:整理 Ensembl 注释供 mkref——去掉 GTF 的 Ensembl ID
版本后缀、按 10x 标准 biotype 白名单过滤(protein_coding、lncRNA、
IG/TR 基因)、移除 Y 染色体 PAR_Y 基因,再 `cellranger mkref` 建参考
(保留 Ensembl 染色体名)。`cellranger count`:barcode 识别/纠错、UMI
归并、STAR 比对、生成 raw/filtered 基因-细胞计数矩阵(filtered 为
cellranger 细胞判别保留的 barcode)。
`cellranger_to_h5ad.py`:用 scanpy `read_10x_mtx` 读过滤矩阵(var 名取
gene symbols 并去重),补 `sample_id` 到 obs,写标准 AnnData。生物学
意义:10x 单细胞 3' 端文库按细胞 barcode 拆分基因表达,过滤矩阵是
下游聚类/注序(Scanpy)分析的起点;h5ad 是单细胞分析的事实标准容器
格式。

## 3. 参考文献

1. Cell Ranger(10x Genomics)—
   https://www.10xgenomics.com/support/software/cell-ranger
2. Scanpy — https://scanpy.readthedocs.io/ ;Wolf FA et al. *Scanpy:
   large-scale single-cell gene expression data analysis*. Genome Biology
   19:15, 2018. DOI: 10.1186/s13059-017-1382-0
