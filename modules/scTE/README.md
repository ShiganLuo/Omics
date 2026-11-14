# scTE — 单细胞 TE 表达定量模块(scTE)

两规则:`scTE_build_index`(建 TE 索引)→ `scTE_quantify`(BAM → TE
计数 h5ad)。引擎为 pip 包 `scte-quant`(CLI `scTE` / `scTE_build`),
定量包装脚本为模块自带 `bin/scTE_quantify.py`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 规则 | 说明 |
|---|---|---|
| `genome.gtf` + `genome.te_bed` | scTE_build_index | 本地基因 GTF + TE BED(两者都给时走资源构建;否则按 genome 名下载) |
| `genome.scTE_index` | scTE_quantify | 预建索引路径(存在则直接用,跳过构建) |
| `{indir}/{sample_id}/{sample_id}.bam` | scTE_quantify | Cell Ranger 比对 BAM(带 CB/UB 标签);`input_mode == "cellranger_dir"` 时取 `{indir}/{sample_id}/outs/possorted_genome_bam.bam` |

### 1.2 输出

```
{outdir}/index/{genome}.{mode}.idx            # scTE 索引(scTE_build 产物)
{outdir}/{sample_id}/{sample_id}_scTE.h5ad    # TE/基因 × 细胞计数矩阵(AnnData)
{outdir}/index/scTE_build_index.sh            # 回执(scTE_quantify.sh 在样本目录)
logdir_index/scTE/scTE_build_index.log        # 构建日志(默认 logs/index)
logdir/{sample_id}/scTE_quantify.log          # 定量日志
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.scTE_build` | 可执行文件 | `scTE_build` | 索引构建 CLI |
| `Procedure.scTE` | `--scte-bin` | `scTE` | 定量 CLI |
| `Params.scTE.genome` | `-g` | `hg38` | 基因组名(下载构建)或 `-g other`(资源构建);也拼索引文件名 |
| `Params.scTE.mode` | `-m` | `exclusive` | 索引模式(TE/基因归属策略,见 scTE 文档);拼索引文件名 |
| `Params.scTE.cb_tag` | `--cb-tag` | `CB` | 细胞 barcode 标签 |
| `Params.scTE.umi_tag` | `--umi-tag` | `UB` | UMI 标签 |
| `genome.gtf` / `genome.te_bed` | `-gene` / `-te` | — | 本地注释(资源构建) |
| `genome.scTE_index` | `--index` | — | 预建索引 |
| `input_mode` | — | `bam` | `cellranger_dir` 时 BAM 取 cellranger outs/ 原始路径 |
| `logdir_index` | — | `logs/index` | 构建日志目录 |
| `indir` / `outdir` / `logdir` / `samples` | — | `input` / `output` / `logs` / [] | 基础键 |

注意:json 把 `gene_gtf`/`te_bed` 放在 `Params.scTE.*`,而 .smk 读取
`genome.gtf`/`genome.te_bed`,接线以 .smk 为准。线程/内存经
`rule_res(...)` 配置:build 1/8000,quantify 8/16000。

直接调用示例(模块 bin/ 引擎脚本,回执原文):

```bash
scTE_build -gene genes.gtf -te te.bed -o out/index/hg38 -m exclusive -g other
python modules/scTE/bin/scTE_quantify.py \
  --input s1.bam --output out/s1/s1_scTE.h5ad --index out/index/hg38.exclusive.idx \
  --cb-tag CB --umi-tag UB --threads 8 --scte-bin scTE --sample-id s1
```

## 2. 算法原理与生物学意义

scTE 按细胞 barcode(CB)与 UMI(UB)标签把 BAM reads 归并到
TE/基因 × 细胞的计数矩阵。`bin/scTE_quantify.py` 在调 scTE 前先用
`samtools -e '[CB] && [UB]'` 预过滤 BAM(Cell Ranger BAM 含未纠错
barcode 的 reads,scTE 的标签检查要求头部 reads 均带标签),再把 scTE
输出转为稀疏 AnnData(h5ad)供 Scanpy 使用。生物学意义:单细胞水平的
TE 重激活与基因组稳定性/细胞状态相关,scTE 矩阵可与基因表达矩阵联合
降维聚类。已知问题:scTE 的 `--expect-cells`(默认 10000)是按 reads
数取 top N barcode 的硬截断,细胞集合可能与 Cell Ranger 不一致——详见
模块内 `KNOWN_ISSUES.md`。

## 3. 参考文献

1. scTE(`scte-quant` pip 包,CLI `scTE`/`scTE_build`,见 `scTE.yaml`)—
   https://pypi.org/project/scte-quant/
