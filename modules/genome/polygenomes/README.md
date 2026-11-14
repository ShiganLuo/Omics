# genome/polygenomes — 多基因组参考准备模块(chrom.sizes + smallRNA 提取)

polygenomes 参考准备:按基因组生成 `chrom.sizes`,并从 GTF 提取
smallRNA 基因 BED + 侧翼 FASTA(供 smallRNA 注释/去歧义等下游使用)。
引擎脚本共享父目录 `modules/genome/bin/extract_smallrna.py`;conda 环境
共享父目录 `../genome.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genome.<g>.fasta` | 每个基因组的参考 FASTA(`<g>` = `{genome}` 通配符取值;`chromosome_sizes` 消费) |
| `genome.gtf` / `genome.fasta`(config) | `extract_smallrna` 读取的注释 GTF 与基因组 FASTA(平铺键;JSON 模板按 `genome.<name>.gtf/fasta` 组织,接线时注意两处键形不同) |
| `chrom.sizes` | 上游 `chromosome_sizes` 产出(`extract_smallrna` 输入:`{outdir}/{genome}/chrom.sizes`) |

### 1.2 输出

```
{outdir}/{genome}/chrom.sizes                        # 染色体名/长度两列
{outdir}/genome/smallrna/smallrna_genes.bed          # smallRNA 基因 BED6+1
{outdir}/genome/smallrna/smallrna_genes_flank.fa     # ± flank 侧翼序列 FASTA
```

### 1.3 参数

config 位于 `Params` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.smallrna_types` | `--types`(空格拆分逐项传) | miRNA snRNA snoRNA rRNA misc_RNA scRNA scaRNA vaultRNA | 计入的 gene_type 集合 |
| `Params.smallrna_flank` | `--flank` | 50 | 基因两侧侧翼长度(bp) |
| `Procedure.samtools` | — | `samtools` | faidx 索引工具 |
| `Procedure.bedtools` | — | `bedtools` | 提取 FASTA 序列(脚本内部调用) |
| — | `samtools faidx` + `cut -f1,2` | 固定 | chrom.sizes 生成方式 |
| threads / mem_mb | — | 1 / 8000(每规则) | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > chrom.sizes

python modules/genome/bin/extract_smallrna.py \
  --gtf gencode.gtf --fasta genome.fa --chrom-sizes chrom.sizes \
  --outdir out/genome/smallrna --flank 50 \
  --types miRNA snRNA snoRNA rRNA misc_RNA scRNA scaRNA vaultRNA \
  --bedtools bedtools
```

## 2. 算法原理与生物学意义

- **chromosome_sizes**:`samtools faidx` 建索引后取 `.fai` 前两列
  (染色体名、长度),得到 `chrom.sizes`,是 bigWig/bedGraph 转换与
  区间运算的标准配套文件。
- **extract_smallrna**:解析 GTF 的 `gene` 行,按 `gene_type` ∈
  `smallrna_types` 挑出 smallRNA 基因,输出 BED6+1
  (chrom start end gene_id score strand gene_name);再以
  `± flank` 扩展区间经 bedtools 取基因组序列,输出侧翼 FASTA。
- **生物学意义**:miRNA/snRNA/snoRNA 等 smallRNA 基因短小、多拷贝、
  与 rRNA/tRNA 等高丰度 RNA 混杂;单独提取其注释区间与侧翼序列,供
  smallRNA 定量、去歧义比对与启动子/侧翼序列分析使用。侧翼序列可辅助
  区分同家族成员与设计特异探针。本模块为工具/脚本封装。

## 3. 参考文献

1. samtools faidx — https://www.htslib.org/(htslib/faidx 文档)。
2. bedtools — https://bedtools.readthedocs.io/
3. GENCODE GTF 格式 — https://www.gencodegenes.org/pages/data_format.html
4. 引擎脚本 `genome/bin/extract_smallrna.py` 为本仓库自研(共享父目录
   `genome/bin/`)。
