# genome — 参考基因组辅助文件构建模块

参考基因组的通用预处理:染色体长度表、smallRNA 基因 BED/FASTA 提取、
FASTA 索引。单模块(`genome.smk`)+ 私有脚本 `bin/extract_smallrna.py`;
同目录 `polygenomes/` 为多基因组版子模块(独立 .smk),本 README 只覆盖
`genome.smk`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genome.fasta` | 参考基因组 FASTA(`faidx` 与提取脚本共用);缺失即报错 |
| `genome.gtf` | GENCODE 风格注释 GTF(smallRNA 提取用);缺失即报错 |

### 1.2 输出

```
{outdir}/genome/chrom.sizes                        # 染色体长度表(fai 前两列)
{outdir}/genome/smallrna/smallrna_genes.bed        # smallRNA 基因 BED
{outdir}/genome/smallrna/smallrna_genes_flank.fa   # smallRNA 序列(±flank)
<genome.fasta>.fai                                 # samtools faidx 索引(写在 FASTA 旁)
```

日志:`{logdir}/genome/chromosome_sizes.log`、`extract_smallrna.log`、
`genome_index.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.fasta` | `samtools faidx` / `--fasta` | 必需 | 参考 FASTA |
| `genome.gtf` | `--gtf` | 必需 | 注释 GTF |
| `Params.smallrna_types` | `--types` | `miRNA snRNA snoRNA rRNA misc_RNA scRNA scaRNA vaultRNA` | 提取的 smallRNA 生物类型 |
| `Params.smallrna_flank` | `--flank` | `50` | 基因两侧扩增(bp) |
| `Procedure.samtools` | `samtools faidx` | `samtools` | 索引/长度表 |
| `Procedure.bedtools` | `--bedtools` | `bedtools` | 序列提取 |
| `Procedure.STAR` | — | `STAR` | json 声明,本 .smk 未消费 |
| 线程/内存 | — | `threads=1`,`mem_mb=8000` | `rule_res(...)` 可覆盖 |

环境:`conda: "genome.yaml"`(模块目录内)。

直接调用示例(脱离流程):

```bash
samtools faidx genome.fa && cut -f1,2 genome.fa.fai > chrom.sizes
python bin/extract_smallrna.py \
  --gtf gencode.annotation.gtf --fasta genome.fa --chrom-sizes chrom.sizes \
  --outdir out/genome/smallrna --flank 50 \
  --types miRNA snRNA snoRNA rRNA misc_RNA scRNA scaRNA vaultRNA \
  --bedtools bedtools
```

## 2. 算法原理与生物学意义

纯工具封装:① `chrom.sizes` 是很多工具(bigWig/bedGraph 转换、Circos、
覆盖度图)的必需输入;② smallRNA 提取脚本按 GTF 的 gene_biotype 筛选
smallRNA 类型,生成基因 BED 并用 bedtools 取 ±flank 序列,供 smallRNA
比对/定量(短读长专用参考);③ `genome_index` 对 FASTA 建 faidx(随机
访问索引)。结果怎么读:`smallrna_genes.bed` 每行一个 smallRNA 基因,
`_flank.fa` 为对应序列,`chrom.sizes` 两列(染色体、长度)。

## 3. 参考文献

1. samtools — https://www.htslib.org/(Li H et al., Bioinformatics
   25:2078–2079, 2009. DOI: 10.1093/bioinformatics/btp352)
2. bedtools — https://bedtools.readthedocs.io/(Quinlan AR, Hall IM,
   Bioinformatics 26:841–842, 2010. DOI: 10.1093/bioinformatics/btq033)
3. GENCODE — https://www.gencodegenes.org/(Frankish A et al.,
   Nucleic Acids Research 49(D1):D916–D923, 2021)
