# RmrRNA — rRNA 序列提取模块

从基因组 FASTA + GTF 注释中提取全部 rRNA(核糖体 RNA)特征序列,
输出 rRNA 参考 FASTA,供下游 rRNA 去除/比对步骤使用。引擎为本模块
自研脚本 `bin/extract_rRNA.py`(pyfaidx + multiprocessing)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genome.fasta`(config) | 参考基因组 FASTA;pyfaidx 随机访问,需有 `.fai` 索引 |
| `genome.gtf`(config) | 基因注释 GTF;脚本从中挑出 rRNA 特征(feature 或 gene_biotype/gene_type/transcript_biotype/transcript_type/biotype/type 为 `rrna`,或 gene/transcript 名称含 `rrna`) |

注:schema/JSON 模板另声明 `genome.rRNA_fasta`、`genome_fasta`、`gtf`
三个键名;规则实际读取 `genome.fasta` 与 `genome.gtf`,接线时以 .smk
键名为准。

### 1.2 输出

```
{outdir}/rRNA.fasta                  # rRNA 序列(60 列换行)
{outdir}/extract_rRNA.sh             # 命令回执
```

### 1.3 参数

无 `Params` 配置段;仅线程/内存资源(经 `resources.<rule>.threads/mem_mb`
配置,默认 threads=2、mem_mb=8000)。

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| — | `--fasta` | — | 基因组 FASTA(= `genome.fasta`) |
| — | `--gtf` | — | GTF 注释(= `genome.gtf`) |
| — | `--output` | — | 输出 rRNA FASTA |
| — | `--threads` | CPU 核数 | 提取进程并行数(流程内取规则 threads) |

直接调用示例(脱离流程):

```bash
python bin/extract_rRNA.py \
  --fasta genome.fa --gtf annotation.gtf \
  --output rRNA.fasta --threads 8
```

## 2. 算法原理与生物学意义

脚本两遍处理:第一遍流式扫 GTF,按判定规则(`feature` 或各 biotype/type
属性为 `rrna`,或 gene_name/transcript_name/gene_id/transcript_id 含
`rrna`)收集 rRNA 记录,以 transcript_id(缺省 gene_id)聚合区间;第二遍
用 pyfaidx 并发抽取每条记录的序列——优先拼接 exon 段(无 exon 用全部
区间),按坐标排序拼接,负链做反向互补。FASTA 头格式为
`<id>|<gene_name>|<transcript_name>|<gene_type>|<chrom>:<start>-<end>(<strand>)`。

生物学意义:rRNA(18S/28S/5.8S/5S 等)占总 RNA 绝大部分,在 RNA-seq
(尤其 totalRNA/ncRNA 建库)中会消耗绝大部分 reads。先提取 rRNA 参考
序列,下游才能把 reads 从 rRNA 上剥离/剔除,让 mRNA、lncRNA 等目标
RNA 的有效测序深度最大化。注意脚本按 GTF 注释认 rRNA,注释缺失的
rRNA 基因不会出现在输出里。

## 3. 参考文献

1. pyfaidx(FASTA 随机访问)—— https://github.com/mdshw5/pyfaidx
2. GTF 格式规范 —— https://github.com/The-Sequence-Ontology/Specifications/blob/master/gtf3.md
