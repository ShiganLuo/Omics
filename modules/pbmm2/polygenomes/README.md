# pbmm2/polygenomes — 多参考基因组 PacBio 长读长比对子模块

pbmm2(PacBio 官方 minimap2 封装)的多参考基因组(polygenomes)变体:
按 `genome` wildcard 将 PacBio 长读长 BAM 分别比对到各参考基因组,输出
`{genome}/{sample_id}` 层级的坐标排序 BAM + 索引。conda 环境共享父目录
`../pbmm2.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `bam_suffix`(config,形态 A) | `{indir}/{sample_id}/{sample_id}{bam_suffix}`,默认后缀 `.bam`:PacBio 长读长(HiFi/CCS)BAM |
| `genome.<G>.fasta` | 每个参考基因组 FASTA;`<G>` 为 genome wildcard(如 GRCh38、GRCm39) |
| `samples` | 样本 ID 列表(供 subworkflow 注册与日志组织) |

### 1.2 输出

```
{outdir}/{genome}/{sample_id}/{sample_id}.sorted.bam      # 坐标排序 BAM
{outdir}/{genome}/{sample_id}/{sample_id}.sorted.bai      # 索引
{logdir}/{sample_id}/{genome}/pbmm2_align.log             # 日志
```

中间产物 `{sample_id}.raw.bam` 及 pbmm2 `--sort` 附带的 `.raw.bam.bai`
在 samtools 二次排序、建索引后被 `rm -f` 删除。收尾规则 `pbmm2_result`
聚合 BAM + BAI。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.pbmm2` | 可执行文件本身 | `pbmm2` | 工具可执行路径 |
| `Procedure.samtools` | 可执行文件本身 | `samtools` | samtools 可执行路径 |
| `bam_suffix` | — | `.bam` | 输入 BAM 文件名后缀(布局参数) |
| `genome.<G>.fasta` | `pbmm2 align` 的参考 FASTA | — | 每基因组参考 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |
| — | `pbmm2 align --log-level` | `INFO` | 固定 `INFO`;pbmm2 默认 `WARN` 级别正常运行零输出,不开则日志只有 samtools 消息 |

规则默认资源:`threads = 16`、`mem_mb = 32000`(经 `rule_res("pbmm2_align")`
可覆盖)。pbmm2 固定附加 `--sort --num-threads N --log-level INFO`;模块随后再用
`samtools sort` 归一化排序并 `samtools index -o` 产出 `.sorted.bai`。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
pbmm2 align --sort --num-threads 16 --log-level INFO GRCh38.fa s1.bam s1.raw.bam
samtools sort -@ 16 s1.raw.bam -o output/GRCh38/s1/s1.sorted.bam
samtools index -@ 16 output/GRCh38/s1/s1.sorted.bam -o output/GRCh38/s1/s1.sorted.bai
rm -f s1.raw.bam s1.raw.bam.bai
```

## 2. 算法原理与生物学意义

pbmm2 是 Pacific Biosciences 官方的 minimap2 封装,针对 PacBio HiFi/CLR
读段固定推荐预设与 BAM I/O(默认按读段类型自动选 preset),输出与
下游 PacBio 工具链兼容的排序 BAM。polygenomes 场景下同一样本分别比对到
多个参考基因组(如人/鼠),用于跨物种样本(移植瘤、共生培养)的 reads
归属与污染排查。结果阅读:`samtools flagstat`/`idxstats` 看各基因组比对
率;同一样本两基因组比对率之比可判别物种来源;`{genome}/{sample_id}`
层级目录使多基因组结果互不覆盖。

## 3. 参考文献

1. pbmm2 —— https://github.com/PacificBiosciences/pbmm2(源码与文档)
2. minimap2:Li H. *Minimap2: pairwise alignment for nucleotide sequences*.
   Bioinformatics 2018, 34(18):3094–3100.
   DOI: 10.1093/bioinformatics/bty191
