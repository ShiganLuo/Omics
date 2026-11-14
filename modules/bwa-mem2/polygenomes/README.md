# bwa-mem2/polygenomes — 多参考基因组短读长比对子模块

bwa-mem2(短读长比对)的多参考基因组(polygenomes)变体:按 `genome`
wildcard 对每样本 × 每参考基因组建索引并比对,输出 `{genome}/{sample_id}`
层级的坐标排序 BAM。conda 环境共享父目录 `../bwa-mem2.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `paired_samples` / `single_samples`(config) | 样本 ID 列表,决定每个样本按双端还是单端解析 FASTQ;两表都未命中时报 `FileNotFoundError` |
| FASTQ(形态 A:indir + 布局) | 双端:`{indir}/{sample_id}/{sample_id}_1.fq.gz` + `_2.fq.gz`;单端:`{indir}/{sample_id}/{sample_id}single.fq.gz`(注意单端命名无下划线,与 `_1` 不同) |
| `genome.<G>.fasta` | 每个参考基因组的 FASTA(建索引输入);`<G>` 即 genome wildcard(如 GRCh38、GRCm39) |
| `genome.<G>.index_prefix` | 可选:预建 bwa-mem2 索引前缀;若 `{prefix}.0123` 存在则直接复用,否则回退到模块自建索引路径 |

### 1.2 输出

```
{outdir}/index/{genome}/{genome}.{0123,amb,ann,bwt.2bit.64,pac}   # bwaMem2_index 产出的索引
{outdir}/{genome}/{sample_id}/{sample_id}.bam                      # 坐标排序 BAM(bwaMem2_alignment)
{logdir}/{index}/{genome}/bwa_mem2-index.log
{logdir}/{sample_id}/{genome}/bwa_mem2-alignment.txt
```

`samtools index` 还会在 BAM 旁生成 `.bam.bai`(副产物,未声明为 rule
output);收尾规则 `bwaMemm2_result` 聚合各 BAM。多参考基因组层级
`{genome}/{sample_id}` 的路径知识由 subworkflow 消费侧按同一约定拼接。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.bwaMem2` | 可执行文件本身 | `bwa-mem2` | 工具可执行路径 |
| `Procedure.samtools` | 可执行文件本身 | `samtools` | samtools 可执行路径 |
| `genome.<G>.fasta` | `bwa-mem2 index` 的输入 FASTA | — | 每基因组参考 FASTA |
| `genome.<G>.index_prefix` | `bwa-mem2 mem` 的索引前缀 | 自建 `{outdir}/index/{genome}/{genome}` | 预建索引前缀,存在 `{prefix}.0123` 时复用 |
| `paired_samples` | — | — | 双端样本列表 |
| `single_samples` | — | — | 单端样本列表 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

规则默认资源:`bwaMem2_index` threads=15、mem_mb=64000;
`bwaMem2_alignment` threads=15、mem_mb=32000(均经 `rule_res()` 可被流程
资源配置覆盖)。比对固定附加 `-T 0`(不按最低比对分数过滤,输出全部命中)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
# 建索引
bwa-mem2 index -p output/index/GRCh38/GRCh38 path/to/GRCh38.fa
# 比对(双端)→ 转 BAM → 排序 → 建索引
bwa-mem2 mem -T 0 -t 15 output/index/GRCh38/GRCh38 s1_1.fq.gz s1_2.fq.gz \
  | samtools view -b - \
  | samtools sort -@ 15 > output/GRCh38/s1/s1.bam
samtools index -@ 15 output/GRCh38/s1/s1.bam
```

## 2. 算法原理与生物学意义

bwa-mem2 是 BWA-MEM 的高速重实现(SIMD/多线程优化,算法同 BWA-MEM:
seeding → 链接(chain)→ 精细比对),对短读长(Illumina 等)做带缺口的
参考基因组比对。polygenomes 场景下同一样本分别比对到多个参考基因组
(如人 GRCh38 / 鼠 GRCm39),用于跨物种样本(移植瘤、人鼠共培养、污染
溯源)的 reads 归属判定:同一 read 在两个基因组上的比对质量差异决定它来
自哪个物种。`-T 0` 保留低分命中以避免跨物种 reads 被静默丢弃。结果阅读:
输出 BAM 按坐标排序、可直接建索引;`samtools flagstat`/`idxstats` 检查
比对率;两基因组 BAM 间比较可判别物种来源。

## 3. 参考文献

1. bwa-mem2 —— https://github.com/bwa-mem2/bwa-mem2(源码与文档;
   BWA-MEM 算法见 https://github.com/lh3/bwa)
2. samtools —— https://www.htslib.org/
