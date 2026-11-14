# star/polygenomes — 多基因组 RNA-seq 比对模块(STAR)

STAR 的 polygenomes 变体:按基因组(`{genome}` 通配符)建索引、比对,
输出坐标排序 BAM + 索引与未比对 reads。子模块共享父目录 conda 环境
(`../star.yaml`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | FASTQ 目录;布局由 `isGenomeSubdir` 决定:`false` 时 `{indir}/{sample_id}/…`,`true` 时 `{indir}/{genome}/{sample_id}/…` |
| `fastq_sample_suffix`(config) | 可选:FASTQ 文件名中缀;命名 `{sample_id}[_{suffix}]_1.fq.gz`/`_2.fq.gz`(双端)、`{sample_id}[_{suffix}].single.fq.gz`(单端) |
| `genome_paired_samples` / `genome_single_samples`(config) | genome → 双端/单端样本列表(per-genome 样本布局);样本不在任一列表则报错 |
| `omics_type`(config) | 默认 `RNA-seq`;`scRNAseq` 时双端输入按 `[R2, R1]` 传给 STAR(R1 为 barcode/UMI) |
| `genome.references.<g>.fasta` | 每基因组参考 FASTA(`star_index` 必填) |
| `genome.references.<g>.gtf` | 可选:注释 GTF(有才传 `--sjdbGTFfile`) |
| `genome.references.<g>.star_index_dir` | 可选:已有 STAR 索引目录,存在则复用,否则用模块自建 `{outdir}/index/{genome}` |

### 1.2 输出

```
{outdir}/index/{genome}/                          # STAR 索引目录(star_index,声明为 directory 输出)
{outdir}/{genome}/{sample_id}/{sample_id}.bam         # 坐标排序 BAM
{outdir}/{genome}/{sample_id}/{sample_id}.bam.bai     # BAM 索引
{outdir}/{genome}/{sample_id}/{sample_id}.Unmapped.out.mate1  # 未比对 R1(可能为占位空文件)
{outdir}/{genome}/{sample_id}/{sample_id}.Unmapped.out.mate2  # 未比对 R2(同上)
```

STAR 的 `Aligned.sortedByCoord.out.bam` 规则内重命名为声明输出;
`outReadsUnmapped` 未开时 mate1/mate2 用 touch 占位。收尾规则
`star_result` 聚合 `{sample_id}.bam`。

### 1.3 参数

`Params.star` 下为 STAR CLI 原生参数名(同名 `--<name>` 传入),索引段:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `index.sjdbOverhang` | `--sjdbOverhang` | 100 | splice junction 数据库 overhang(需 GTF 才生效) |
| `index.genomeSAindexNbases` | `--genomeSAindexNbases` | 自动 | 显式配置优先;缺省按 STAR 建议 `min(14, log2(基因组大小)/2−1)` 自动算(小基因组必需) |

比对段(`Params.star.*`):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `twopassMode` | `--twopassMode` | 无 | 两遍比对(配置了才传) |
| `alignEndsType` | `--alignEndsType` | `Local` | 端比对模式 |
| `outFilterMismatchNoverReadLmax` | 同名 | 1.0 | 每 read 错配比例上限 |
| `outFilterMismatchNmax` | 同名 | 10 | 每 read 错配数上限 |
| `outFilterMultimapNmax` | 同名 | 10 | 最大多重比对位点数 |
| `outSAMtype` | 同名 | `BAM SortedByCoordinate` | 输出类型 |
| `winAnchorMultimapNmax` | 同名 | 50 | 锚窗口多重比对上限 |
| `genomeLoad` | 同名 | `NoSharedMemory` | 索引加载模式 |
| `limitBAMsortRAM` | 同名 | 0 | 排序内存上限(0=按可用内存) |
| `outReadsUnmapped` | 同名 | 无 | 未比对 reads 输出方式 |
| `outFilterMismatchNoverLmax` | 同名 | 0.3 | 比对长度错配比例上限 |
| `outFilterMatchNminOverLread` | 同名 | 0.66 | 最小匹配长度比例 |
| `alignSJoverhangMin` / `alignSJDBoverhangMin` | 同名 | 5 / 3 | 剪接位点 overhang |
| `alignSJstitchMismatchNmax` | 同名 | `0 -1 0 0` | 剪接位点拼接错配 |
| `chimSegmentMin` / `chimOutType` / `chimSegmentReadGapMax` | 同名 | 0 / `Junctions` / 0 | 嵌合(融合)比对参数 |
| `chimJunctionOverhangMin` / `chimScoreMin` / `chimScoreDropMax` / `chimScoreJunctionNonGTAG` / `chimScoreSeparation` | 同名 | 20 / 0 / 20 / -1 / 10 | 嵌合打分参数 |
| `outSAMstrandField` | 同名 | 无 | SAM 链字段 |
| `outSAMattributes` | 同名 | `NM` | SAM 标签 |
| `outMultimapperOrder` | 同名 | `Old_2.4` | 多重比对输出顺序 |
| `runRNGseed` | 同名 | 777 | 随机种子(结果可复现) |
| `outSAMmultNmax` | 同名 | -1 | 多重比对输出上限 |
| `soloType` / `soloCBwhitelist` / `soloBarcodeReadLength` | 同名 | 无 / 无 / 1 | STARsolo 单细胞参数 |
| `limitSjdbInsertNsj` | 同名 | 1000000 | 插入剪接位点上限 |
| `Procedure.STAR` / `Procedure.samtools` | — | `STAR` / `samtools` | 可执行路径 |
| threads / mem_mb | — | index 12/64000;align 12/32000 | `rule_res` 配置 |

直接调用示例(脱离流程):

```bash
cd out/index/GRCh38
STAR --runMode genomeGenerate --runThreadN 12 \
  --genomeDir out/index/GRCh38 --genomeFastaFiles genome.fa \
  --genomeSAindexNbases 14 \
  --sjdbGTFfile gencode.gtf --sjdbOverhang 100

STAR --runThreadN 12 --genomeDir out/index/GRCh38 \
  --readFilesCommand zcat --readFilesIn S1_1.fq.gz S1_2.fq.gz \
  --genomeLoad NoSharedMemory --outSAMtype "BAM SortedByCoordinate" \
  --outFileNamePrefix out/GRCh38/S1/S1.
mv out/GRCh38/S1/S1.Aligned.sortedByCoord.out.bam out/GRCh38/S1/S1.bam
samtools index -@ 12 out/GRCh38/S1/S1.bam
```

## 2. 算法原理与生物学意义

STAR 用后缀数组找最大精确匹配(seed)再做剪接感知延伸,一步式完成
RNA-seq 比对与剪接位点发现;`--sjdbGTFfile` 把已知剪接位点编入索引
提高灵敏度,`--genomeSAindexNbases` 对小基因组(ecDNA/质粒/线粒体)
必须下调(默认 14 会段错误或丢 reads,规则按基因组大小自动算)。
chimeric 参数组支持融合基因检出,STARsolo 参数组支持单细胞
barcode/UMI。生物学意义:RNA-seq 分析的第一步,剪接感知比对的准确度
决定表达定量与融合/可变剪接检出的可靠性;polygenomes 布局让同一样本
对多个参考基因组分别比对(如宿主/污染源)。结果读法:Log.final.out 的
唯一比对率、多重比对率、剪接位点数是标准质控。

## 3. 参考文献

1. STAR — https://github.com/alexdobin/STAR
2. Dobin A, Davis CA, Schlesinger F, et al. *STAR: ultrafast universal
   RNA-seq aligner*. Bioinformatics 29(1):15-21, 2013.
   DOI: 10.1093/bioinformatics/bts635
