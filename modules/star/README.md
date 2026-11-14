# star — RNA-seq 比对模块(STAR)

STAR 基因组索引 + 转录本感知比对:为参考基因组建索引,把每样本 FASTQ
(单/双端,scRNA-seq 时按 solo 约定换序)比对为坐标排序 BAM 并建索引,
附带嵌合比对(融合基因)参数组。同目录子模块:`polygenomes/`(多基因组版,
star_index/star_align/star_result)、`star_3pass/`(3' 端基因特异比对,
star_3pg_gene_specific/star_3p_align),各有独立 .smk,本 README 只覆盖
`star.smk`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `fastq_sample_suffix`(可选) | 形态 A 拼接:双端 `{indir}/{sample_id}/{sample_id}_{<suffix>_}1.fq.gz`、`_2.fq.gz`;单端 `{indir}/{sample_id}/{sample_id}[_<suffix>].single.fq.gz` |
| `paired_samples` / `single_samples` | 样本列表,决定单/双端;不在任一列表即报错 |
| `genome.fasta` | 参考基因组 FASTA;parse 期校验存在(`genome.is_fasta_intermediate=true` 时跳过校验,允许上游流程内产物) |
| `genome.gtf` | 参考注释 GTF(可选);提供时启用 `--sjdbGTFfile`/`--sjdbOverhang` 注释辅助建索引 |
| `genome.index_dir` | 可选:预建 STAR 索引目录,存在则跳过自建 `{outdir}/index` |
| `omics_type` | `RNA-seq`(默认)或 `scRNAseq`;scRNAseq 时 `get_alignment_input` 返回 `[R2, R1]`(R1 含 cell barcode+UMI、R2 为转录本,匹配 STARsolo 的 cDNA+barcode 顺序) |

### 1.2 输出

```
{outdir}/index/                                 # STAR 索引目录(directory())
{outdir}/{sample_id}/{sample_id}.bam            # 坐标排序 BAM(自 Aligned.sortedByCoord.out.bam 改名)
{outdir}/{sample_id}/{sample_id}.bam.bai        # samtools index
{outdir}/{sample_id}/{sample_id}.Unmapped.out.mate1   # 未比对 R1(STAR 未产出时 touch 空文件)
{outdir}/{sample_id}/{sample_id}.Unmapped.out.mate2   # 未比对 R2(同上)
```

日志:`{logdir_index}/star/{log_index_substring}.log`(默认子串 `star_index`)、
`{logdir}/{sample_id}/star_align.log`。

### 1.3 参数

config 位于 `Params.star` 下(`Params.star.index.*` 为建索引参数),CLI 一一对应:

| 参数 | CLI | 默认(.smk 回退值) | 说明 |
|---|---|---|---|
| `indir`/`outdir`/`logdir` | — | `input`/`output`/`log` | 输入/输出/日志目录 |
| `logdir_index` | — | `logdir` | 建索引日志目录 |
| `log_index_substring` | — | `star_index` | 建索引日志文件名片段 |
| `fastq_sample_suffix` | — | `None` | 文件名中样本与 read 号之间的后缀段 |
| `Procedure.STAR` | `STAR` | `STAR` | STAR 可执行路径 |
| `Procedure.samtools` | `samtools index` | `samtools` | 索引 BAM |
| `Params.star.index.sjdbOverhang` | `--sjdbOverhang` | `100` | 仅在提供 `genome.gtf` 时传入 |
| `Params.star.index.genomeSAindexNbases` | `--genomeSAindexNbases` | 自动 | 未配置时按 STAR 建议 `min(14, log2(G)/2 − 1)`(夹到 [2,14]),小基因组必需 |
| `Params.star.twopassMode` | `--twopassMode` | `None` | 非空才传 |
| `Params.star.alignEndsType` | `--alignEndsType` | `Local` | 末端比对模式 |
| `Params.star.outFilterMismatchNoverReadLmax` | `--outFilterMismatchNoverReadLmax` | `1.0` | 每读段错配比例上限 |
| `Params.star.outFilterMismatchNmax` | `--outFilterMismatchNmax` | `10` | 错配绝对上限 |
| `Params.star.outFilterMultimapNmax` | `--outFilterMultimapNmax` | `10` | 最大多重比对位点数 |
| `Params.star.winAnchorMultimapNmax` | `--winAnchorMultimapNmax` | `50` | 锚点窗口多重比对上限 |
| `Params.star.genomeLoad` | `--genomeLoad` | `NoSharedMemory` | 索引载入模式 |
| `Params.star.limitBAMsortRAM` | `--limitBAMsortRAM` | `0` | 排序内存上限(0 = STAR 默认) |
| `Params.star.outReadsUnmapped` | `--outReadsUnmapped` | `None` | 非空才传 |
| `Params.star.outFilterMismatchNoverLmax` | `--outFilterMismatchNoverLmax` | `0.3` | 每碱基错配比例上限 |
| `Params.star.outFilterMatchNminOverLread` | `--outFilterMatchNminOverLread` | `0.66` | 最小匹配长度占比 |
| `Params.star.alignSJoverhangMin` | `--alignSJoverhangMin` | `5` | 剪接位点最小跨距 |
| `Params.star.alignSJDBoverhangMin` | `--alignSJDBoverhangMin` | `3` | 注释剪接位点最小跨距 |
| `Params.star.chimSegmentMin` | `--chimSegmentMin` | `0` | 嵌合片段最小长度(>0 启用融合检测,arriba 需要) |
| `Params.star.chimOutType` | `--chimOutType` | `Junctions` | 嵌合输出类型 |
| `Params.star.chimJunctionOverhangMin` | `--chimJunctionOverhangMin` | `20` | 嵌合接头最小跨距 |
| `Params.star.outSAMstrandField` | `--outSAMstrandField` | `None` | 非空才传(intronMotif 等) |
| `Params.star.chimScoreMin` | `--chimScoreMin` | `0` | 嵌合得分阈值 |
| `Params.star.chimScoreDropMax` | `--chimScoreDropMax` | `20` | 嵌合得分最大落差 |
| `Params.star.chimScoreJunctionNonGTAG` | `--chimScoreJunctionNonGTAG` | `-1` | 非 GT-AG 接头罚分 |
| `Params.star.chimScoreSeparation` | `--chimScoreSeparation` | `10` | 嵌合位点得分间隔 |
| `Params.star.alignSJstitchMismatchNmax` | `--alignSJstitchMismatchNmax` | `"0 -1 0 0"` | 剪接拼接错配允许 |
| `Params.star.chimSegmentReadGapMax` | `--chimSegmentReadGapMax` | `0` | 嵌合片段读段内最大间隙 |
| `Params.star.outSAMattributes` | `--outSAMattributes` | `NM` | SAM 标签列表 |
| `Params.star.outMultimapperOrder` | `--outMultimapperOrder` | `Old_2.4` | 多重比对输出顺序 |
| `Params.star.runRNGseed` | `--runRNGseed` | `777` | 随机种子(确定性) |
| `Params.star.outSAMmultNmax` | `--outSAMmultNmax` | `-1` | 每读段输出比对数上限(-1 不限) |
| `Params.star.soloType` | `--soloType` | `None` | STARsolo 类型(如 `CB_UMI_Simple`) |
| `Params.star.soloCBwhitelist` | `--soloCBwhitelist` | `None` | cell barcode 白名单 |
| `Params.star.soloBarcodeReadLength` | `--soloBarcodeReadLength` | `1` | barcode 读段长度(0/1=不校验;json 示例给 98) |
| `Params.star.limitSjdbInsertNsj` | `--limitSjdbInsertNsj` | `1000000` | 注释剪接位点数上限 |
| 线程/内存(index) | `--runThreadN` | `threads=12`,`mem_mb=64000` | `rule_res("star_index")` 可覆盖 |
| 线程/内存(align) | `--runThreadN` | `threads=12`,`mem_mb=32000` | `rule_res("star_align")` 可覆盖 |

注:表中默认为 .smk 的 `or` 回退值;json 示例按 scRNA-seq 场景给了另一组推荐值
(`soloType=CB_UMI_Simple`、`soloBarcodeReadLength=98`、`outSAMattributes`
含 `CR CY UR UY`、`outSAMmultNmax=1`、`outMultimapperOrder=random`、
`limitSjdbInsertNsj=3000000`)。固定行为:`--readFilesCommand zcat`、
`--outSAMtype BAM SortedByCoordinate`、`--outFileNamePrefix {outdir}/{sample_id}/{sample_id}.`。
环境:`conda: "star.yaml"`(模块目录内)。

直接调用示例(等价引擎命令):

```bash
STAR --runMode genomeGenerate --runThreadN 12 \
  --genomeDir output/index --genomeFastaFiles genome.fa \
  --genomeSAindexNbases 14 --sjdbGTFfile genome.gtf --sjdbOverhang 100
STAR --runThreadN 12 --genomeDir output/index \
  --readFilesCommand zcat --genomeLoad NoSharedMemory \
  --alignEndsType Local --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix output/s1/s1. --readFilesIn s1_1.fq.gz s1_2.fq.gz
samtools index -@ 12 output/s1/s1.bam
```

## 2. 算法原理与生物学意义

纯工具封装(STAR:Spliced Transcripts Alignment to a Reference):基于
后缀数组/SA 索引的剪接感知比对器,分两阶段——种子搜索(最长公共前后缀定位
锚点)+ 合成剪接比对(动态规划跨越内含子),是 ENCODE 标准 RNA-seq 比对器。
本模块默认参数组同时为融合基因检测留口:`chim*` 参数组控制嵌合读段输出
(`chimSegmentMin > 0` 时开启),下游可接 arriba 等融合调用工具;
`solo*` 参数组开启 STARsolo 单细胞定量。`genomeSAindexNbases` 自动按基因组
大小取 `min(14, log2(G)/2 − 1)` 是关键工程细节——默认 14 会让小基因组
(病毒、支原体)索引构建失败。结果怎么读:BAM 坐标排序、NH/HI/AS 等标签
记录多重比对;`Unmapped.out.mate*` 为未比对读段;Log.final.out 中
`Uniquely mapped reads %` 是主要质控指标。

## 3. 参考文献

1. STAR — https://github.com/alexdobin/STAR ,
   https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
2. Dobin A et al. *STAR: ultrafast universal RNA-seq aligner*.
   Bioinformatics 29(1):15–21, 2013. DOI: 10.1093/bioinformatics/bts635
