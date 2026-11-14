# star/polygenomes/star_3pass — 小 RNA 三遍比对子模块(canonical + gene-specific)

对小非编码 RNA reads 做多轮 STAR 比对的 polygenomes 变体。目录含两个 stem:

- `star_3pass`:canonical 三遍比对(`star_3p_align`);
- `star_3pass_gene`:在 canonical 结果之上做基因特异局部重比对 + Tailer
  尾部修饰分析(`star_3pg_gene_specific`)。

两个 stem 共享**父级非 polygenomes 模块**的引擎脚本 `three_pass_align.py` /
`gene_specific_align.py` / `prepare_gene_inputs.py`(位于
`modules/star/star_3pass/bin/`,经 `config["ROOT_DIR"]` 拼接引用,不在本
目录内);conda 环境共享 `../../star.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

star_3pass(`fq_inputs` 动态解析,形态 A):

| 输入 | 说明 |
|---|---|
| `paired_samples` | 双端样本列表;命中则取 `{indir}/{sample_id}/{sample_id}_1.fq.gz` + `_2.fq.gz` |
| (其余样本) | 单端:`{indir}/{sample_id}/{sample_id}.single.fq.gz` |
| `genome.<G>.genome_index` | STAR 全基因组索引目录(`<G>` 为 genome wildcard) |
| `genome.<G>.smallrna_index` | STAR smallRNA 专一索引目录 |
| `genome.<G>.smallrna_bed` | smallRNA 基因注释 BED |

star_3pass_gene:

| 输入 | 说明 |
|---|---|
| `final_bam_dir` | canonical star_3pass 的最终 BAM 目录:`{final_bam_dir}/{genome}/{sample_id}/{sample_id}.bam(.bai)` |
| `genome.<G>.smallrna_bed` | 基因注释 BED(界定基因特异窗口) |
| `genome.<G>.genome_fasta` | 参考基因组 FASTA |

### 1.2 输出

star_3pass(`{outdir}/{genome}/{sample_id}/` 下):

```
{sample_id}_pass1.bam(.bai)                    # pass1 比对
{sample_id}_pass1_extract_1.fq.gz / _2.fq.gz   # pass1 提取的 reads
{sample_id}_pass2.bam(.bai)                    # pass2 比对(smallRNA 索引)
{sample_id}_pass2_unmapped_fq_1.fq.gz / _2.fq.gz  # pass2 未比对 reads
{sample_id}_pass3a.bam(.bai) / {sample_id}_pass3b.bam(.bai)  # pass3 两路比对
{sample_id}.bam(.bai)                          # 最终合并 BAM
{logdir}/{sample_id}/{genome}/star3p.log
```

star_3pass_gene(`{outdir}/{genome}/{sample_id}/` 下):

```
{sample_id}.bam / {sample_id}.bam.bai    # 基因特异重比对合并 BAM
{sample_id}_tail.csv                     # Tailer 尾部修饰结果(各基因拼接)
{logdir}/{sample_id}/{genome}/star3p_gene_specific.log
```

### 1.3 参数

公共工具路径(两 stem 相同):`Procedure.STAR`(默认 `STAR`)、
`Procedure.samtools`(`samtools`)、`Procedure.bedtools`(`bedtools`);
star_3pass_gene 另有 `Procedure.tailer`(默认 `Tailer`)。资源默认:
star_3p_align threads=12、mem_mb=32000;star_3pg_gene_specific
threads=2、mem_mb=8000(均经 `rule_res()` 可覆盖)。

各 pass 的 STAR 参数:config 键 `Params.star_3pass.passN.<key>`(gene 版为
`Params.star_3pass_gene.passN.<key>`)↔ 引擎 CLI `--passN-<key>`,由 .smk
白名单转发(pass1:outFilterMultimapNmax、outFilterMultimapScoreRange、
outFilterMismatchNoverLmax、alignIntronMin、alignMatesGapMax;pass2 追加
outFilterMismatchNoverReadLmax、alignEndsType、clip5pNbases、clip3pNbases;
pass3 追加 alignEndsType)。引擎默认值(`three_pass_align.py` /
`gene_specific_align.py` 的 `PASSn_STAR_DEFAULTS`):

| 参数 | pass1 | pass2 | pass3 | 说明 |
|---|---|---|---|---|
| `outFilterMultimapNmax` | 1000 | 1000 | 1000 | 允许的最大多重比对位点数 |
| `outFilterMultimapScoreRange` | 1 | 0 | 0 | 多重比对分数窗口 |
| `outFilterMismatchNoverLmax` | 0.2 | 0.2 | 0.025 | 错配数/比对长度上限 |
| `outFilterMismatchNoverReadLmax` | — | 0.05 | — | 错配数/read 长上限(pass2) |
| `alignIntronMin` | 9999999 | 9999999 | 9999999 | 等效关闭内含子(小 RNA 无剪接) |
| `alignMatesGapMax` | — | 500 | 500 | 双端最大间隔 |
| `alignEndsType` | Local | EndToEnd | Local | 末端比对方式 |
| `clip5pNbases` / `clip3pNbases` | — | `'20 0'`/`'0 20'`(PE)、`'20'`/`'0'`(SE,three_pass 为 `'20'`/`'0'`;gene 版默认 None) | — | 两端裁剪碱基数 |

注:值为 `None` 表示引擎默认不传给 STAR,但可由 config 指定;引擎另支持
`--passN-extra-args` 追加任意 STAR 参数。

star_3pass_gene 专属参数:

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.star_3pass_gene.ambiguous` | `--ambiguous` | `exclude` | 落在多个基因的歧义 reads:`exclude` 丢弃 / `duplicate` 复制到各基因 |
| `Params.star_3pass_gene.flank` | `--flank` | 50 | 基因窗口两侧延伸 bp |
| `Params.star_3pass_gene.index.genomeSAindexNbases` | `--index-genome-sa-index-nbases` | 3 | 小索引 genomeGenerate 参数 |
| `Params.tailer.threshold` | `--tailer-threshold` | 50 | Tailer 尾部判定阈值 |
| `Params.tailer.rev_comp` | `--tailer-rev-comp` | false | 尾部序列取反向互补 |
| `final_bam_dir` | `--input-bam` 的来源 | — | canonical BAM 目录 |

直接调用示例(引擎脚本在 `modules/star/star_3pass/bin/`):

```bash
# canonical 三遍比对
python $ROOT_DIR/modules/star/star_3pass/bin/three_pass_align.py \
  --fastq s1_1.fq.gz s1_2.fq.gz \
  --genome-index idx/GRCh38_genome --smallrna-index idx/GRCh38_smallrna \
  --smallrna-bed GRCh38_smallrna.bed \
  --pass1-bam s1_pass1.bam --pass1-fq1 s1_pass1_extract_1.fq.gz \
  --pass1-fq2 s1_pass1_extract_2.fq.gz \
  --pass2-bam s1_pass2.bam --pass2-unmapped1 s1_pass2_unmapped_fq_1.fq.gz \
  --pass2-unmapped2 s1_pass2_unmapped_fq_2.fq.gz \
  --pass3a-bam s1_pass3a.bam --pass3b-bam s1_pass3b.bam \
  --bam s1.bam --log star3p.log --threads 12 \
  --star STAR --samtools samtools --bedtools bedtools

# 基因特异重比对 + Tailer
python $ROOT_DIR/modules/star/star_3pass/bin/gene_specific_align.py \
  --input-bam s1.bam --input-bed GRCh38_smallrna.bed \
  --input-genome-fasta GRCh38.fa \
  --output-bam s1.bam --output-tail s1_tail.csv \
  --input-dir work/ --helper .../prepare_gene_inputs.py \
  --log star3p_gene_specific.log --threads 2 \
  --ambiguous exclude --flank 50 --index-genome-sa-index-nbases 3
```

## 2. 算法原理与生物学意义

小非编码 RNA reads(miRNA、piRNA、tRNA 片段等)短、多重比对严重、末端
修饰(甲基化、加尾)引入错配,单一一轮比对难以兼顾灵敏度与特异性。本模块
的三遍策略(具体拆分见 `three_pass_align.py`):

- **pass1**(全基因组索引,Local 末端,允许较大多重比对分数窗口):抓住
  与基因组同源的 reads 并提取(pass1_extract);
- **pass2**(smallRNA 专一索引,EndToEnd + 两端裁剪 20 nt):在 smallRNA
  注释尺度上精确指派,收集未比对 reads(pass2_unmapped)进入下一轮;
- **pass3**(更严错配阈值 0.025,Local 末端):对剩余 reads 收尾
  (pass3a/pass3b 两路),最终与前轮合并为 `{sample_id}.bam`。

各 pass 关闭内含子(`alignIntronMin=9999999`)以匹配小 RNA 无剪接的特性。
**gene 版**:以 `smallrna_bed` 的每个基因为单位截取 ±`flank` 窗口,做基因
特异的局部重比对(小索引 `genomeSAindexNbases=3` 适配短参考),歧义 reads
按 `ambiguous` 策略处理(`exclude` 提高特异、`duplicate` 提高灵敏度);
最后 Tailer 对 3' 端做尾部(poly(A)/尿嘧啶化等)判定,输出
`{sample_id}_tail.csv`。结果阅读:最终 BAM 的多重比对分布反映 smallRNA
家族共享序列;tail.csv 给出逐基因尾部修饰长度/类型,是 tRNA/piRNA 代谢
与稳定性的观测指标。

## 3. 参考文献

1. STAR:Dobin A, et al. *STAR: ultrafast universal RNA-seq aligner*.
   Bioinformatics 2013, 29(1):15–21.
   DOI: 10.1093/bioinformatics/bts635 —— https://github.com/alexdobin/STAR
2. 本模块引擎脚本:`modules/star/star_3pass/bin/`
   (three_pass_align.py、gene_specific_align.py、prepare_gene_inputs.py)。
