# DRclass — DNA(WES somatic) × RNA-seq 联合证据分类模块

把 WES somatic 变异与 RNA-seq 变异/表达证据联合,对每个候选位点给出 D×R
证据标签(D1R1、D2R1 …),回答「DNA 层面的体细胞突变在 RNA 层面是否有
转录本支持、RNA 层面的候选是否只是伪影」。引擎为本模块自研脚本
`bin/drclass.py`(仅 Python 标准库),设计文档见 `3.2.2.3节写作构思.md`
第三节(脚本 docstring 中注明路径)。

## 1. 模块输入输出与参数

### 1.1 输入

所有输入为具体文件路径(`samples.<id>.xxx` 显式声明,无目录扫描):

| 输入 | 说明 |
|---|---|
| `samples`(config) | sample_id → 输入字典:`somatic_vcf`(必填,somatic 变异 VCF,仅取 FILTER=PASS 行)、`rna_vcf`(必填,RNA-seq 变异 VCF)、`rna_bam`(必填,RNA 比对 BAM,需有 `.bai` 索引,规则按 `rna_bam` 路径 `.bam→.bai` 推导)、`germline_vcf`(可选,自身 WES germline VCF)、`strain_vcf`(可选,VCF 列表,品系位点)、`splice_bed`(可选,BED 列表,剪接位点区)、`editing_bed`(可选,BED 列表,RNA 编辑位点)、`rna_bam_index`(可选) |
| `genome` | `fasta` + `fai_index`(必填)。兼容两种结构:扁平 `{fasta:...}`(模块直跑)或嵌套 `{default:..., references:{<g>:{fasta,fai_index}}}`(subworkflow 传入,取 `default` 或首个 reference) |
| VCF 兼容性 | 引擎用 gzip 流式读取,普通 VCF 与 bgzip VCF 均可;FORMAT 需含 AD/DP(AF 可选,缺失时由 AD/DP 派生) |
| BED 兼容性 | `splice_bed`/`editing_bed` 为 0-based BED 区间或 VCF 单点,自动识别 |

### 1.2 输出

```
{outdir}/{sample_id}/anchors.tsv                  # M4 锚表(chrom pos ref alt anchor)
{outdir}/{sample_id}/anchors.bed                  # 锚位点 BED(pileup -R 用)
{outdir}/{sample_id}/rna.pileup.tsv               # M3 RNA 定量(bcftools query 输出)
{outdir}/{sample_id}/{sample_id}.dna_anchored.tsv # O1:DNA 锚逐位点判定
{outdir}/{sample_id}/{sample_id}.rna_anchored.tsv # O2:RNA 锚逐位点判定
{outdir}/DRclass_summary.tsv                      # O3:全样本标签计数/穿透率/质控
```

### 1.3 参数

config 位于 `Params.DRclass` 下,与 `bin/drclass.py` CLI / `bcftools mpileup`
CLI 一一对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `min_af` | `--min-af` | 0.02 | RNA(DNA 锚)/DNA(RNA 锚)支持判定的最小等位基因频率 |
| `min_ad` | `--min-ad` | 4 | 支持判定(HIT-strong)的最小 ALT 深度 |
| `min_dp` | `--min-dp` | 20 | 判定所需的最小位点深度(≤ 判 NOINFO) |
| `germline_af` | `--germline-af` | 0.20 | germline VCF 计入伪影的最小 AF |
| `germline_dp` | `--germline-dp` | 10 | germline VCF 计入伪影的最小 DP |
| `splice_window` | `--splice-window` | 5 | 剪接位点伪影的窗口扩展(bp) |
| `min_mq` | `bcftools mpileup --min-MQ` | 20 | RNA pileup 最小比对质量 |
| `min_bq` | `bcftools mpileup --min-BQ` | 30 | RNA pileup 最小碱基质量 |
| `max_r4_fraction` | `summary --max-r4` | 0.30 | R4 占 RNA 检出比例超过此值时 summary 报 WARN |

`Procedure` 可覆盖可执行路径:`python3`、`bcftools`。

直接调用示例(脱离流程):

```bash
# 1) 锚表
python3 bin/drclass.py anchors \
  --somatic-vcf somatic.vcf.gz --rna-vcf rna.vcf \
  --out-tsv anchors.tsv --out-bed anchors.bed

# 2) RNA 定量(模块内由 bcftools 完成)
bcftools mpileup -f ref.fa -R anchors.bed -d 1000 \
  --min-MQ 20 --min-BQ 30 -a FORMAT/AD,FORMAT/DP -Ou -q 0 -Q 0 rna.bam \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n' > rna.pileup.tsv

# 3) 分类
python3 bin/drclass.py classify --sample S1 \
  --anchors anchors.tsv --pileup rna.pileup.tsv \
  --somatic-vcf somatic.vcf.gz --rna-vcf rna.vcf \
  --min-af 0.02 --min-ad 4 --min-dp 20 \
  --germline-af 0.20 --germline-dp 10 --splice-window 5 \
  --out-dna S1.dna_anchored.tsv --out-rna S1.rna_anchored.tsv

# 4) 汇总
python3 bin/drclass.py summary --tables S1.dna_anchored.tsv S1.rna_anchored.tsv \
  --max-r4 0.30 --out DRclass_summary.tsv
```

## 2. 算法原理与生物学意义

### 2.1 生物学意义

WES somatic 检出的突变是否真的在肿瘤 RNA 中转录?RNA-seq 变异检出的候选
有多少是剪接、RNA 编辑、品系多态或 germline 残留的伪影?本模块用同一
状态机对两个方向的证据做系统分类:DNA 锚位点看 RNA 支持(D1Rx),RNA
锚位点看 DNA 支持(DxR1),并把伪影标记为 R4/R1-ANCHOR-REJ,最终以
穿透率(有转录支持的 DNA 突变占比)与 R4 占比给出样本质控。

### 2.2 算法原理

- **锚表(M4)**:somatic VCF 中 FILTER=PASS 的位点为 DNA 锚;RNA VCF
  检出但不在 DNA 锚内的位点为 RNA 锚(`anchors.tsv` 的 `anchor` 列),
  `anchors.bed` 为 pileup 提供目标区间。
- **RNA 定量(M3)**:`bcftools mpileup -R anchors.bed -d 1000` 对 RNA BAM
  在锚位点计数,取 FORMAT/AD、FORMAT/DP,VAF = ALT_AD / DP。
- **状态机优先级**:REJECTED > HIT-pass > HIT-strong > HIT-weak > ABSENT >
  NOINFO。
- **表 1(DNA 锚,D 轴恒 HIT-pass,判 R 轴)**:
  - RNA VCF 检出且无伪影 → `D1R1`;检出但带伪影标记(germline/strain/
    splice/editing/multiallelic)→ `D1R4`(REJECTED);
  - 未检出按 pileup 判:DP>min_dp 且 AF≥min_af 且 ALT_AD>min_ad → `D1R2`
    (HIT-strong);DP>min_dp 且 AF≥min_af → `D1R3`(HIT-weak,弱支持);
    仅 DP>min_dp → `D1R6`(ABSENT);否则 `D1R5`(NOINFO)。
- **表 2(RNA 锚,R 轴恒 HIT-pass,判 D 轴,锚级伪影先剔除)**:
  - germline/strain 命中 → `D0R1`;splice/editing/multiallelic →
    `R1-ANCHOR-REJ`;
  - 其余按 somatic VCF 判:AF≥min_af 且 DP>min_dp 且 ALT_AD>min_ad →
    `D2R1`;AF≥min_af 且 DP>min_dp → `D3R1`;仅 DP>min_dp → `D6R1`;
    否则或不在 callset → `D5R1`。
- **汇总(O3)**:穿透率 = (D1R1+D1R2)/(D1R1+D1R2+D1R6);
  r4_fraction = D1R4/(D1R1+D1R4),超过 `max_r4_fraction` 报
  `WARN:r4_fraction>…`,否则 QC=OK。
- **结果怎么读**:`*.dna_anchored.tsv`/`*.rna_anchored.tsv` 每行一个锚位点,
  含 DNA/RNA 两侧的 call status、FILTER、AF/AD/DP、`artifact_flags`、
  `germline_hit_source` 与 `dr_label`;`DRclass_summary.tsv` 每样本
  一行给出 12 类标签计数、穿透率与 QC。

## 3. 参考文献

1. bcftools — https://samtools.github.io/bcftools/(mpileup/query 子命令,
   模块内 RNA 定量引擎)。
2. Danecek P et al. *Twelve years of SAMtools and BCFtools*. GigaScience
   10(2):giab008, 2021. DOI: 10.1093/gigascience/giab008。
3. 分类规则的内部设计文档:`3.2.2.3节写作构思.md` 第三节(路径见
   `bin/drclass.py` docstring)。
