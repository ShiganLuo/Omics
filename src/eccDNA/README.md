# eccDNA 表达分析 (ecDNA expression evidence)

判断 ecDNA/eccDNA 是否被转录的证据链脚本与方案文档,以及断点/整体重复序列分析。

## 内容

| 文件 | 说明 |
|---|---|
| `METHOD.md` | 方案文档:三种证据层级、混杂控制、判读标准、参考文献;断点/整体重复分析 |
| `match_arriba_decoil.py` | 证据 1:Arriba 融合断点 与 decoil ecDNA 断点交叉 |
| `align_ecdna_junction.py` | 证据 2:ecDNA 重建序列 STAR 比对 + junction-spanning reads 统计(含剪接 motif 区分 circRNA) |
| `rna_dna_ratio.py` | 证据 3:外显子归一 RNA/线性背景比值(`--gtf` 外显子交集归一,屏蔽侧翼,`flank_empty`/`no_exonic`/`low_conf` 标记) |
| `breakpoint_families.py` | 跨样本断点共享家族聚类(TSV:家族 × 断点明细) |
| `plot_breakpoint_families.py` | 断点共享家族图(每家族一行,共享断点虚线标注) |
| `annotate_family_genes.py` | 家族基因座基因注释(GRCm39 GENCODE,full/partial) |
| `repeat_utils.py` | 重复分析公共层:rmsk GTF 索引(带 pickle 缓存)、fai 随机读 FASTA、decoil 环/链接加载 |
| `breakpoint_repeat.py` | 断点重复含量:窗口命中率 + 最近距离 + 家族组成,3 零模型置换 |
| `ecdna_repeat_enrichment.py` | 整体重复富集:class/family bp 覆盖 vs 同染色体同长度随机区间 |
| `junction_microhomology.py` | junction 微同源(断点序列级短重复)vs 随机断点对 |
| `plot_eccdna_repeats.py` | 上述三项 2×2 汇总图 |

## 快速开始

```bash
# 证据 1:Arriba fusion vs decoil 断点交叉(已实测)
python match_arriba_decoil.py \
    --decoil-dir output/PacVar/eccdna/PlaB_P20/decoil \
    --arriba-tsv output/RNAseq/fusion/GRCm39/E14-2/E14-2_passed_fusions.tsv \
    --arriba-discarded output/RNAseq/fusion/GRCm39/E14-2/E14-2_discarded_fusions.tsv \
    --out output/eccdna_expression/PlaB_P20_vs_E14-2/

# 证据 2:ecDNA junction-spanning reads(需先建 STAR 索引,见 METHOD.md)
python align_ecdna_junction.py \
    --ecdna-fasta output/PacVar/eccdna/PlaB_P20/decoil/reconstruct.ecDNA.filtered.fasta \
    --ecdna-bed output/PacVar/eccdna/PlaB_P20/decoil/reconstruct.ecDNA.filtered.bed \
    --bam <ecDNA-index STAR 比对产物>.bam \
    --out output/eccdna_expression/PlaB_P20/

# 证据 3:RNA/DNA 比值
python rna_dna_ratio.py \
    --ecdna-bed output/PacVar/eccdna/PlaB_P20/decoil/reconstruct.ecDNA.filtered.bed \
    --rna-bam output/RNAseq/common/11_featureCounts_bam/GRCm39/<sample>/<sample>.bam \
    --out output/eccdna_expression/PlaB_P20/
```

三个脚本互相独立,可单独运行;汇总判读标准见 `METHOD.md`。

## 快速开始:重复序列分析(断点 + 整体)

输入参数为 `SAMPLE=PATH` 显式指向文件(文件可位于任意位置,不依赖目录布局;
样本名仅用于输出标注):

```bash
REPEAT=workflow/Omics/src/eccDNA
DB=/data/pub/zhousha/Reference/mouse/GENCODE/GRCm39
EC=output/PacVar/eccdna
BEDS=(--ecdna-bed DMSO_P06=$EC/DMSO_P06/decoil/reconstruct.ecDNA.filtered.bed
      --ecdna-bed DMSO_P20=$EC/DMSO_P20/decoil/reconstruct.ecDNA.filtered.bed
      --ecdna-bed PlaB_P06=$EC/PlaB_P06/decoil/reconstruct.ecDNA.filtered.bed
      --ecdna-bed PlaB_P20=$EC/PlaB_P20/decoil/reconstruct.ecDNA.filtered.bed)
LINKS=(--ecdna-links DMSO_P06=$EC/DMSO_P06/decoil/reconstruct.links.ecDNA.filtered.txt
       --ecdna-links DMSO_P20=$EC/DMSO_P20/decoil/reconstruct.links.ecDNA.filtered.txt
       --ecdna-links PlaB_P06=$EC/PlaB_P06/decoil/reconstruct.links.ecDNA.filtered.txt
       --ecdna-links PlaB_P20=$EC/PlaB_P20/decoil/reconstruct.links.ecDNA.filtered.txt)

# Q1 断点重复含量(3 零模型:局部/环内部/全基因组,1000 次置换)
python $REPEAT/breakpoint_repeat.py "${BEDS[@]}" \
    --rmsk $DB/GRCm39_GENCODE_rmsk_TE.sorted.gtf.gz \
    --chrom-sizes $DB/chrom.sizes \
    --cache-dir output/eccdna_repeats/cache --out-dir output/eccdna_repeats

# Q2 整体重复富集(同染色体同长度随机区间 1000 次置换)
python $REPEAT/ecdna_repeat_enrichment.py "${BEDS[@]}" \
    --rmsk $DB/GRCm39_GENCODE_rmsk_TE.sorted.gtf.gz \
    --chrom-sizes $DB/chrom.sizes \
    --cache-dir output/eccdna_repeats/cache --out-dir output/eccdna_repeats

# Q1 补充:junction 微同源
python $REPEAT/junction_microhomology.py "${BEDS[@]}" "${LINKS[@]}" \
    --genome-fasta $DB/GRCm39.primary_assembly.genome.fa \
    --out-dir output/eccdna_repeats

# 汇总图(需 matplotlib,如 DNA conda env)
python $REPEAT/plot_eccdna_repeats.py --out-dir output/eccdna_repeats
```

三个分析脚本纯标准库(不依赖 pysam/pandas);rmsk 索引解析一次后 pickle
缓存(`--cache-dir`),重跑秒级。实测结果见
`../../output/eccdna_repeats/RESULTS.md`(要点:断点精确位置回避重复
(8/28,耗竭 p=0.044)但近邻 ERV1 富集 4x(FDR=0.016);整体重复含量≈背景、
LTR class 1.5x(FDR=0.044);junction 2 bp 微同源 5/14(p=0.002))。

依赖:pysam(证据 2/3,`RNA` conda env 已装)、pandas 可选。证据 1 已在
E14-2 vs PlaB_P20 实测跑通;证据 2/3 已在合成 BAM 与 E14-2 真实 RNA BAM 上验证。

⚠️ 实测结论:passed fusions 与 PlaB_P20 ecDNA 断点零命中;discarded fusions
(28 万条低置信)虽大量落位但属噪声,已标记 `discarded_only_noise_prone`。
表达判定请以证据 2 的 junction-spanning nonsplice reads 为准。
