# eccDNA 表达分析 (ecDNA expression evidence)

判断 ecDNA/eccDNA 是否被转录的证据链脚本与方案文档。

## 内容

| 文件 | 说明 |
|---|---|
| `METHOD.md` | 方案文档:三种证据层级、混杂控制、判读标准、参考文献 |
| `match_arriba_decoil.py` | 证据 1:Arriba 融合断点 与 decoil ecDNA 断点交叉 |
| `align_ecdna_junction.py` | 证据 2:ecDNA 重建序列 STAR 比对 + junction-spanning reads 统计(含剪接 motif 区分 circRNA) |
| `rna_dna_ratio.py` | 证据 3:RNA/DNA 覆盖度比值(RNA-seq vs decoil fragment coverage) |

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

依赖:pysam(证据 2/3,`RNA` conda env 已装)、pandas 可选。证据 1 已在
E14-2 vs PlaB_P20 实测跑通;证据 2/3 已在合成 BAM 与 E14-2 真实 RNA BAM 上验证。

⚠️ 实测结论:passed fusions 与 PlaB_P20 ecDNA 断点零命中;discarded fusions
(28 万条低置信)虽大量落位但属噪声,已标记 `discarded_only_noise_prone`。
表达判定请以证据 2 的 junction-spanning nonsplice reads 为准。
