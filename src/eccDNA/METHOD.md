# ecDNA/eccDNA 表达判定方案

目标:回答"某个 ecDNA 是否被转录"。核心原则 — **覆盖度富集不是证据,
junction-spanning 转录本才是**。拷贝数、比对偏向、DNA 污染都会造成覆盖假象;
只有跨环形连接点的 RNA reads 才能区分"存在"与"表达"。

| 证据 | 回答的问题 | 证明力 | 脚本 |
|---|---|---|---|
| 1. Arriba × decoil 断点交叉 | 环断点处有无融合转录本 | 中(间接) | `match_arriba_decoil.py` |
| 2. RNA 比对 ecDNA 重建序列 | 是否存在环化特异 junction reads | seam 决定性,fj 中,bs 阴性对照 | `align_ecdna_junction.py` |
| 3. RNA/线性背景比值 | 转录水平相对拷贝数高低 | 弱(辅助) | `rna_dna_ratio.py` |

方法归属:证据 2 为本项目自行设计(判定逻辑自 DNA 层 circle-supporting
chimeric reads,文献 4、6,迁移至 RNA 层);证据 1 对应 Zhu 2023 的 Arriba
流程;证据 3 原法为覆盖度相关分析(文献 1),"比值 vs 屏蔽侧翼"为本项目改进。

---

## 证据 1:Arriba × decoil 断点交叉

**做法**:decoil `filtered.bed` 的 fragment 边界坐标 × Arriba
`*_passed/discarded_fusions.tsv` 的断点坐标,断点落在边界 ±1 kb 内即命中。
**定位**:廉价初筛,复用主流程产物。**局限**:线性重排(易位/缺失)产生
完全相同的 read 模式,命中只能作线索,必须回到证据 2 确认。reported 断点用
fragment 边界坐标;discarded fusions 的命中标为 `discarded_only_noise_prone`,
不单独作为证据。

## 证据 2:RNA 比对 ecDNA 重建序列(核心,本项目设计)

**原理**:环化点在环分子上物理连续,以环为模板的转录物必然包含跨环化点的
连续序列;线性转录本无法产生这种比对。**ring seam**(线性化序列首尾连接点)
只存在于环上,是环化转录的决定性证据;**fragment junction**(环内相邻
fragment 连接)与融合基因断点在序列层面不可区分,只作结构支持。执行必须
分层报告。

**流程**(子集重比对,已验证;全转录组直接比对小参考需 5–13 h/样本且无增益):

```bash
# 1. 取子集(junction reads 必然来自 ecDNA 区域,子集化无损失)
samtools view -b -L regions.bed <sample>.bam > subset.bam
# 2. 还原 fastq(⚠️ 四个文件必须全部喂入 STAR:-0 只写 READ_OTHER,
#    丢弃 r1/r2 会静默漏检)
samtools fastq -1 r1.fq.gz -2 r2.fq.gz -0 r0.fq.gz -s rs.fq.gz subset.bam
# 3. 建索引(genomeSAindexNbases 按参考大小自动计算)
STAR --runMode genomeGenerate --genomeDir ecdna_index \
     --genomeFastaFiles reconstruct.ecDNA.filtered.fasta \
     --genomeSAindexNbases <auto>
# 4. 比对(PE 按 4 个独立单端输入,解除配对约束反而提高灵敏度)
STAR --genomeDir ecdna_index --readFilesIn r1.fq.gz,r2.fq.gz,r0.fq.gz,rs.fq.gz \
     --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
     --chimSegmentMin 10 --chimOutType WithinBAM SoftClip \
     --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
     --outFilterMismatchNoverReadLmax 0.04 --seedSearchStartLmax 50
```

完成后检查 `Log.final.out` 的 input reads 非零。

**read 分类与判定**:

| 类别 | 证明力 | 判定 |
|---|---|---|
| nonsplice(ring seam) | 决定性 | seam ≥3 且重复方向一致 → **表达阳性** |
| nonsplice(fragment junction) | 支持 | fj ≥3 但 seam 不足 → 弱证据("连续转录支持,环化转录未证实") |
| backsplice(GT-AG) | 阴性对照 | backsplice-only → 疑似 circRNA,不算 ecDNA 表达 |
| internal | 流程质控 | 为零 = 流程失败,重跑 |

**seam 统计 = 主脚本 pos=0 位点 + soft-clip 核查的并集**(缺一不可):STAR
local 模式常把跨 seam read 的另一半 soft-clip 掉,主脚本统计不到。核查标准:
clip ≥10 bp 与环另一侧精确匹配(容 0–1 mismatch),MAPQ ≥10。

**环间比较**:不同环 junction 位点数不同(1/2/3 片段环 = 1/2/3 个位点),
raw count 不能直接跨环比较;per-junction 归一用于排序,seam 每环恒 1 个天然
归一(判定阈值放在 seam 层的原因)。极短 fragment(<15 kb)侧的位点产出低,
跨环比较时注明。两重复不一致时取并集计数并注明差异;单重复达阈值且另一重复
方向一致可判阳性。

## 证据 3:RNA/线性背景比值

**做法**:对每个 fragment,取其左右各 ±50 kb 的基因组窗口,减去该样本全部
ecDNA 区间(屏蔽),得到线性背景;外显子区间取 GENCODE(GRCm39 vM38)全部
exon 合并。分别统计 fragment 与背景窗口内 MAPQ≥10 的 RNA reads(按
query_name 去重),并计算各自与外显子区间的交集长度,然后

```
ratio = (片段 reads / 片段∩外显子 kb) / (背景 reads / 背景∩外显子 kb)
```

分母用外显子交集长度而非区间总长,因为 RNA-seq 信号只落在外显子上,裸
区间 RPK 会把两侧外显子密度差(实测 fragment 0–14.4%、背景 1.0–10.9%)混入
比值。标记:背景零 reads → `flank_empty`;任一侧外显子交集为 0 →
`no_exonic`(环携带非编码/基因间区内容);背景外显子交集 <5 kb →
`low_conf`(泊松噪声主导,不参与强弱排序);否则按 ≥2x / 1–2x / <1x 记
`rna_above/equal/below_flank`。

**边界**:比值是相对量,只比较 fragment/环之间;reads 无法区分来自环分子
还是线性同座,本证据不能单独判定环化转录,用于刻画区段相对转录水平并与
证据 2 交叉印证。

## 断点共享家族(环间同源判定)

**判定**:两环 fragment 边界 ±200 bp 内记为一个共享断点,共享 ≥2 个断点归
入同一家族。脚本 `breakpoint_families.py`(聚类+断点表 TSV)、
`plot_breakpoint_families.py`(图,`--genes-tsv` 可加基因模型)、
`annotate_family_genes.py`(基因注释,full/partial)。数据源:
`reconstruct.ecDNA.filtered.bed`。

**实测(4 样本 × 26 环 → 9 家族)**:

| 家族 | locus | 成员(样本:环) | 结构 |
|---|---|---|---|
| 1 | chr17:27.455/27.639Mb | 104(P06), 81(P20), 99(P06), 81(P20) | 单片段,断点 0 bp |
| 2 | chr17:30.666/31.143Mb | 105(P06), 82(P20), 100(P06), 82(P20) | 单片段,1 bp |
| 3 | chr3:106.030/106.219Mb | 24(P06), 16(P20), 21(P06), 18(P20) | 单片段,0 bp |
| 4 | **chr4:145.156–145.323Mb(同源环)** | 33(P06), 22(P20), 35(P06), 31(P20) | 6 断点部分共享,结构分化 |
| 5 | chr9:46.782/47.071Mb | 66(P06), 55(P20), 69(P06), 58(P20) | 单片段,≤1 bp |
| 6 | chr5:148.975/149.316Mb | 47(P06), 43(P20) | 单片段,1 bp |
| 7 | chr1:143.415/143.555Mb | 9(P06), 9(P20) | 单片段,0 bp |
| 8 | chr4:146.27–146.57Mb | 32(P20) | 单样本环 |
| 9 | chr17:36.370/36.503Mb | 83(P20) | 单样本环(PlaB_P20:83 被 decoil 过滤) |

家族 4 的断点分层:外端 145,156/145,322 为 33+35+31 共享(除 22 的对应端);
中央 145,295/145,300 为 33+35+31 共享;145,276/277 对为 35+31 独有。

**结论口径**:家族 1/2/3/5 是"同一环的跨样本拷贝"(单片段,断点 0–1 bp,
无结构分化);家族 4 是唯一结构分化的家族(同一祖先环重排产生的结构变体)。
"同源环"指共享断点的环家族;fragment 级逐段比较只在结构相同的环间成立
(当前仅 35 vs 31)。家族基因座均有明确基因内容(家族 9 = H2-T MHC-I 簇,
家族 3 = Chil3/4 几丁质酶簇,家族 6 含 Hmgb1/Alox5ap,家族 4 = Platr16/
Gm13230/Zfp990/Zfp268),可与证据 3 转录比值直接对应。

## 决策树

```
证据1 命中(融合断点≈环化断点)?
├─ 是 → 证据2 seam reads 支持该断点?
│        ├─ 是 → 表达阳性(融合转录本 + junction counts)
│        └─ 否 → backsplice 高 → 疑似 circRNA
└─ 否 → 证据2 任意环 seam ≥3?
         ├─ 是 → ecDNA 表达(无融合)
         └─ 否 → 证据3 ratio > 2x → 弱表达/待验证(建议长读长)
                                 否则 → 无表达证据
```

## 升级路径(现有数据不足时)

- **长读长 RNA 测序**(Iso-Seq/ONT):全长跨连接点转录本,最终证明;样本选高拷贝组。
- **capture RNA-seq**:环区域+连接点两侧探针富集,成本低,适合环少的场景。
- **单细胞平行测序**(scEC&T-seq):环与转录组同细胞关联,回答异质性。

## 实测结果汇总

**证据 2**(SE,2026-10-22;4 样本 × 2 重复,26 环全覆盖):

| 样本 | 环 | fj nonsplice | seam | 判定 |
|---|---|---|---|---|
| DMSO_P06 rep1/rep2 | circ33 | 2 / **8** | 0 / **4** | **表达阳性**(rep2 达阈值,rep1 方向一致) |
| DMSO_P06 rep2 | circ105 | — | 1 | 单 read 辅助 |
| DMSO_P20 | 全部 7 环 | 0 | 0 | 无表达 |
| PlaB_P06 rep1/rep2 | circ35 | 6/5 | 1/0 | 弱证据 |
| PlaB_P20 rep1/rep2 | circ31 | 8/5 | 1/1 | 弱证据 |
| 其余 22 环 | — | 0 | 0 | 无表达(internal 非零,真实阴性) |

详见 `../../output/iteration/eccDNA_index/junction_results/RESULTS.md`。

**证据 3**(详见)
`../../output/iteration/eccDNA_index/rna_dna_ratio/RESULTS.md`):

(2026-10-24 外显子归一口径:分子分母均以 GENCODE exon 交集长度为分母,
旧区间口径的数字已废弃)

- family 6(chr5,Hmgb1/Alox5ap):242–419x,全图最强;family 2(chr17:30Mb,
  Glo1/Dnah8):5.7–19.3x,四样本一致活跃;
- 同源环(家族 4):片段异质 0.2–5.4x(旧口径 L 6–7x 系片段外显子密度
  3 倍于背景所致,已修正);family 3(Chil3/Chil4)归一后 0.3–2.5x,
  旧"6–8x 活跃"为外显子密度假象;
- 沉默:family 1(Mir 簇)0.0–0.3x(零/极低 reads,口径无关的稳健结论)、
  family 5 0.7–1.0x、family 7 0.9–1.3x、family 9 0.9x;
  family 8 两片段无外显子交集(no_exonic),RNA 趋零 — 有基因装载 ≠ 转录,
  非编码环内容可直接识别;
- **高 RNA 丰度 ≠ 环化转录**:family 2/6 RNA 最强但 junction reads 为 0;
  唯一同时具备转录富集 + seam 证据的是 DMSO_P06 circ33。

## 参考文献

1. Zhu Y, et al. **Parallel sequencing of extrachromosomal circular DNAs and
   transcriptomes in single cancer cells**. *Nat Genet* 55, 880–890 (2023).
   doi:10.1038/s41588-023-01386-y — scEC&T-seq;证据 1 参照、证据 3 原法。
2. Turner KM, et al. **Extrachromosomal oncogene amplification drives tumour
   evolution and genetic heterogeneity**. *Nature* 543, 122–125 (2017).
   doi:10.1038/nature21356 — ecDNA 基因表达与融合转录本。
3. Gao L, et al. **Integrated CIRCLE-seq with RNA-seq...** *Front Oncol* 15,
   1701989 (2025). doi:10.3389/fonc.2025.1701989 — 证据 3 比值思路。
4. Wang Y, et al. **scCircle-seq unveils the diversity...** *Nat Commun* 15,
   3574 (2024). doi:10.41467-024-45972-y — circle-supporting chimeric reads
   判定逻辑(证据 2 的源头)。
5. Chen X, et al. **A comparative analysis of eccDNA and circRNA tools...**
   (2026) — eccDNA/circRNA 工具混杂与互补。
6. Møller HD, et al. **Circle-Map**. *Nucleic Acids Res* 46, e83 (2018).
   doi:10.1093/nar/gky352 — chimeric reads 判定统计学基础。
7. Picard C, et al. **decoil** — ecDNA fasta/bed/links 来源工具。
8. Dobin A, et al. **STAR**. *Bioinformatics* 29, 15–21 (2013).
   doi:10.1093/bioinformatics/bts635 — chimeric 参数依据。
