# ecDNA/eccDNA 表达判定方案

目标:回答"某个 ecDNA 是否被转录"。核心原则 — **覆盖度富集不是证据,junction-spanning
转录本才是**。ecDNA 拷贝数高、比对偏向、DNA 污染都会造成覆盖假象;只有跨越环形连接点
的 RNA reads(或由其拼接出的融合转录本)才能区分"存在"与"表达"。

三份证据按证明力从弱到强排列,建议全部执行、交叉印证:

| 证据 | 数据来源 | 证明力 | 脚本 |
|---|---|---|---|
| 1. Arriba 断点 × decoil 断点交叉 | 已有 RNAseq + PacVar 输出 | 中(间接:同源断点的融合转录本) | `match_arriba_decoil.py` |
| 2. ecDNA 参考序列 junction reads | ecDNA fasta + RNA-seq fastq 重比对 | 强(直接:junction-spanning RNA) | `align_ecdna_junction.py` |
| 3. RNA/DNA 覆盖度比值 | RNA BAM + decoil bed coverage | 弱(辅助,排除拷贝数混杂) | `rna_dna_ratio.py` |

---

## 证据 1:Arriba × decoil 断点交叉

decoil 的 `reconstruct.ecDNA.filtered.bed` 给出 ecDNA 各 fragment 的基因组坐标,
`reconstruct.links.txt` 给出 fragment 间连接;Arriba 的 `*_passed_fusions.tsv`
(以及 `*_discarded_fusions.tsv` — 低置信也会保留断点坐标)给出 RNA 层面检测到的
融合断点 `breakpoint1/2`。

判定:Arriba 断点落在某 ecDNA fragment 边界 ±window 内(默认 1 kb,因为 RNA 剪接
断点与 DNA 环化断点不必精确重合 — 剪接可以跨越环化点附近少量序列),即认为该 ecDNA
存在融合转录本表达。

注意:
- 这是间接证据 — 融合转录本可能来自线性基因组重排而非环形分子。需要证据 2 确认
  断点两侧序列在 RNA 上**共价连续**(线性 deletion 的两端不会同时出现在同一条
  ecDNA 环上)。
- reported 断点要用 bed 的 fragment **边界**坐标,不是 fragment 内部。

## 证据 2:ecDNA 参考序列比对(核心方法)

`decoil reconstruct` 已输出环形重建序列 `reconstruct.ecDNA.filtered.fasta`
(每条 ecDNA 一条 record,fragment 首尾拼接)。把 RNA-seq reads 比对到该参考:

1. 建 STAR 索引(小型,秒级):
   ```bash
   STAR --runMode genomeGenerate --genomeDir ecdna_index \
        --genomeFastaFiles reconstruct.ecDNA.filtered.fasta \
        --sjdbOverhang 99 --genomeSAindexNbases 8   # 小参考必须调小
   ```
2. 比对(保留 chimeric 检测参数与主流程一致):
   ```bash
   STAR --genomeDir ecdna_index --readFilesIn <R1> <R2> \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
        --chimSegmentMin 10 --outFilterMultimapNmax 100
   ```
3. `align_ecdna_junction.py` 统计每条 ecDNA 的三类 reads(实现注意:STAR 对跨
   fasta 首尾连接点的 read 输出 split alignment(N)或 supplementary 对齐;脚本
   以"单条 alignment 内连续跨越"与"N 分割"两类统计,supplementary 链接的
   跨连接点 reads 需结合 `--chimSegmentMin` 输出的 CHIMERIC 区段另行核对):
   - **junction_spanning_nonsplice**:跨 fasta 首尾连接点的连续比对(无 N 剪接)。
     这是 ecDNA 模板转录的**直接证据** — 线性基因组的 reads 无法这样比对;
   - junction_spanning_backsplice:跨连接点且为 canonical GT-AG 剪接。⚠️ 这类更
     可能是 **circRNA**(剪接体 back-splicing),不能作为 ecDNA 表达证据,单独列出;
   - internal:ecDNA 内部比对(背景,用于归一化)。
4. 判定:`junction_spanning_nonsplice >= 3` 且显著高于阴性对照样本
   (不含该 ecDNA 的样本,估计背景 chimeric 率)即判定表达。
   报告时给出 raw count 和相对背景的富集倍数。

实测提示(2026-10,E14-2 vs PlaB_P20 decoil):Arriba 的
`*_passed_fusions.tsv` 与 PlaB_P20 ecDNA 断点无一处命中;所有匹配来自
`*_discarded_fusions.tsv`(28 万条 min_support 拒绝记录,几乎全是噪声 chimera)。
脚本已把这类匹配标记为 `discarded_only_noise_prone`,不能单独作为表达证据 —
这正是证据 2(junction reads)必须做的原因。

混杂控制(此方法天然排除大部分):
- circRNA 假阳性 → nonsplice/backsplice 分类;
- 线性 readthrough → 无法跨 fasta 首尾连接点;
- 剩余主要混杂是 mapping 错误,通过 `-outFilterMismatchNoverReadLmax 0.04` 严控。

## 证据 3:RNA/DNA 覆盖度比值

对每个 ecDNA fragment:
```
ratio = RNA-seq 该 fragment 的 RPKM / decoil 该 fragment 的 coverage
```
并与(同一 fragment 的)侧翼线性区域 ±50kb 比值、以及 ecDNA 阴性样本同基因座
比值比较。`reconstruct.ecDNA.filtered.bed` 的 coverage 列即 DNA 层覆盖。

用途:排除"高覆盖是拷贝数造成"的质疑;表达活跃的 ecDNA 其 RNA/DNA ratio 应
高于侧翼本底。单凭此不能定性,作为图表辅助。

---

## 决策树

```
证据1 命中(融合断点≈环化断点)?
├─ 是 → 证据2 nonsplice reads 支持该断点?
│        ├─ 是 → 强证据:ecDNA 表达,报告融合转录本 + junction counts
│        └─ 否 → 检查证据2 backsplice:若高,疑似 circRNA 而非 ecDNA 表达
└─ 否 → 证据2 任意 ecDNA nonsplice >= 3?
         ├─ 是 → ecDNA 表达(无融合,可能转录本未跨断点连接基因)
         └─ 否 → 证据3 ratio > 侧翼 2 倍以上 → 弱表达/待验证(建议长读长确认)
                                  否则 → 无表达证据
```

## 升级路径(现有数据不足时)

- **长读长 RNA 测序**(PacBio Iso-Seq / ONT direct RNA):单分子水平跨越连接点
  的全长转录本,最终证明;还能给出异构体结构。样本可选 ecDNA 高拷贝组(如 PlaB_P20)。
- **capture RNA-seq**:针对 ecDNA 区域 + 连接点两侧 500bp 设计探针,富集后短读长
  测序,成本低、灵敏度高,适合 ecDNA 数量少的场景。
- **单细胞平行测序**(scEC&T-seq):ecDNA 与转录组同细胞关联,可回答异质性。

## 参考文献

1. Zhu Y, et al. **Parallel sequencing of extrachromosomal circular DNAs and
   transcriptomes in single cancer cells**. *Nat Genet* 55, 1366–1377 (2023).
   doi:10.1038/s41588-023-01386-y — scEC&T-seq;ecDNA SV junction → fusion
   transcript 的论证范式;junction-spanning reads 计数方法。
2. Wang Y, et al. **EccDNA-mediated oncogene amplification in cancer** /
   ecDNA 转录综述见:Wu S, et al. *Circular DNA in cancer*;核心参照:
   Turner KM, et al. **Extrachromosomal oncogene amplification drives tumour
   evolution and genetic heterogeneity**. *Nature* 543, 122–125 (2017).
   doi:10.1038/nature21356 — ecDNA 上的基因表达与融合转录本 (e.g. EGFRvIII)。
3. Gao L, et al. **Integrated CIRCLE-seq with RNA-seq to decipher the quantity,
   localization, and functional features of eccDNA in AML**. *Front Oncol* 15,
   1701989 (2025). doi:10.3389/fonc.2025.1701989 — CIRCLE-seq + RNA-seq 整合,
   DNA 环化与转录水平一致性分析(证据 3 的比值思路)。
4. Wang Y, et al. **scCircle-seq unveils the diversity and complexity of
   extrachromosomal circular DNAs in single cells**. *Nat Commun* 15, 3574 (2024).
   doi:10.41467-024-45972-y — chimeric junction reads 作为 circle-supporting
   evidence 的判定逻辑(CPR + chimeric reads 模式)。
5. Chen X, et al. **A comparative analysis of eccDNA and circRNA tools shows
   increased accuracy of tool combination**. (2026) — eccDNA (DNA 层)与 circRNA
   (RNA 层)工具的混杂与互补;back-splice motif 区分依据。
6. Møller HD, et al. **Circular DNA elements of chromosomal repeat regions —
   Circle-Map**. *Nucleic Acids Res* 46, e83 (2018). doi:10.1093/nar/gky352 —
   circle-supporting chimeric reads 判定的统计学基础。
7. Decoil:Picard C, et al. **Reconstruction of ecDNA structures from
   long-read sequencing** (decoil, https://github.com/Decoil-Tools/decoil) —
   本项目 ecDNA fasta/bed/links 的来源工具。
8. Dobin A, et al. **STAR: ultrafast universal RNA-seq aligner**. *Bioinformatics*
   29, 15–21 (2013). doi:10.1093/bioinformatics/bts635 — chimeric 检测参数
   (chimSegmentMin 等)依据。

> 注:文献 5 为 2026 年 gigascience 比较研究(eccDNA ATAC-seq/Circle-seq vs
> circRNA RNase± 工具组合),引用时以正式卷期为准。
