# fiberseq（Fiber-seq 下游分析家族模块，布局 C：父目录 + 子模块）

单工具链（python + numpy/scipy/matplotlib/pysam + UCSC bedGraphToBigWig），
所有子模块共享父目录 `fiberseq.yaml`（一个 SIF）。分析方法来源：
Cell Rep Methods 2024 [1]、Genome Res 2024 [2]、Nat Genet 2025 [3]。

## 子模块（每个 .smk 一个原子分析，独立可复用）

| 子模块 | 规则 | 输出（`<outdir>/{sample_id|contrast}/`） |
|---|---|---|
| `annotate/` | fiberseq_annotate_run | peaks_annotated.tsv / peaks_summary.tsv |
| `motif/` | fiberseq_motif_scan, fiberseq_motif_profile | motif_sites.bed / motif_profile.tsv / 谱图 |
| `footprint/` | fiberseq_footprint_run | footprint_{categories,profiles,scores}.tsv + PDF |
| `nucpos/` | fiberseq_nucpos_run | nuc_offset_{per_anchor,summary}.tsv + PDF |
| `quant/` | fiberseq_quant_run, fiberseq_quant_hap, fiberseq_quant_merge | region_quant.tsv / matrix/ |
| `diff/` | fiberseq_diff_run | diff_*.tsv + volcano + 基因表（供 modules/function） |
| `coactuation/` | fiberseq_coactuation_run | coactuation*.tsv |
| `haplotype/` | fiberseq_haplotype_diff | haplotype_diff.tsv |
| `censat/` | fiberseq_censat_run | censat_*.tsv + PDF |
| `track/` | fiberseq_track_run | *.bw + trackDb.txt |

共享脚本在 `bin/`，共享库 `src/common/util/FiberlibUtil.py`。

## 跨子模块输入约定（config 的 indir = results 根）

- 元素文件：`{indir}/{sample_id}/{sample_id}.{m6a,nuc,msp,cpg}.bed.gz`、`.fire.bed`、`.fire_peaks.bed`
- motif 位点：`{indir}/fiberseq_motif/{sample_id}/{sample_id}.motif_sites.bed`
- 定量矩阵：`{indir}/fiberseq_quant/matrix/`；hap 定量：`{indir}/fiberseq_quant/{sample_id}/`
- peak 注释：`{indir}/fiberseq_annotate/{sample_id}/`（diff 用，annot_sample 可指定）

## 环境

conda/SIF 共用 `fiberseq.yaml`（SIF：`fiberseq.sif`）。子模块 smk 以
`include: "../../common/common.smk"` 直接引入共享工具（与 modules/openms 一致），
`conda: "../fiberseq.yaml"`。

## Sources

- [1] https://doi.org/10.1016/j.crmeth.2024.100911
- [2] https://doi.org/10.1101/gr.279095.124
- [3] https://doi.org/10.1038/s41588-025-02324-w
