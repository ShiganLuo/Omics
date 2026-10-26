# 示例结果（真实数据渲染验证）

| 文件 | 图 | 数据 |
|------|----|------|
| `example_manhattan.png` | 全基因组 TE 位点散点（ciTotiSC vs mESC） | **真实**：GRCm39 rmsk TE 3,799,629 位点坐标 + 11 个真实 BAM 的位点级计数（featureCounts）算出的 log2FC；7,176 个位点 log2FC>1.5 |
| `example_family_box.png` | 家族表达分布（violin + Mann-Whitney p） | **真实** GRCm39_TEcount.tsv（TLSC2 按约定剔除） |
| `example_locus_coverage.png` | MT2_Mm 位点覆盖度轨道 | **真实** BAM（chr1:3146900-3147900；TLSC1 强覆盖、mESC/mFibroblast 沉默）+ 真实 rmsk TE 注释 |

`samples_mouse.tsv` 为示例用 sample→group 注释（order 列定组序）。

## 复现命令

```bash
cd /data/pub/zhousha/Totipotent20251031/workflow/Omics/src/TE
PY=/data/pub/zhousha/miniforge3/envs/DNA/bin/python
D=/data/pub/zhousha/Totipotent20251031/output/pluripotent2totipotent/RNAseq
B=$D/common/3_TEtranscripts_bam/GRCm39
TEGTF=/data/pub/zhousha/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.gtf

# 1) 家族分布（真实数据）
$PY te_family_boxplot.py \
    --matrix $D/counts/TEtranscripts/GRCm39/GRCm39_TEcount.tsv \
    --annotation examples/samples_mouse.tsv \
    --family MT2 --family L1Md --family IAP --family B2 \
    --kind violin --ncols 2 --test mannwhitney --test-ref-group mESC \
    --ylabel "log2(CPM+1)" -o examples/out -f png

# 2) 位点覆盖度（真实数据）
$PY te_locus_coverage.py \
    --track "mESC1=$B/mESC1/mESC1.bam" \
    --track "ciTotiSC1-1=$B/ciTotiSC1-1/ciTotiSC1-1.bam" \
    --track "TLSC1=$B/TLSC1/TLSC1.bam" \
    --track "mFibroblasts1=$B/mFibroblasts1/mFibroblasts1.bam" \
    --annotation examples/samples_mouse.tsv \
    --region "MT2_Mm_chr1=chr1:3146900-3147900" \
    --features-gtf $TEGTF --bin-size 25 -o examples/out -f png

# 3) Manhattan（真实数据）：位点级计数 -> log2FC -> 绘图
# 3a. featureCounts 位点级计数（SIF 需 bind /data/pub；BAM 为双端）
apptainer exec -B /data/pub:/data/pub \
    /data/pub/zhousha/Database/env/featureCounts/featureCounts.sif featureCounts \
    -T 10 -p --countReadPairs -t exon -g transcript_id -a $TEGTF \
    -o te_locus_counts.txt $B/mESC/mESC.bam $B/mESC1/mESC1.bam $B/mESC2/mESC2.bam \
    $B/ciTotiSC1-1/ciTotiSC1-1.bam $B/ciTotiSC1-2/ciTotiSC1-2.bam \
    $B/ciTotiSC2-1/ciTotiSC2-1.bam $B/ciTotiSC2-2/ciTotiSC2-2.bam \
    $B/ciTotiSC3-1/ciTotiSC3-1.bam $B/ciTotiSC3-2/ciTotiSC3-2.bam \
    $B/ciTotiSC4-1/ciTotiSC4-1.bam $B/ciTotiSC4-2/ciTotiSC4-2.bam
# 3b. 每位点 log2FC（组均值 CPM）+ Mann-Whitney p / BH padj
$PY examples/featurecounts_to_values.py --counts te_locus_counts.txt \
    --group-a ciTotiSC1-1,ciTotiSC1-2,ciTotiSC2-1,ciTotiSC2-2,ciTotiSC3-1,ciTotiSC3-2,ciTotiSC4-1,ciTotiSC4-2 \
    --group-b mESC,mESC1,mESC2 -o locus_values.tsv
# 3c. 绘图（keyed 模式：transcript_id 经 TE GTF 解析坐标）
$PY te_manhattan.py --values locus_values.tsv --te-gtf $TEGTF \
    --locus-col transcript_id --value-col value --signif-col padj \
    --value-thresh 1.5 --burst-sign up \
    --highlight-family MT2 --highlight-family IAP \
    --label-top 8 --hline 1.5 --ylabel "log2FC (ciTotiSC vs mESC)" \
    --title "TE locus burst: ciTotiSC vs mESC" --figsize 16,5.5 -o examples/out -f png
```

## 说明

- 覆盖度示例看点：MT2_Mm（MERVL）位点在 mESC 沉默、ciTotiSC 微弱、
  TLSC1 强转录，即"TE 爆发"的位点级证据形态。
- 家族分布示例看点：MT2/IAP 在 ciTotiSC、TLSC 组整体上移（p 标注 vs mESC）。
- Manhattan 示例看点：爆发位点（红）在全基因组范围成片抬升。
- 统计注意：Mann-Whitney 在 n=3 vs 8 下 p 有约 0.006 下限，26 万受检位点
  过不了 BH（min padj≈0.30），故图中 burst 标注只用 log2FC 阈值；如需
  位点级显著性请用 DESeq2/edgeR 等借力模型。
