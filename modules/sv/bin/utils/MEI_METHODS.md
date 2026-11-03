# MEI (Mobile Element Insertion) 分析方法文档

针对 SURVIVOR 合并 SV VCF（pbsv 来源）的转座事件分析管线。
脚本:本目录 `mei_analysis.py`(分析)+ `mei_plot.py`(可视化)。
运行环境:`sv.sif` 容器(RepeatMasker / pysam / pandas / matplotlib)。

## 一、输入契约

| 输入 | 要求 | 说明 |
| --- | --- | --- |
| `-i` merged VCF | SURVIVOR merge(pbsv 单样本 calls) | INS 的 ALT 携带完整插入序列(pbsv AAL 字段被 SURVIVOR 保留),这是分类的原始证据 |
| `--te-lib` | TE 一致序列 FASTA | 本机:`/data/16T/Database/Reference/famdb/rhesus_te_lib.fa` |
| `--gtf` | 基因注释 GTF(.gz) | 用于 gene/exon/intron 语境注释 |
| `--te-bed` | 参考 TE 区域 BED | rmsk lift-over,用于 TE-into-TE 判定 |
| `--fasta` | 已索引参考 FASTA | 用于 TSD 检测(需 .fai) |

SUPP_VEC 位数自适应:对照组在前、实验组在后的 merge 顺序下,
`10`=ctrl 特异、`01`=exp 特异、`11`=共享;多样本时逐位对应。

## 二、分析流程(五步)

1. **INS 提取**:解析 VCF,取 SVTYPE=INS 且 ALT 带序列的记录
   (REF/ALT 共享前缀裁剪后得到纯插入序列)。symbolic `<INS>` 无序列
   的记录丢弃并计数。
2. **TE 分类**:插入序列 ≥100bp 写 FASTA → RepeatMasker
   `-s -no_is -lib <te_lib>`;每条取覆盖跨度最大的命中,
   记录亚家族/覆盖度/链方向。家族归桶:
   Alu(SINE)、L1(LINE)、SVA、LTR、DNA、simple_repeat、other、unclassified。
   <100bp 标记 short_unclassified(多为测序噪声,保留不删)。
3. **MEI 指纹**(L1 介导插入的标准特征):
   - **poly(A) 尾**:末 60bp 内连续 A≥12 或 A 占比≥90%;
   - **TSD(靶位点重复)**:断点两侧紧邻的 7–20bp 完全相同 k-mer
     (参考基因组经 pysam 取侧翼)。注:参考基因组该位点若已带祖先
     插入则无 TSD 信号——TSD 缺失≠非转座,仅阳性有确证意义。
4. **语境注释**:GTF gene/exon → exon/intron/intergenic;
   rmsk BED → 是否落入已有 TE 区域(TE-into-TE vs 独特位点)。
5. **全长 L1 候选**:SVLEN≥6000 且家族=L1(潜在有活性元件清单)。

## 三、输出

| 文件 | 内容 |
| --- | --- |
| `mei_events.tsv` | per-event 全表:坐标/长度/组别/家族/亚家族/覆盖度/polyA/TSD/语境/基因/TE区域/全长L1标记 |
| `mei_family_by_group.tsv` | 家族 × 组别计数矩阵 |
| `mei_context.tsv` | 家族 × 组别 × 语境计数 |
| `mei_full_length_L1.tsv` | 全长 L1 候选清单 |
| `mei_summary.json` | 汇总指标 |
| `mei_fig1..6.png` | 图组(家族×组别 / 长度分布 / 语境 / TE-into-TE / 全长L1 / 指纹比例) |
| `mei_report.html` | 自包含报告(嵌图+表+方法+统计警告) |
| `repeatmasker/ins.fa.out` | RepeatMasker 原始输出(可复查亚家族) |

## 四、复现命令

```bash
# 分析(约 25 分钟,RepeatMasker 为主体;12 线程)
apptainer exec --bind /data/16T/Database --bind <output_dir> --bind <workflow_dir> \
  /data/16T/Database/env/sv/sv.sif python \
  <workflow_dir>/modules/sv/bin/utils/mei_analysis.py \
  -i <comparison>_merged.vcf \
  -o <comparison>_mei \
  --te-lib /data/16T/Database/Reference/famdb/rhesus_te_lib.fa \
  --gtf /data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/GCF_049350105.2_T2T-MMU8v2.0_genomic.sorted.gtf.gz \
  --te-bed /data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/rheMac10_rmsk_TE.lifted.bed \
  --fasta /data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/GCF_049350105.2_T2T-MMU8v2.0_genomic.fa \
  --threads 12

# 可视化 + HTML 报告(秒级)
apptainer exec ... /data/16T/Database/env/sv/sv.sif python \
  <workflow_dir>/modules/sv/bin/utils/mei_plot.py \
  --mei-dir <comparison>_mei --title "<label>"
```

调试:`--max-records 500` 先小样本试跑。

## 五、本数据集关键结果(luanchao_mulatta, youth_vs_aged, 1v1)

- 76,322 INS 全量分类:Alu 45,545 / L1 5,589 / LTR 1,530 / SVA 20
- 家族分布在两组间几乎对称(Alu 16,424 vs 16,829;L1 2,011 vs 1,988)
  ——**个体多态为主,无年龄方向性偏移**
- 全长 L1 候选 1,486 条(ctrl 610 / exp 611 / shared 265)
- 外显子插入 49 条(功能破坏候选,events 表按 context=exon 过滤)
- 语境:intergenic 54,288 / intron 21,985 / exon 49

## 六、已知限制

1. **1v1 设计**:每组 1 样本,"组特异"插入主体是个体多态而非组效
   应;组间定量比较需每组 ≥3 样本(脚本可直接复用,SUPP_VEC 自适应)。
2. **polyA 检出偏低**(26%):SURVIVOR/pbsv 对 ALT 序列有截断,
   末尾常缺失;全长 L1 子集内的检出率更有参考价值。
3. **TSD 仅阳性可确证**:参考位点已有祖先元件时断点两侧无重复序列。
4. RepeatMasker 分类不是比对到基因组拷贝,对嵌合/重排插入只报
   最长命中;`te_cov<0.6` 的记录分类可信度低。
5. 全长 L1 候选未做 ORF 完整性检查(内部缺失的假全长会混入);
   需要时可在子集上跑 ORFfinder 交叉验证。

## 七、后续扩展方向

- 多样本:merge 时每组多放几个样本,SUPP_VEC 变长,全部统计自动适应
- Fiber-seq 独有角度:插入位点的 m6a/FIRE 表观状态(活化 vs 沉默 TE)
- 外显子插入清单 → sv_go 富集(modules/function/sv_go,macaque 支持)
