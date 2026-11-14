# disambiguate — 异种/多基因组比对归属消歧模块(ngs_disambiguate)

对同一 read 集在两个基因组(物种)上的比对结果做归属消歧:每读段唯一归属
Species A / Species B,或判 ambiguous,产出按基因组命名的拆分 BAM 与
QC 汇总。典型用途:人-鼠异种移植(xenograft)、宿主-微生物混合样本。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + `genome_pairs` | 形态 A 拼接,多基因组布局:`{indir}/{genomeA}/{sample_id}/{sample_id}.bam` 与 `{indir}/{genomeB}/{sample_id}/{sample_id}.bam`(两个基因组的比对 BAM) |
| `genome_pairs` | 长度恰为 2 的基因组名列表 `[genomeA, genomeB]`(parse 期解包);A/B 顺序决定输出命名 |
| `paired_samples` / `single_samples` | 样本列表(QC 汇总聚合用) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.disambiguatedSpecies_{genomeA}.bam   # 唯一归属 A(+ .bai)
{outdir}/{sample_id}/{sample_id}.disambiguatedSpecies_{genomeB}.bam   # 唯一归属 B(+ .bai)
{outdir}/{sample_id}/{sample_id}.ambiguousSpecies_{genomeA}.bam       # 歧义读段(A 侧)
{outdir}/{sample_id}/{sample_id}.ambiguousSpecies_{genomeB}.bam       # 歧义读段(B 侧)
{outdir}/{sample_id}/{sample_id}_summary_renamed.tsv                  # 每样本归属计数(表头换成基因组名)
{outdir}/disambiguate_qc.tsv                                          # 全样本 QC 汇总
```

中间产物(`.disambiguatedSpeciesA/B.bam`、`.ambiguousSpeciesA/B.bam`、
`*_summary.tsv`、name-sort 临时 BAM)均为 temp,消费后清理。
日志:`{logdir}/{sample_id}/ngs_disambiguate.log`、`sort_rename.log`、
`{logdir}/disambiguate_report.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志目录 |
| `genome_pairs` | — | `[]` | `[genomeA, genomeB]`,必需 |
| `paired_samples` / `single_samples` | — | `[]` / `[]` | 样本列表 |
| `Procedure.ngs_disambiguate` | `ngs_disambiguate -s ... -o ... -a ...` | `ngs_disambiguate` | 消歧引擎 |
| `Procedure.samtools` | `samtools sort -n/-@ ...`、`samtools index` | `samtools` | 排序/索引 |
| `Params.ngs_disambiguate.aligner` | `-a` | `hisat2` | 上游比对器(可选 Tophat/Hisat2/STAR/BWA mem,见 json 注释) |
| `ROOT_DIR` | — | `.` | 拼汇总脚本路径(见下注) |
| 线程/内存 | `-@` | `threads=4`,`mem_mb=16000` | `rule_res(...)` 可覆盖 |

注:`disambiguate_report` 规则的组合脚本路径为
`ROOT_DIR/modules/disambiguate/combineDisambiguateQC.py`,而脚本实际位于
`bin/combineDisambiguateQC.py` —— 接线/拷贝模块时注意该路径差异。
环境:`conda: "disambiguate.yaml"`(模块目录内);模块另有 `tests.sh`。

直接调用示例(等价引擎命令):

```bash
samtools sort -n -@ 4 -o s1.bamA.sortN.bam in/GRCm39/s1/s1.bam
samtools sort -n -@ 4 -o s1.bamB.sortN.bam in/GRCh38/s1/s1.bam
ngs_disambiguate -s s1 -o out/s1 -a hisat2 s1.bamA.sortN.bam s1.bamB.sortN.bam
python bin/combineDisambiguateQC.py -i s1_summary_renamed.tsv s2_summary_renamed.tsv \
  -o out/disambiguate_qc.tsv
```

## 2. 算法原理与生物学意义

ngs_disambiguate(AMDIS 套件)按比对得分比较同一读段在两个基因组上的
最优比对:显著更优者归属该物种,两者相当则判 ambiguous。异种移植样本中
人源(肿瘤)与鼠源(基质)读段必须分开定量,否则表达/突变分析会被交叉
污染;拆分后 `.disambiguatedSpecies_<genome>.bam` 即"纯净"物种 BAM,可直接
进各自基因组的下游分析。结果怎么读:`*_summary_renamed.tsv` 行为归属计数
(unique A/B、ambiguous),`disambiguate_qc.tsv` 汇总全样本;ambiguous 比例
过高提示两基因组同源区混淆或比对参数过松。

## 3. 参考文献

1. ngs_disambiguate — https://github.com/AstraZeneca-NGS/NGS_Disambiguate
2. Ahdesmäki MJ et al. *Disambiguate: An open-source application for
   augmented counting of xenografted tumour samples in RNA-Seq data*.
   GigaScience 5(suppl 1):s13058-016-0754-9, 2016(如引用请以原文 DOI 为准)
