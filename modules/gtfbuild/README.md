# gtfbuild（基因组 GTF 构建：注释投射）

为没有官方注释的组装构建基因注释 GTF。当前实现为同物种注释投射
（LiftOff + minimap2 + gffread），与 Ensembl Rapid Release / RefSeq 新组装
的投射做法同源；后续可扩展 BRAKER3（RNA-seq + 蛋白证据）等路线。

## 规则

| 规则 | 输出（`<outdir>/<genome_name>.`） |
|---|---|
| gtfbuild_run | 聚合目标（gtf + qc） |
| gtfbuild_liftoff | .lifted.gff3 / .unmapped.gff3 / .gtf |
| gtfbuild_qc | .gtf_qc.tsv（映射率/完整模型/多余拷贝） |
| gtfbuild_bed_run | 聚合目标（lifted bed + stats） |
| gtfbuild_wga | .wga.paf（目标 vs 来源全基因组比对） |
| gtfbuild_liftover_bed | `<bed名>.lifted.bed` / `<bed名>.lifted.stats.tsv` |

区间文件（如 RepeatMasker BED）走 `gtfbuild_wga` + `gtfbuild_liftover_bed`
坐标投射（paftools.js liftover，经共线性长比对块携带），不能按序列直接映射——
TE 拷贝间近乎相同，会落到旁系同源位点上。

## config

`genome.fasta`（目标组装）、`genome.source_fasta`（来源组装）、
`genome.source_gtf`（来源注释，如 GENCODE）；`genome_name` 用于输出命名，
请在文件名中标注来源（如 `GCA_xxx_vM38_liftoff`）避免与官方注释混淆。
`Params.gtfbuild.*`：threads / copies（搜旁系同源拷贝）/ coverage（最低覆盖度）。

## 环境

`gtfbuild.yaml`（liftoff + minimap2 + gffread + samtools），SIF: `gtfbuild.sif`。

## 独立运行

```
snakemake -s modules/gtfbuild/gtfbuild.smk --configfile gtfbuild.json \
    --cores 48 --sdm apptainer
```
