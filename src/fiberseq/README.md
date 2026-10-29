# fiberseq src scripts（Fiber-seq 下游专题分析，配合 modules/fiberseq 模块）

`modules/fiberseq/` 的规则级模块（annotate/motif/footprint/nucpos/quant/
diff/coactuation/haplotype/censat/track）覆盖常规 fire 下游；本目录放
一次性/专题脚本：变异 circos、TE 转位事件、TE 染色质景观、eccDNA 环化
read 检测、端粒长度/组成、区域（中心粒/端粒/线粒体）剖析、fiber×scRNA
关联。全部跑在 `fiberseq.sif`（numpy/pandas/scipy/matplotlib/pysam）；
例外：`fiber_scrna_link.py` 需要 anndata，跑在 `scanpy.sif`。

数据（本次验证用全能干细胞，小鼠 T2T_mhaESC_v1.5，DMSO=pluripotency vs
PlaB=totipotency；猴子对照运行见文末）：

| 输入 | 路径 |
|---|---|
| 元件 BED（每样本） | `/data/16T/luosg/output/Totipotency/Fiberseq/results/{DMSO,PlaB}/{sample}.{msp,nuc,cpg,m6a}.bed.gz` |
| fire BAM | `.../results/{DMSO,PlaB}/{sample}.fire.sorted.bam` |
| SNV/indel VCF | `.../variation/germline_snv_indel/{sample}/{sample}.vcf.gz` |
| SV VCF | `.../variation/germline_sv/{sample}/{sample}.sv.vcf.gz` |
| TE 注释（GRCm39→T2T WGA 投影） | `/data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GRCm39_GENCODE_rmsk_TE.lifted.bed` |
| RepeatMasker .out（censatbuild 产物） | `/data/16T/luosg/output/censatbuild/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5/RepeatMasker/asm.bp.p_ctg.fa.out` |
| fai | `/data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_genomic.fa.fai` |

输出统一在 `/data/16T/luosg/output/Totipotency/Fiberseq/analysis/<topic>/`。

## 脚本与一条命令用法

环境前缀（下文省略）：
```
SIF="apptainer exec --bind /data/16T:/data/16T /data/16T/Database/env/fiberseq/fiberseq.sif python"
cd /data/16T/luosg/workflow/Omics
A=/data/16T/luosg/output/Totipotency/Fiberseq
```

### 1. build_region_bed.py — 中心粒/端粒/线粒体区域 BED

从 RepeatMasker .out + fai 构建三类区域（col4=class）：
tel_array_5p/3p + subtel_flank、majorSat_array/minorSat_core/otherSat +
flank_10kb、MT_genome。小鼠卫星名与 Dfam 不同：GSAT_MM=major，
MMSAT4/IMPB=minor，用 `--major-pattern/--minor-pattern` 匹配。

```
$SIF src/fiberseq/build_region_bed.py \
  --rm-out /data/16T/luosg/output/censatbuild/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5/RepeatMasker/asm.bp.p_ctg.fa.out \
  --fai /data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_genomic.fa.fai \
  -o $A/analysis/fiberseq_regions --mode all \
  --major-pattern GSAT --minor-pattern MMSAT
# -> asm.telomere.bed asm.censat.bed asm.mito.bed (merge: asm.all_regions.bed)
```
已验证输出：majorSat 118 / minorSat 409 / otherSat 17824 阵列 +
36702 flank；端粒 23 个末端阵列；mito=CM166269.1 (16,299 bp)。

### 2. region_profile.py — 区域类开放性/甲基化（中心粒/端粒/线粒体）

对每个区域类算 MSP(>=50bp) 开放性（pooled + per-fiber）、核小体占位、
m6A/mCpG 每 kb read overlap 密度，类级 rank-sum 检验 + 归一化 metaprofile。

```
$SIF src/fiberseq/region_profile.py \
  -r $A/analysis/fiberseq_regions/asm.all_regions.bed \
  -s DMSO=$A/results/DMSO/DMSO.msp.bed.gz,$A/results/DMSO/DMSO.nuc.bed.gz,$A/results/DMSO/DMSO.cpg.bed.gz,$A/results/DMSO/DMSO.m6a.bed.gz \
  -s PlaB=$A/results/PlaB/PlaB.msp.bed.gz,$A/results/PlaB/PlaB.nuc.bed.gz,$A/results/PlaB/PlaB.cpg.bed.gz,$A/results/PlaB/PlaB.m6a.bed.gz \
  -o $A/analysis/region_profile -l $A/analysis/region_profile/region_profile.log
# -> region_metrics.tsv region_class_summary.tsv metaprofile_*.tsv *.png/pdf
```
注意：流式扫 4 个元件文件（每样本 ~2-17 GB gzip），首次运行约 1-2 h。

### 3. te_landscape.py — TE 家族染色质景观（开放性 + 5mC + 6mA）

每 TE 拷贝算 msp_frac/nuc_frac/mcpg/m6a per kb，聚合 class/family/
subfamily 三级 + rank-sum，输出 heatmap。回答"重复序列上开放性变化"
"5mC/6mA 在 TE 上的分布"。

```
$SIF src/fiberseq/te_landscape.py \
  --te-bed /data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GRCm39_GENCODE_rmsk_TE.lifted.bed \
  -s DMSO=$A/results/DMSO/DMSO.msp.bed.gz,...（同上 4 件套） \
  -s PlaB=$A/results/PlaB/PlaB.msp.bed.gz,... \
  -o $A/analysis/te_landscape --level family
# -> te_copy_metrics.tsv te_{class,family,subfamily}_summary.tsv te_family_heatmap.png
```

### 4. var_circos.py — 变异/信号环状图（纯 matplotlib circos）

四环密度（SNV、FIRE peak、m6A、mCpG，1 Mb bin）+ SV 弧（按类型着色）。
每个样本一张。

```
$SIF src/fiberseq/var_circos.py --fai <fai> \
  --vcf $A/variation/germline_snv_indel/DMSO/DMSO.vcf.gz \
  --sv-vcf $A/variation/germline_sv/DMSO/DMSO.sv.vcf.gz \
  --fire-peaks $A/results/DMSO/DMSO.fire_peaks.bed \
  --m6a $A/results/DMSO/DMSO.m6a.bed.gz --cpg $A/results/DMSO/DMSO.cpg.bed.gz \
  --sample DMSO -o $A/analysis/circos/DMSO -f png pdf
# -> DMSO.circos.{png,pdf} circos_bin_density.tsv  （两样本均已跑完）
```

### 5. sv_te_events.py — TE 相关 SV / 转位事件

SV 断点 ±100bp 与 TE 注释 intersect；SVtype×TEclass/family 计数矩阵；
interchrom 事件表；TE-class 断点富集（置换检验，--n-perm 0 可跳过）。

```
$SIF src/fiberseq/sv_te_events.py \
  --sv-vcf $A/variation/germline_sv/DMSO/DMSO.sv.vcf.gz \
  --te-bed <TE lifted bed> --fai <fai> --sample DMSO \
  -o $A/analysis/sv_te --n-perm 100
# -> DMSO.sv_te_events.tsv DMSO.interchrom_events.tsv
#    DMSO.svtype_x_{class,family}.tsv DMSO.te_class_enrichment.tsv
```
已验证计数：DMSO 77,647 SV 中 45,693 TE 关联、382 interchrom；
PlaB 77,110 SV 中 45,709 TE 关联、304 interchrom。

### 6. eccdna_scan.py — eccDNA 环化 read 检测 + 组成 + 甲基化

在 fire BAM 找 concatemer read（周期性大 internal deletion，>=3 个
周期、周期一致 CV<25%），聚类成环位点，TE 组成（class:family 覆盖率），
并用 tabix 提取位点上的 m6A/mCpG calls（eccDNA 是否甲基化）。这说明
Fiber-seq 的长读确实能直接观测滚环扩增的 eccDNA（decoil 依赖 SV caller
的 BND，会漏掉小的/低丰度环）。

```
$SIF src/fiberseq/eccdna_scan.py \
  --bam $A/results/DMSO/DMSO.fire.sorted.bam --te-bed <TE lifted bed> \
  -s DMSO=$A/results/DMSO/DMSO.msp.bed.gz,...（4 件套） \
  --sample-name DMSO -o $A/analysis/eccdna
# -> DMSO.eccdna_loci.tsv DMSO.eccdna_size_hist.tsv （每样本 ~9 min）
```
已验证：DMSO 613 loci / PlaB 628 loci；多数环 200-300 bp（小卫星类），
部分 kb 级环 TE 覆盖 50-87%（L1/ERV 为主）。

### 7. telomere_profile.py — 端粒长度 + 重复变异体组成

fire BAM 直接扫读末端（TTAGGG|CCCTAA 滑窗，60% 阈值），长度分布 +
motif 组成（TTAGGG/TCAGGG/TTGGGG/... 比例）+ 双样本对比。

```
$SIF src/fiberseq/telomere_profile.py \
  -s DMSO=$A/results/DMSO/DMSO.fire.sorted.bam \
  -s PlaB=$A/results/PlaB/PlaB.fire.sorted.bam \
  -o $A/analysis/telomere -f png pdf --max-reads 2000000
# -> telomere_per_read.tsv telomere_summary.tsv telomere_composition.tsv violin
```

### 8. fiber_scrna_link.py — Fiber-seq × scRNA(scTE) 关联（scanpy.sif）

TE 家族层：te_landscape 的 copy 表 × scTE h5ad 的 per-cell-type TE 表达
（raw 层 gene_type=='TE'），Spearman；可选基因 promoter 层（--gtf --msp）。

```
apptainer exec --bind /data/16T:/data/16T /data/16T/Database/env/scanpy/scanpy.sif \
  python src/fiberseq/fiber_scrna_link.py \
  --h5ad /data/16T/luosg/output/luanchao_mulatta/scRNAseq/common/5_combine_h5ad/ovaries/ovaries_scTE_auto.h5ad \
  --te-metrics $A/analysis/te_landscape/te_copy_metrics.tsv \
  --te-level family -o /path/to/link_out
```

## 已知边界

- **5hmC 测不到**：fibertools fire 产物只有 m6A 和 mCpG 两类 mark；5hmC
  需要 PacBio kinetics 专用 caller（或 ONT+Remora），后续如有 modkit
  pileup 可接入 region_profile 的 `--pileup` 思路（同 censat_profile）。
- eccdna_scan 的 --min-period 200 起，<200 bp 的环不可见；周期性 internal
  deletion 检测在高错误率读上有假阳性，看结果时优先 period_cv 低、
  n_reads 高的位点。
- var_circos 的 SV 弧是折线近似（polar plot 限制），数量大时 alpha 降低。

## 猴子（luanchao_mulatta）对照运行

参考与数据（注意 T2T 参考无 MT contig，线粒体分析做不了）：

```
MA=/data/16T/luosg/output/luanchao_mulatta/Fiberseq
MF=/data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/GCF_049350105.2_T2T-MMU8v2.0_genomic.fa
# 元件: $MA/results/{luanchao-11236,luanchao-21310}/<sample>.{msp,nuc,cpg,m6a}.bed.gz
# SV/SNV: 无 —— 先跑 modules/deepvariant + pbsv（430-460 GB BAM，数小时/样本）
# TE 注释: 无 —— 先跑 modules/gtfbuild gtfbuild_bed_run:
#   source_bed=/data/16T/Database/Reference/mulatta/ENSEMBL/Mmul_10/rheMac10_rmsk_TE.bed
#   (chrN -> 需先重命名为 1..20,X,Y 匹配 Mmul_10 fa 的 Ensembl 命名，再 WGA 投影)
# 端粒: telomere_profile.py 直接可用（不需要任何注释）
# 区域: build_region_bed.py 需要 RepeatMasker .out —— FamDB 缺失，
#   替代: 用 gtfbuild 投影 mouse/human alphaSat，或装 FamDB 后跑 censatbuild
```

猴样本 11236 vs scRNA 的 11238 编号不一致，配对分析前需确认是否同个体；
11310/21310 可直接配 fiber_scrna_link。
