# gtfbuild（基因组注释构建：GTF 投射 + 区间坐标投射）

为没有官方注释的组装构建基因注释 GTF，并把区间注释（如 RepeatMasker TE BED）
投射到目标组装坐标。与 Ensembl Rapid Release / RefSeq 新组装的投射做法同源；
后续可扩展 BRAKER3（RNA-seq + 蛋白证据）等路线。

## 规则

| 规则 | 输出（`<outdir>/`） |
|---|---|
| gtfbuild_run | 聚合目标（gtf + qc） |
| gtfbuild_liftoff | `<genome>.lifted.gff3` / `.unmapped.gff3` / `.gtf` |
| gtfbuild_qc | `<genome>.gtf_qc.tsv`（映射率/完整模型/多余拷贝） |
| gtfbuild_bed_run | 聚合目标（lifted bed + stats） |
| gtfbuild_wga | `<genome>.wga.paf`（目标 vs 来源全基因组比对，可复用） |
| gtfbuild_liftover_bed | `<bed名>.lifted.paftools.bed` / `.lifted.bed` / `.lifted.stats.tsv` |

## 方法原理

两条路线按注释类型选择：

- **基因 GTF → LiftOff 基因投射**（序列+基因结构）。同物种（如同为 C57BL/6）
  外显子边界基本一致，投射质量高；UTR/可变剪接继承来源注释的转录本结构。
- **区间 BED（TE 等）→ WGA 共线坐标投射**。TE 拷贝间序列近乎相同，按序列映射
  会落到旁系同源位点；必须先做全基因组比对，让区间由 >=50kb 的共线性长比对块
  携带到直系同源位置（paftools.js liftover 的 -l/-q 过滤）。

## 一键运行

```
snakemake -s /data/16T/luosg/workflow/Omics/modules/gtfbuild/gtfbuild.smk \
    --configfile <outdir>/gtfbuild.json --cores 48 \
    --sdm apptainer --singularity-args "--bind /data/16T" \
    --rerun-triggers=mtime gtfbuild_run gtfbuild_bed_run
```

目标用聚合规则（`gtfbuild_run` / `gtfbuild_bed_run`）或具体输出文件
（注意用绝对路径匹配规则输出）。config 示例见 `gtfbuild.json`，字段：

```
genome.fasta        目标组装 fasta（同时投 GTF 和 BED 时共用）
genome.source_fasta 来源组装 fasta
genome.source_gtf   来源注释（如 GENCODE）
genome.source_bed   来源区间 BED（如 RepeatMasker TE）
genome_name         输出命名前缀（建议标注来源，如 GCA_xxx_vM38_liftoff）
Params.gtfbuild.threads / copies / coverage / wga_preset /
        liftover_min_mapq / liftover_min_aln / extra_args
```

## 规则执行的命令（模板）

GTF（gtfbuild_liftoff + gtfbuild_qc）：

```
liftoff -g <source.gtf> -o <out>.lifted.gff3 -u <out>.unmapped.gff3 \
    -p <threads> -a <coverage> -dir <out>.liftoff_intermediate \
    <target.fa> <source.fa>
gffread <out>.lifted.gff3 -T -o <out>.gtf
python bin/gtfbuild_add_gene_features.py --gtf <out>.gtf --out tmp \
    --genome <genome_name> --source-gtf <source.gtf> && mv tmp <out>.gtf
python bin/gtfbuild_qc.py --src-gtf <source.gtf> --lifted-gff3 <out>.lifted.gff3 \
    --lifted-gtf <out>.gtf --unmapped <out>.unmapped.gff3 --out <out>.gtf_qc.tsv
```

- `-a` 为映射覆盖度阈值（默认 0.5）；gffread 转 GTF 后不产 gene 行，
  `gtfbuild_add_gene_features.py` 按转录本 span 补 gene 特征
  （下游 fire_annotate / nucpos 按 gene+exon 分类和取 TSS）。

BED（gtfbuild_wga + gtfbuild_liftover_bed）：

```
minimap2 -cx asm5 -t <threads> --secondary=no <target.fa> <source.fa> > <out>.wga.paf
paftools.js liftover -q 5 -l 50000 <out>.wga.paf <source.bed> > <out>.lifted.paftools.bed
python bin/gtfbuild_restore_bed.py --lifted <out>.lifted.paftools.bed \
    --source <source.bed> --out <out>.lifted.bed --genome <genome_name>
```

- paftools.js 输出丢原注释列（只留 old-locus name）；`gtfbuild_restore_bed.py`
  按 old-locus（name 列 `chrom_start_end`）join 回源行，输出 BED8：
  新坐标 / name(TE名) / . / strand / family / class（name 列为下游约定的
  "重复类型"列，family/class 保留在 7/8 列）。
- `.wga.paf` 可复用：再投其他区间文件（censat、其他 repeat 库）无需重跑比对。

## 实例：T2T_mhaESC_v1.5 ← GRCm39 的完整命令

输入文件：

| 用途 | 路径 |
|---|---|
| 目标基因组（T2T_mhaESC_v1.5） | `/data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_genomic.fna` |
| 来源基因组（GRCm39） | `/data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa` |
| 来源注释 GTF（GENCODE vM38 basic） | `/data/16T/Database/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf` |
| 来源 TE BED（RepeatMasker） | `/data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.bed` |

输出文件（安装进 Reference 的最终产物）：

| 产物 | 路径 |
|---|---|
| 投射 GTF | `/data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.gtf` |
| 投射 TE BED | `/data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GRCm39_GENCODE_rmsk_TE.lifted.bed` |

实际使用的 config（`/data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/gtfbuild.json`）：

```json
{
    "ROOT_DIR": "/data/16T/luosg/workflow/Omics",
    "indir": "",
    "outdir": "/data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro",
    "logdir": "/data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/log",
    "samples": [],
    "genome_name": "GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff",
    "Procedure": {},
    "Params": {
        "gtfbuild": {
            "threads": 48,
            "copies": false,
            "coverage": 0.5,
            "extra_args": ""
        }
    },
    "genome": {
        "fasta": "/data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_genomic.fna",
        "source_fasta": "/data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa",
        "source_gtf": "/data/16T/Database/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf",
        "source_bed": "/data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.bed"
    },
    "env": {
        "env_dir": "/data/16T/Database/env",
        "gtfbuild": "/data/16T/Database/env/gtfbuild/gtfbuild.sif"
    }
}
```

运行命令：

```
snakemake -s /data/16T/luosg/workflow/Omics/modules/gtfbuild/gtfbuild.smk \
    --configfile /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/gtfbuild.json \
    --cores 48 --sdm apptainer --singularity-args "--bind /data/16T" \
    --rerun-triggers=mtime gtfbuild_run gtfbuild_bed_run
```

展开后各规则实际执行（容器 gtfbuild.sif 内，-p 48）：

```
# gtfbuild_liftoff
liftoff -g /data/16T/Database/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
    -o /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.lifted.gff3 \
    -u /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.unmapped.gff3 \
    -p 48 -a 0.5 \
    -dir /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.liftoff_intermediate \
    /data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_genomic.fna \
    /data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa

gffread /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.lifted.gff3 \
    -T -o /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.gtf

python /data/16T/luosg/workflow/Omics/modules/gtfbuild/bin/gtfbuild_add_gene_features.py \
    --gtf  /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.gtf \
    --out  /data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.gtf.tmp
mv ...gtf.tmp ...gtf

# gtfbuild_qc
python /data/16T/luosg/workflow/Omics/modules/gtfbuild/bin/gtfbuild_qc.py \
    --src-gtf /data/16T/Database/Reference/mouse/GENCODE/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf \
    --lifted-gff3 ...lifted.gff3 --lifted-gtf ...gtf --unmapped ...unmapped.gff3 \
    --genome GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff --out ...gtf_qc.tsv

# gtfbuild_wga
minimap2 -cx asm5 -t 48 --secondary=no \
    /data/16T/Database/Reference/mouse/NCBI/T2T_mhaESC_v1.5/GCA_056825265.1_T2T_mhaESC_v1.5_genomic.fna \
    /data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa \
    > ...GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.wga.paf

# gtfbuild_liftover_bed
paftools.js liftover -q 5 -l 50000 ...wga.paf \
    /data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.bed \
    > ...GRCm39_GENCODE_rmsk_TE.lifted.paftools.bed

python /data/16T/luosg/workflow/Omics/modules/gtfbuild/bin/gtfbuild_restore_bed.py \
    --lifted ...lifted.paftools.bed \
    --source /data/16T/Database/Reference/mouse/GENCODE/GRCm39/GRCm39_GENCODE_rmsk_TE.bed \
    --out ...GRCm39_GENCODE_rmsk_TE.lifted.bed
```

（`...` 为 `/data/16T/luosg/output/gtfbuild/T2T_mhaESC_v1.5_repro/` 前缀；
每次运行的完整脚本留档在 outdir 的 `gtfbuild_*.sh`。）

结果（2026-10-10，安装进 Reference 的产物）：

- GTF：gencode.vM38.basic → `GCA_056825265.1_T2T_mhaESC_v1.5_vM38_liftoff.gtf`
  （77,666 gene / 127,112 transcript / 671,648 exon；gene 映射率 99.15%）
- BED：GRCm39_GENCODE_rmsk_TE.bed（3,799,630 条）→ `GRCm39_GENCODE_rmsk_TE.lifted.bed`
  （lift_rate 99.32%；SINE 148 万 / LTR 110 万 / LINE 97 万 / DNA 16 万）
- 复现验证（module 全流程重跑 vs 安装产物）：BED 逐字节一致；GTF 99.996% 一致，
  唯一差异是 1 个旁系同源假基因位点（Cwc22rt2/Cwc22rt4）上 LiftOff 平局裁决
  的不确定性，QC 指标完全一致。
- fire_annotate 正例验证：peak 落在投射 TE 内 → repeat_flag=yes、repeat_names
  输出 TE 名；peak 在基因 TSS 附近 → promoter + gene_name。

## 环境

`gtfbuild.yaml`（python + liftoff + minimap2 + gffread + samtools）。
SIF 构建（与 Database/env 其他 env 同约定）：

```
cp gtfbuild.yaml /data/16T/Database/env/gtfbuild/
cd /data/16T/Database/env/gtfbuild && apptainer build gtfbuild.sif gtfbuild.def
```

## 目录约定

**outdir 指向工作目录**，Database 只放最终参考文件：运行完成后把 `<gtf>`
与 `<bed>` 拷入对应的 `Database/Reference/...` 目录，并在 workflow config
（如 `config/Fiberseq.json` 的 `genome.references.<组装>`）登记路径。
日志/脚本/PAF/中间 gff3/liftoff intermediate 都留在工作目录。

## QC 指标含义

- `.gtf_qc.tsv`：gene_mapping_rate（投射基因占比）、partial_transcripts、
  extra_copy_transcripts、unmapped_features。
- `.lifted.stats.tsv`：lift_rate = out_intervals / in_intervals；
  被 -q/-l 过滤掉的多为结构变异区域。

## 输出头注释

两个最终产物顶部都有 `#` 注释头（GTF/BED 解析器均跳过 `#` 行，下游
fire_annotate/nucpos/read_bed 不受影响）：

- GTF：`#!genome-build`（目标组装）、`#!annotation-source`（来源注释）、
  `#!annotation-method`（投射方法）、`#!note`（提示为投射注释）。
- BED：`#target-genome`、`#lifted-from`（来源 BED + 投射方法）、`#note`、
  列名行 `#chrom start end name score strand family_id class_id`。

## 局限

- 投射注释 ≠ 官方注释：无 RNA-seq 证据、无人工校对；输出文件名必须标注来源
  （如 `_vM38_liftoff`）避免与官方注释混淆。官方 GTF 可用后替换 config 重跑即可。
- 结构变异区域的基因/区间会被阈值过滤（本实例 TE 约 0.7% 未投上）。
- 同物种投射假设共线性；跨物种使用请调低 wga_preset（asm10/asm20）并复核 QC。
- LiftOff 在旁系同源重复位点（如假基因多拷贝）的取舍顺序不确定，
  重跑可能在个别位点换一个等价拷贝（坐标/QC 不变）；需要完全确定可开
  `Params.gtfbuild.copies=true` 保留全部拷贝。
