# censatbuild（着丝粒/卫星区域集构建：censat_bed）

为组装构建着丝粒区域分类 BED（`genome.censat_bed`，供 modules/fiberseq/censat
使用）：Dfam 物种库卫星注释（复用 modules/centromere 的 repeatmasker 规则）
→ 阵列聚类/分类 → 加侧翼。另有可选的 CpG 甲基化 pileup 规则
（`genome.mcpg_pileup`）。

小鼠着丝粒结构：minor satellite（~120bp 单体，动粒核心）+ major satellite
（~234bp 单体，周着丝粒大阵列），与人类 alphaSat HOR 不同，按 repeat 名
分类而非 HOR。

## 规则

| 规则 | 输出（`<outdir>/`） |
|---|---|
| censatbuild_run | 聚合目标（censat bed + stats） |
| censatbuild_link_genome | `{genome}/assembly/asm.bp.p_ctg.fa`（复用布局的软链） |
| censatbuild_repeatmasker_run | `{genome}/RepeatMasker/asm.bp.p_ctg.fa.out/.tbl`（复用 centromere） |
| censatbuild_arrays | `<genome>.sat_arrays.bed` / `.sat_hits.bed`（聚类+分类） |
| censatbuild_flank | `<genome>.censat.bed` / `.censat.stats.tsv`（+侧翼） |
| censatbuild_mcg_pileup | `<genome>.mcpg_pileup.bedGraph`（可选，需 Params.censatbuild.bam） |

RepeatMasker 步骤直接 `use rule repeatmasker_run from centromere`（Dfam
物种库 + 既有 env），不重复造轮子；`censatbuild_link_genome` 把目标基因组
放到该规则约定的 hifiasm 布局路径上。

## FamDB（Dfam 数据库 + 环境版本配套，必须先准备好）

RepeatMasker 的 `-species` 查询依赖 FamDB（Dfam 的 h5 库）。两个配套问题：

1. **版本配套**：FamDB 文件有格式版本（Dfam 3.9 = format 2.0；Dfam 4.0 =
   format 3.0）。RepeatMasker 4.2.4 配 famdb.py 3.0.0（读 format 3.0）；
   RepeatMasker 4.2.3 自带旧 famdb.py（读 format 2.0）。**换环境比换库快**——
   本模块约定 `centromere` env 钉 `repeatmasker=4.2.3`，配 Dfam 3.9。
2. **库文件**：按物种分类分片（format 2.0），至少 root 片 + 目标物种所在片：

```
# 下载源 https://www.dfam.org/releases/Dfam_3.9/families/FamDB/
dfam39_full.0.h5     # Partition 0: root（必需）
dfam39_full.7.h5     # Partition 7: Mammalia（小鼠/人/猴），解压后约 56GB
```

放置于 `/data/16T/Database/env/famdb/`。RepeatMasker 从
`share/RepeatMasker/Libraries/famdb` 找库（容器内 sif 只读），运行时绑定：

```
--singularity-args "--bind /data/16T,/data/16T/Database/env/famdb:\
/opt/conda/envs/centromere/share/RepeatMasker/Libraries/famdb"
```

（apptainer 会自动创建容器内挂载点；PacVar 跑 RepeatMasker 时同样加这个 bind。）

## 方法与命令（模板）

```
# 卫星注释（centromere 模块的 repeatmasker_run 展开）
RepeatMasker -species 'Mus musculus' -pa <threads> \
    -dir <outdir>/<genome>/RepeatMasker -gff <outdir>/<genome>/assembly/asm.bp.p_ctg.fa

# 聚类+分类（只留 Satellite 类；minor→minorSat_core, major→majorSat_array）
python bin/censatbuild_rm2bed.py --rm-out <...>/asm.bp.p_ctg.fa.out \
    --arrays-out <genome>.sat_arrays.bed --hits-out <genome>.sat_hits.bed \
    --gap <gap> --major-pattern major --minor-pattern minor \
    --class-pattern Satellite

# 侧翼 + 汇总（flank 类别名 flank_<N>kb）
bedtools flank -i <genome>.sat_arrays.bed -g <target.fa.fai> -b <flank> \
    | awk 'BEGIN{OFS="\t"}{$4="flank_<N>kb"; print}' > flank.bed
cat <genome>.sat_arrays.bed flank.bed | sort -k1,1 -k2,2n > <genome>.censat.bed

# 可选：CpG 甲基化 pileup（bedGraph，供 censat 的 mcpg_pileup）
modkit pileup --cpg --combine-strands <bam> -o tmp.bedMethyl
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3,$11/100}' tmp.bedMethyl > mcpg_pileup.bedGraph
```

## config

```
genome.fasta        目标组装 fasta
genome.fai          目标组装 fai（bedtools flank 需要）
Params.censatbuild.species（RepeatMasker -species，默认 "Mus musculus"）/
        gap（阵列内最大间隙，默认 10kb）/ flank（侧翼，默认 10kb）/
        major_pattern / minor_pattern / class_pattern（默认 Satellite）/
        bam（可选，mcg_pileup 的输入 BAM）
```

## 环境

- RepeatMasker 步骤用 `centromere` env（`/data/16T/Database/env/centromere/`）。
- 聚类/侧翼/pileup 用 `censatbuild` env（repeatmasker + bedtools + samtools +
  ont-modkit + python），SIF 同 Database/env 约定构建。

## 运行（独立）

```
snakemake -s modules/censatbuild/censatbuild.smk --configfile censatbuild.json \
    --cores 32 --sdm apptainer \
    --singularity-args "--bind /data/16T,/data/16T/Database/env/famdb/famdb.conf:/opt/conda/envs/centromere/share/famdb-3.0.0/famdb.conf" \
    --rerun-triggers=mtime censatbuild_run
```

## CDR_core（数据依赖，后续补充）

CDR（着丝粒可及性凹陷区）是阵列内部的可及性凹陷，需用 Fiber-seq 数据定义
（Nat Genet 2025 的做法）。初版可用 `minorSat_core` 当 proxy；等主流程产出
`msp.bed.gz` 后在 majorSat 阵列内部滑窗找 MSP 覆盖凹陷段，把区间以
`CDR_core` 类别并入 censat.bed 即可（censat_profile 按 name 列分组，
类别可随时增补）。

## 目录约定

outdir 指向工作目录，Database 只放参考文件与数据库资源（famdb h5 属于
数据库资源放 `Database/env/famdb/`；最终 `censat.bed` 拷入 Reference 目录后
在 workflow config 登记 `genome.censat_bed`）。
