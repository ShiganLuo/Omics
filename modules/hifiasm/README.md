# hifiasm — HiFi 基因组组装模块

对 PacBio HiFi BAM 从头组装,产出主 contig FASTA。独立成模块供多个下游
消费(centromere 着丝粒分析、telomere `assembly_scan`/`tidk` 等)。
输出布局契约:`{outdir}/{sample_id}/assembly/asm.bp.p_ctg.fa`,下游按
`assembly_dir` 配置键或共享 outdir 常量接线。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `inputs.hifiasm.bam` / `inputs.hifiasm.pbi` | dict-mode(优先):样本 ID → 完整 BAM / `.bam.pbi` 路径字典(Fiberseq 的 `{raw}/{sample}/{sample}.bam` staging 布局) |
| `indir` + `bam_substring` | 形态 A(缺省):`{indir}/{sample_id}[.{bam_substring}].bam` + 同名 `.bam.pbi` 索引(PacBio HiFi BAM) |
| `samples` | 样本列表 |

### 1.2 输出

```
{outdir}/{sample_id}/assembly/asm.bp.p_ctg.fa        # 组装主 contig(下游契约路径)
{outdir}/{sample_id}/assembly/asm.bp.p_ctg.gfa       # hifiasm 主 contig 图(工具产物)
{outdir}/{sample_id}/assembly/asm.hifiasm.stamp      # GFA 复用指纹(BAM 路径+大小+mtime)
```

日志:`{logdir}/{sample_id}/hifiasm_assemble.log`。工具自身输出(hifiasm
`[M::...]` 阶段统计、samtools 提示)直接追加在**同一规则日志**里;注意
hifiasm 在读段加载/纠错期几乎不输出,统计集中在运行结束时成批打印。引擎
另在日志里补时间线:工具版本、Stage 1/2(流式组装)与 Stage 2/2(GFA→FASTA)
的边界与耗时、contig 数/总长/N50,便于中途判断进度。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `logs` | 目录 |
| `bam_substring` | 文件名插入段 | `""` | 见输入(形态 A) |
| `inputs.hifiasm.bam` / `.pbi` | — | 无 | 见输入(dict-mode,优先于形态 A) |
| `Procedure.hifiasm` / `samtools` / `python` | — | 同名命令 | 可执行路径 |
| `--force-reassemble` | 引擎旗标 | 关 | 忽略已有 GFA 与 hifiasm 断点缓存,强制重新组装;缺省时 GFA 带指纹 stamp 可跨重跑复用(阶段 2 失败/中断后重跑不重装) |
| 线程/内存 | — | 32 / `mem_mb=64000` | 真实入口是 `resources.hifiasm_assemble.threads/.mem_mb`(`rule_res()` 解析并按宿主预算钳制;hifiasm 线程扩展到 ~32 收益明显,再往上递减) |

环境:`conda: "hifiasm.yaml"`。**SIF 复用**:容器模式下显式映射
`"hifiasm": "<env_dir>/centromere/centromere.sif"`(探针验证该镜像完整覆盖
本模块 yaml 依赖;centromere.sif 由含 hifiasm/samtools 的 centromere.yaml
治理,升级或重建后须对本模块重跑探针:import 覆盖 + `--help` 冒烟)。

直接调用示例(模块引擎命令,与回执脚本一致):

```bash
# samtools fastq 直接管道进 hifiasm,无中间 FASTQ;GFA→FASTA 由引擎完成
python modules/hifiasm/bin/hifiasm_assemble.py \
  --bam s1.bam --out-prefix out/s1/assembly/asm \
  --out-fasta out/s1/assembly/asm.bp.p_ctg.fa --threads 32
```

## 2. 算法原理

hifiasm 用串图(BOG)组装 HiFi 读段,纠错后输出单倍型解析
(`bp.hap1/hap2.p_ctg`)与主 contig(`bp.p_ctg`)图;本模块取主 contig 的
GFA `S` 记录转为 FASTA(即 `asm.bp.p_ctg.fa`)。引擎把 `samtools fastq`
的输出流直接接 hifiasm 的 `/dev/stdin`(实测 0.25.0 支持管道输入),
不落中间 FASTQ:省一次全量写读、免单线程 gzip 压缩成本。

## 3. 参考文献

1. hifiasm — https://github.com/chhylp123/hifiasm(Cheng H et al.,
   Nature Methods 18:608–610, 2021. DOI: 10.1038/s41592-020-01056-5)
