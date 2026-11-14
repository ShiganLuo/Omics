# fiberseq/track — 基因组浏览器轨道(bigWig)生成子模块

把 Fiber-seq 元素 BED(MSP / m6A / 核小体)汇总为分箱 bedGraph 并转
bigWig,附 UCSC trackDb 片段,供 IGV/UCSC 浏览器展示。引擎脚本为父模块
共享的 `modules/fiberseq/bin/fiberseq_track.py`,conda 环境共享父目录
`../fiberseq.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 上游元素目录,默认 `common/4_elements`;input 函数拼 `{indir}/{sample_id}/{sample_id}.msp.bed.gz`、`.m6a.bed.gz`、`.nuc.bed.gz`(BED12,需按染色体分组排序 `sort -k1,1`) |
| `genome.fasta` + `.fai` | 染色体长度取 `genome.fasta` 路径 + `.fai`(`-c` 接受 .fai 前两列;JSON 模板里为 `genome.fai` 键,规则实际拼 `fasta + ".fai"`) |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.accessibility.bw   # MSP 可及性轨道 bigWig
{outdir}/{sample_id}/{sample_id}.m6a_density.bw     # m6A 密度轨道 bigWig
{outdir}/{sample_id}/trackDb.txt                    # UCSC trackDb 片段
```

(转换成功后中间 bedGraph 即删;`bedGraphToBigWig` 不在 PATH 时保留
bedGraph 并告警。)收尾规则 `fiberseq_track_result` 聚合 accessibility
bigWig。

### 1.3 参数

config 位于 `Params.fiberseq_track` 下,与 `fiberseq_track.py` CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `bin` | `-b/--bin` | 10 | 分箱宽度(bp) |
| — | `-s/--msp`, `-a/--m6a`, `-u/--nuc` | — | 三类元素 BED12(.bed.gz;nuc 规则内必传) |
| — | `-c/--chrom-sizes` | — | 染色体长度(.fai 或 .chrom.sizes) |
| — | `-o` / `--sample` | — | 输出目录 / 样本名前缀 |
| — | `-l` | — | 日志文件 |
| — | `-n/--dry-run` | off | 只打印计划不产出 |
| threads / mem_mb | — | 4 / 65000 | `rule_res("fiberseq_track")` 配置 |

直接调用示例(脱离流程):

```bash
python modules/fiberseq/bin/fiberseq_track.py \
  -s S1.msp.bed.gz -a S1.m6a.bed.gz -u S1.nuc.bed.gz \
  -c genome.fa.fai -o out/S1 --sample S1 -b 10 -l track.log
```

## 2. 算法原理与生物学意义

按染色体流式读入三类 BED(同一遍推进,峰值内存 = 最大染色体的记录),
以 `bin`(默认 10 bp)分箱聚合两个信号:

- **accessibility**:MSP 纤维占比(`block_signal(uniq_fibers=True)`,
  每 bin 内被 MSP 覆盖的唯一流数 / 覆盖流数),值 > 0 的 bin 写入
  bedGraph;
- **m6A density**:m6A 位点密度(按纤维计数,`uniq_fibers=False`)除以
  该 bin 的覆盖流数(msp+nuc),即每纤维每 bin 的 m6A 密度。

bedGraph 经 `bedGraphToBigWig` 转 bigWig(成功删中间 bedGraph),
trackDb.txt 给出 `track type=bigWig` 条目。生物学意义:可及性轨道把
Fiber-seq 的单分子开放染色质证据聚合成连续信号,叠加 m6A 密度可直观
看到调控元件(高可及 + 高 m6A)、核小体占位与样本间差异,是 IGV 检查
上游元素调用质量与展示结果的标准轨道。

## 3. 参考文献

1. FIRE analysis pipeline — https://github.com/fiberseq/FIRE(元素 BED
   约定)。
2. Kent WJ et al. *The human genome browser at UCSC*. Genome Research
   12(6):996-1006, 2002. DOI: 10.1101/gr.229102(bigWig/trackDb 所属
   浏览器体系);bedGraphToBigWig — http://hgdownload.soe.ucsc.edu/admin/exe/
3. 引擎脚本 `bin/fiberseq_track.py` 为本仓库自研(共享 `fiberseq/bin/`)。
