# fiberseq/annotate — FIRE peaks 基因组上下文注释子模块

对上游 Fiber-seq FIRE peak(`{sample_id}.fire_peaks.bed`)做基因组上下文
注释:把每个 peak 分类为 promoter / exon / intron / intergenic,并按 repeat
注释(如 RepeatMasker BED)标记重复序列重叠。引擎脚本为父模块共享的
`modules/fiberseq/bin/fire_annotate.py`,conda 环境共享父目录
`../fiberseq.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | 上游元素目录,默认 `common/4_elements`;input 函数拼 `{indir}/{sample_id}/{sample_id}.fire_peaks.bed`(FIRE peak BED) |
| `genome.gtf`(config) | 基因注释 GTF(必填,缺失 input 函数报错) |
| `genome.repeat_bed`(config) | 可选:repeat 注释 BED(如 RepeatMasker);缺失或不存在时跳过 repeat 重叠并告警 |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.peaks_annotated.tsv  # 逐 peak 注释表
{outdir}/{sample_id}/{sample_id}.peaks_summary.tsv    # 分类计数汇总
```

`peaks_annotated.tsv` 列:`chrom start end name score region_class
nearest_gene repeat_flag repeat_frac repeat_names`;
`peaks_summary.tsv` 列:`category key count`(region_class / repeat_flag /
total)。收尾规则 `fiberseq_annotate_result` 聚合全部样本两表。

### 1.3 参数

config 位于 `Params.fiberseq_annotate` 下,与 `fire_annotate.py` CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `tss_flank` | `-t/--tss-flank` | 1000 | TSS 两侧 promoter 窗口(bp) |
| `repeat_frac` | `-f/--repeat-frac` | 0.25 | peak 与 repeat 重叠碱基占比 ≥ 此值才标记 `repeat_flag=yes` |
| — | `-p/--peaks` | — | FIRE peaks BED(= `indir` 输入) |
| — | `-g/--gtf` | — | 基因注释 GTF |
| — | `-r/--repeat-bed` | 无 | repeat BED(可选) |
| — | `-o` / `--sample` | — | 输出目录 / 样本名前缀 |
| — | `-l` | — | 日志文件 |
| — | `-n/--dry-run` | off | 只打印计划不产出 |
| threads / mem_mb | — | 4 / 8000 | `rule_res("fiberseq_annotate")` 配置 |

直接调用示例(脱离流程):

```bash
python modules/fiberseq/bin/fire_annotate.py \
  -p S1.fire_peaks.bed -g gencode.gtf -r repeat.bed \
  -o out/S1 --sample S1 -t 1000 -f 0.25 -l fire_annotate.log
```

## 2. 算法原理与生物学意义

- **基因组分类**(优先级 promoter > exon > intron > intergenic):取 peak
  中点;若落在任一基因 TSS ± `tss_flank` 窗口内 → `promoter`(基因归属取
  距中点最近的 TSS);否则与 exon 区间相交 → `exon`;否则与 gene body
  相交 → `intron`;其余 → `intergenic`。无归属基因时按 TSS 距离补
  `nearest_gene`。
- **repeat 重叠**:peak 与 repeat BED 的重叠碱基数 / peak 长度 =
  `repeat_frac` 列;≥ `repeat_frac` 阈值标记 `repeat_flag=yes`,
  `repeat_names` 列出 repeat 名。
- **汇总**:按 region_class、repeat_flag 计数,附总 peak 数。
- **生物学意义**:FIRE peak 是 Fiber-seq 判定的开放染色质候选调控元件;
  注释回答"这些元件落在基因组什么功能上下文、是否重复序列衍生"——
  promoter 富集支持启动子调控解释,intergenic 多提示远端增强子候选,
  repeat 高重叠提示 TE 来源的调控序列或需要谨慎解读的区域。
- **结果怎么读**:`peaks_annotated.tsv` 逐 peak 一行;先看 `region_class`
  分布(`peaks_summary.tsv`),再挑 `repeat_flag=no` 的 promoter/intergenic
  peak 作为高置信调控候选。

## 3. 参考文献

1. FIRE analysis pipeline — https://github.com/fiberseq/FIRE
2. GTF/GFF 格式规范 —
   https://github.com/The-Sequence-Ontology/Specifications/blob/master/gtf3.md
3. 引擎脚本 `bin/fire_annotate.py` 为本仓库自研(共享 `fiberseq/bin/`,
   见父模块 `fiberseq/README.md`)。
