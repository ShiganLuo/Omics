# seqprep — 双端 reads 合并模块(SeqPrep merge)

封装 SeqPrep:对双端 FASTQ 做重叠区合并(merge),输出合并 reads、
未合并双端、被拒绝 reads 及统计报告。单规则 `seqprep_merge` + 收尾
`seqprep_result`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir`(config) | FASTQ 目录;input 函数拼 `{indir}/{sample_id}/{sample_id}_1.fq.gz` 与 `_2.fq.gz`(双端) |
| `samples` | 样本 ID 列表 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_merge.fastq.gz      # 合并 reads(SeqPrep -s)
{outdir}/{sample_id}/{sample_id}_unmerge_1.fastq.gz  # 未合并 R1(-1)
{outdir}/{sample_id}/{sample_id}_unmerge_2.fastq.gz  # 未合并 R2(-2)
{outdir}/{sample_id}/{sample_id}_reject_1.fastq.gz   # 被拒 R1(-3)
{outdir}/{sample_id}/{sample_id}_reject_2.fastq.gz   # 被拒 R2(-4)
{outdir}/{sample_id}/{sample_id}_merge.report        # seqkit stats 报告
```

收尾规则 `seqprep_result` 聚合 merge.fastq.gz + report。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.SeqPrep` | — | `SeqPrep` | SeqPrep 可执行路径(注意 config 键为 `SeqPrep`,大小写敏感) |
| `Procedure.seqkit` | — | `seqkit` | 统计工具 |
| — | `SeqPrep -f/-r/-1/-2/-3/-4/-s` | 固定 | 双端入、四路出 + 合并出 |
| — | `seqkit stats <merge>` | 固定 | 合并 reads 统计(写 `.report`) |
| mem_mb | — | 8000 | `rule_res("seqprep_merge")` 配置 |

直接调用示例(脱离流程):

```bash
SeqPrep -f S1_1.fq.gz -r S1_2.fq.gz \
  -1 S1_unmerge_1.fastq.gz -2 S1_unmerge_2.fastq.gz \
  -3 S1_reject_1.fastq.gz -4 S1_reject_2.fastq.gz \
  -s S1_merge.fastq.gz
seqkit stats S1_merge.fastq.gz > S1_merge.report
```

## 2. 算法原理与生物学意义

SeqPrep 对双端 reads 在重叠区做序列比对与共识合并(顺带质量校正与
接头修剪),插入片段短于 2×read 长时两条 read 可拼成一条更长的
consensus;不能合并的保留为 unmerge 双端,质量过差的进 reject。合并
后单条 read 覆盖整个插入片段,提高短片段(如 cfDNA/小 RNA/ChIP 短峰)
的比对与变异灵敏度。`seqkit stats` 给出合并库的条数/长度统计。纯工具
封装。

## 3. 参考文献

1. SeqPrep — https://github.com/jstjohn/SeqPrep
2. seqkit — https://bioinf.shenwei.me/seqkit/
