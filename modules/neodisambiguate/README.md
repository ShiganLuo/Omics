# neodisambiguate — 双基因组比对歧义读段消歧模块

把同一批读段分别比对到两个基因组(如宿主/共生、人/鼠、两个单倍型)后
得到的两份 BAM 做消歧:输出各自基因组的干净 BAM 与歧义比对 BAM。
纯工具封装模块(引擎 `neodisambiguate` 由 `Procedure.neodisambiguate`
注入,模块内无 bin/ 脚本)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | 双基因组布局 BAM 目录(上游 config 的 `outdir`):`{indir}/{genomeA}/{sample_id}.bam` 与 `{indir}/{genomeB}/{sample_id}.bam` |
| wildcards `genomeA` / `genomeB` / `sample_id` | 两个基因组 ID 与样本 ID;genomeA/genomeB 顺序不固定,输出以 A/B 后缀区分 |
| `outdir` / `logdir` (config) | 输出/日志目录 |
| `Procedure.neodisambiguate` | 消歧引擎可执行路径 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.{genomeA}.neodisambiguatedA.bam          # 基因组 A 干净 BAM
{outdir}/{sample_id}/{sample_id}.{genomeB}.neodisambiguatedB.bam          # 基因组 B 干净 BAM
{outdir}/{sample_id}/ambiguous-alignments/{sample_id}.{genomeA}.ambiguousA.bam   # 歧义比对(A 参考)
{outdir}/{sample_id}/ambiguous-alignments/{sample_id}.{genomeB}.ambiguousB.bam   # 歧义比对(B 参考)
{outdir}/{sample_id}/neodisambiguate.sh                                   # 回执
```

日志:`{logdir}/{sample_id}/neodisambiguate.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.neodisambiguate` | 引擎可执行 | `neodisambiguate` | — |
| (prefix) | `-s` | `{outdir}/{sample_id}.neodisambiguated` | 输出前缀(由规则拼接) |
| (outdir) | `-o` | config `outdir` | 输出目录 |
| (threads) | — | 8 | `rule_res("neodisambiguate")` 可覆盖 |
| (mem_mb) | — | 16000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
neodisambiguate \
    -s outdir/sample1.neodisambiguated \
    -o outdir \
    input/genomeA/sample1.bam \
    input/genomeB/sample1.bam
```

## 2. 算法原理与生物学意义

当同一份测序数据可能来自两个基因组(宿主 vs 移植/共生、人 vs 鼠
异种移植、两个物种共培养)时,分别比对会产生一批在两个基因组上都
比对得上的歧义读段。消歧引擎按两份比对的比对质量/唯一性把每个读段
判给更可信的基因组或判为歧义:

- **clean BAM**(neodisambiguatedA/B):归属明确的读段,供各基因组
  下游分析(定量、变异检测)使用,避免跨基因组污染;
- **ambiguous BAM**:两个基因组都无法明确判别的读段,留作审计与
  评估污染程度;若占比过高提示参考基因组选择或样本纯度有问题。

输出命名中 A/B 后缀与 genomeA/genomeB wildcards 对应(wildcards 顺序
不固定),消费方按后缀而非位置取文件。

## 3. 参考文献

1. `neodisambiguate` 引擎(经 `Procedure.neodisambiguate` 注入)——
   模块文件内未附公开文献;同仓库另有同类封装 `disambiguate/`
   (AstraZeneca Disambiguate)与 `XenofilteR/` 可参考。
