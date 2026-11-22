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
| wildcards `genomeA` / `genomeB` / `sample_id` | 两个基因组 ID 与样本 ID;genomeA/genomeB 顺序不固定,输出文件名内嵌基因组 ID,消费方按基因组名而非位置取文件 |
| `outdir` / `logdir` (config) | 输出/日志目录 |
| `Procedure.neodisambiguate` | 消歧引擎可执行路径 |

引擎契约(实测 neodisambiguate 1.1.1):两份输入 BAM 的模板名集合必须
一致(同一批数据比对到两个参考,真实场景天然满足),否则报
`Templates with different names found!` 退出。

### 1.2 输出

```
{outdir}/{sample_id}/{genomeA}_{genomeB}/
├── {sample_id}.neodisambiguated.{genomeA}.bam          # 基因组 A 干净 BAM
├── {sample_id}.neodisambiguated.{genomeB}.bam          # 基因组 B 干净 BAM
├── ambiguous-alignments/
│   ├── {sample_id}.{genomeA}.ambiguous.bam             # A 输入未定归属的记录
│   └── {sample_id}.{genomeB}.ambiguous.bam             # B 输入未定归属的记录
├── inputs/                                             # 供引擎读取的改名暂存软链
└── neodisambiguate.sh                                  # 回执
```

日志:`{logdir}/{sample_id}/{genomeA}_{genomeB}/neodisambiguate.log`。

暂存软链的必要性:引擎的 ambiguous 文件按"输入文件名 stem"命名
(`{stem}.ambiguous.bam`),两份输入同名(都是 `{sample_id}.bam`)时两个
ambiguous 文件互相覆盖、先写的一份被静默丢弃;故规则先把两份输入软链为
`inputs/{sample_id}.{genomeA}.bam` / `.{genomeB}.bam` 再喂给引擎,保证
ambiguous 输出按输入各一份。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.neodisambiguate` | 引擎可执行 | `neodisambiguate` | 回执以 `bash "$(command -v …)"` 调用:conda 包的 launcher 用 bash 数组语法却声明 `#!/usr/bin/env sh`,在 dash 为 sh 的容器里直接执行会 Syntax error |
| (inputs) | `-i` | 规则暂存的两份软链 | 待消歧的两份 BAM |
| (prefix) | `-o` | `{pair_dir}/{sample_id}.neodisambiguated` | 输出前缀(引擎输出 `{prefix}.{genome}.bam`) |
| (names) | `-n` | `genomeA genomeB` | 两个参考名,决定干净 BAM 后缀 |
| (threads) | — | 8 | `rule_res("neodisambiguate")` 可覆盖 |
| (mem_mb) | — | 16000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
bash "$(command -v neodisambiguate)" \
    -i outdir/sample1/genomeA_genomeB/inputs/sample1.genomeA.bam \
       outdir/sample1/genomeA_genomeB/inputs/sample1.genomeB.bam \
    -o outdir/sample1/genomeA_genomeB/sample1.neodisambiguated \
    -n genomeA genomeB
```

## 2. 算法原理与生物学意义

当同一份测序数据可能来自两个基因组(宿主 vs 移植/共生、人 vs 鼠
异种移植、两个物种共培养)时,分别比对会产生一批在两个基因组上都
比对得上的歧义读段。消歧引擎按两份比对的比对质量/唯一性把每个读段
判给更可信的基因组或判为歧义:

- **clean BAM**(`.neodisambiguated.{genome}.bam`):归属明确的读段,
  供各基因组下游分析(定量、变异检测)使用,避免跨基因组污染;
- **ambiguous BAM**:该输入中未能明确判归的读段(按输入各一份),
  留作审计与评估污染程度;若占比过高提示参考基因组选择或样本纯度
  有问题。

输出命名内嵌 genomeA/genomeB 基因组 ID(wildcards 顺序不固定),消费方
按文件名中的基因组名取文件,不依赖位置。

## 3. 参考文献

1. `neodisambiguate` 引擎(fulcrumgenomics Disambiguate 1.1.1,经
   `Procedure.neodisambiguate` 注入)——同仓库另有同类封装
   `disambiguate/`(AstraZeneca Disambiguate)与 `XenofilteR/` 可参考。
