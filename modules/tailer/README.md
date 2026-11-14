# tailer — RNA 3' 端尾部分析模块(Tailer)

Tailer 从比对 BAM 提取 RNA 3' 端尾部信息(如 poly(A) 尾),输出逐样本
`{sample_id}_tail.csv`。目录下另有 `polygenomes/` 子模块(多参考基因组
场景),不在本文档范围;conda 环境 `tailer.yaml`(容器同源)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` + 样本 ID 布局(形态 A) | `{indir}/{sample_id}/{sample_id}.bam` + `.bam.bai`:坐标排序 BAM 及索引 |
| `genome.gtf` | 基因注释 GTF(Tailer `-a`,尾部读段的基因归属) |
| `paired_samples` / `single_samples` | 样本列表(config 声明,规则本身不消费;供 subworkflow 注册) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_tail.csv    # 尾部分析结果
{logdir}/{sample_id}/tailer.log
```

注意:Tailer 把 `{sample_id}_tail.csv` 写在**输入 BAM 同目录**,规则随后
`shutil.copy2` 到输出路径(拷贝前先清理旧文件;若工具未生成该文件则抛
`FileNotFoundError`)。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.tailer` / `Procedure.Tailer` | 可执行文件本身 | `Tailer` | Tailer 可执行路径(两个键名都接受) |
| `Params.tailer.read` | `-read` | 2 | 分析的 read 端(双端数据取哪一端) |
| `Params.tailer.threshold` | `-t` | 100 | 尾部判定阈值 |
| `Params.tailer.rev_comp` | `--rev_comp` | false | 尾部序列取反向互补 |
| `genome.gtf` | `-a` | — | 基因注释 GTF |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

资源默认:`threads = 1`、`mem_mb = 8000`(经 `rule_res("tailer_global")`
可覆盖)。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的命令形态):

```bash
Tailer -a annotation.gtf -read 2 -t 100 input/s1/s1.bam
# 产物 s1_tail.csv(写在 BAM 同目录)被拷贝为 output/s1/s1_tail.csv
```

## 2. 算法原理与生物学意义

Tailer 解析比对 reads 的 3' 端,识别非模板加尾(poly(A)、oligo(U) 等),
按 GTF 注释把 read 归属到基因,输出尾长/尾型统计 CSV。生物学意义:mRNA
poly(A) 尾长决定稳定性与翻译效率,small RNA 的尿嘧啶化/腺苷化调控其稳定
性与功能;尾部谱是转录组后调控的重要观测(如与 star_3pass_gene 联用做
small RNA 尾部)。结果阅读:`{sample_id}_tail.csv` 按基因/尾型给出计数与
尾长分布;`-t` 阈值控制尾部判定灵敏度。

## 3. 参考文献

1. Tailer 尾部分析工具(经 `tailer.yaml` conda 环境安装;工具主页见
   环境内 `Tailer` 所属软件包说明)
