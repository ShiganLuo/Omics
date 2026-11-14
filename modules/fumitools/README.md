# fumitools — UMI 提取与 UMI 去重模块(fumi_tools)

把 read 中的 UMI 拷入 FASTQ 标题(`copy_umi`,双端/单端)并按
UMI + 比对位置去重(`dedup`);`bin/detect_umi.py` 检测标题中是否已有
UMI,已有则符号链接直通。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 规则 | 说明 |
|---|---|---|
| `{indir}/{sample_id}_1.fq.gz` + `_2.fq.gz` | fumitools_copy_umi_paired | 双端 FASTQ(indir 平铺,无样本子目录) |
| `{indir}/{sample_id}.single.fq.gz` | fumitools_copy_umi_single | 单端 FASTQ |
| `{outdir}/{sample_id}/{sample_id}.sorted.bam` | fumitools_dedup | 坐标排序 BAM(上游比对产物) |

`detect_umi.has_umi_in_header` 检查 R1 首条 read 名:`_<ACGTN>+$` 后缀即
认为 UMI 已提取(umi_tools extract / fumi_tools copy_umi 的命名约定),
此时输出直接 symlink 输入、跳过提取。

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}_1.umi.fq.gz        # UMI 拷入标题的 R1(或 symlink)
{outdir}/{sample_id}/{sample_id}_2.umi.fq.gz        # UMI 拷入标题的 R2
{outdir}/{sample_id}/{sample_id}.umi.single.fq.gz   # 单端版本
{outdir}/{sample_id}/{sample_id}.dedup.bam          # UMI 去重 BAM
{outdir}/{sample_id}/fumitools_copy_umi.sh          # 回执(fumitools_dedup.sh 同目录)
logdir/{sample_id}/fumitools_copy_umi.log           # 运行日志(dedup 为 fumitools_dedup.log)
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.fumitools` | 可执行文件 | null(回退 `fumi_tools`) | fumi_tools 路径 |
| `Params.fumitools.umi_length` | `--umi-length` | 12 | UMI 长度;copy_umi 实际提取时必填(缺失抛 `ValueError`) |
| `Params.fumitools.tag_umi` | `--tag-umi` | false | 以标签形式附加 UMI |
| `Params.fumitools.start_only` | `--start-only` | false | 只按 read 起点判重(dedup) |
| `Params.fumitools.memory` | `--memory` | `3G` | dedup 内存上限 |
| `Params.fumitools.threads` | `--threads` | 4 | dedup 线程(该规则不走 rule_res) |
| `paired_samples` | `--paired` | [] | 样本在列表中时 dedup 加 `--paired` |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录基路径 |

copy_umi 规则线程经 `rule_res(...)` 配置,默认 4。

引擎命令(回执脚本原文):

```bash
fumi_tools copy_umi -i R1.fq.gz -I R2.fq.gz -o {s}_1.umi.fq.gz -O {s}_2.umi.fq.gz --threads 4 --umi-length 12
fumi_tools dedup -i {s}.sorted.bam -o {s}.dedup.bam --threads 4 --memory 3G [ --paired ] [ --start-only ]
```

## 2. 算法原理与生物学意义

纯工具封装:`copy_umi` 把 read 内的 UMI(建库时加在 read 头部的随机
分子标签)切下并以 `_` 连接符追加到 read 名(如
`@...:1407:1064_ANGGCTTCCCTG`),供后续按分子归并;`dedup` 把"UMI 相同
且比对位置相同(可选仅起点)"的 reads 视为同一原始分子的 PCR 拷贝,只留
一条。生物学意义:UMI 去重消除扩增偏好,让计数回到原始分子层面,是
低起始量/单细胞/3' 端计数类文库定量准确性的关键;标题检测避免重复提取
破坏已有 UMI。

## 3. 参考文献

1. fumi_tools 0.19.0(conda channel `ccb-sb`,见 `fumitools.yaml`)
