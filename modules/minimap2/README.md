# minimap2 — 长读长比对模块

minimap2 长读长比对(`-ax` 预设)+ samtools 转 BAM/排序/建索引,输出坐标
排序 BAM。适用于 ONT/PacBio 等长读长单端数据(也接受双端 R1)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `genome.fasta`(config,必填) | 参考基因组 FASTA;缺失时 parse 期直接抛 `ValueError` |
| FASTQ(形态 A:indir + 布局参数) | 模块 input 函数 `get_minimap2_input` 动态解析:优先 `{indir}/{sample_id}/{sample_id}{fastq_suffix}`(默认后缀 `.single.fq.gz`);若不存在但存在 `{indir}/{sample_id}/{sample_id}_1.fq.gz` 则取双端 R1;否则回退单端命名(长读长默认) |
| `samples` | 样本 ID 列表(供 subworkflow 注册驱动 DAG) |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.sorted.bam       # 坐标排序 BAM
{outdir}/{sample_id}/{sample_id}.sorted.bam.bai   # 索引
{outdir}/{sample_id}/minimap2_align_{sample_id}.sh  # 回执(整条流水命令)
{logdir}/{sample_id}/minimap2_align.log           # 日志
```

中间产物 `{sample_id}.raw.bam(.bai)` 由回执脚本 `mv` 为最终的
`.sorted.bam(.bai)`。收尾规则 `minimap2_result` 聚合 BAM + BAI。

### 1.3 参数

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.minimap2` | 可执行文件本身 | `minimap2` | 工具可执行路径 |
| `Procedure.samtools` | 可执行文件本身 | `samtools` | samtools 可执行路径 |
| `Params.minimap2.preset` | `-ax` | `map-ont` | minimap2 预设(ONT 噪声长读长);PacBio HiFi 用 `map-hifi`,CLR 用 `map-pb` |
| `Params.minimap2.fastq_suffix` | — | `.single.fq.gz` | 输入 FASTQ 后缀(形态 A 布局参数) |
| `Params.minimap2.extra` | — | `""` | 预留的额外参数字段,当前未拼入命令行 |
| `genome.fasta` | 参考序列参数 | — | 参考基因组 FASTA,必填 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 输入/输出/日志根目录 |

规则默认资源:`threads = 16`、`mem_mb = 32000`(经 `rule_res("minimap2_align")`
可被流程资源配置覆盖)。命令固定附加 `--eqx -L`、`-t N`。

直接调用示例(模块无 bin/ 引擎脚本;规则拼出的回执命令形态):

```bash
minimap2 -ax map-ont --eqx -L -t 16 ref.fa sample.single.fq.gz \
  | samtools view -@ 4 -bS - \
  | samtools sort -@ 4 -o sample.raw.bam
samtools index -@ 4 sample.raw.bam
mv sample.raw.bam  output/{sample_id}/{sample_id}.sorted.bam
mv sample.raw.bam.bai output/{sample_id}/{sample_id}.sorted.bam.bai
```

## 2. 算法原理与生物学意义

minimap2 是基于 minimizer 的长读长比对器:对参考与 query 各取窗口最小
k-mer(minimizer)作种子,动态规划链接(chain)成链,再做碱基级比对
(带缺口扩展),输出 PAF/SAM。本模块固定:

- `-ax <preset>`:按数据类型套用推荐参数组合(`map-ont` 面向 ~10% 错误率
  的 ONT 噪声长读长;`map-hifi` 面向 HiFi;`map-pb` 面向 CLR);
- `--eqx`:CIGAR 中以 `=`/`X` 显式区分匹配/错配位点,便于下游变异与
  单倍型分析直接读取比对差异;
- `-L`:把超过 65535 个操作的 CIGAR 写入 CG 标签,保证长读长超长比对能被
  旧版 BAM 工具正确转换(长读长常规需要);
- 比对结果经 `samtools view -bS` 转 BAM、`samtools sort` 坐标排序、
  `samtools index` 建索引,成为下游分析通用形态。

生物学意义:长读长读段跨越重复区与结构变异断点,是基因组组装、结构变异、
单倍型分析的基础数据;比对质量直接决定下游结论。结果阅读:`samtools
flagstat`/`idxstats` 查比对率与参考覆盖;CIGAR 中 `=`/`X` 直接可读出匹配
与错配分布;日志与回执 `.sh` 保留完整命令可复现。

## 3. 参考文献

1. minimap2 —— https://github.com/lh3/minimap2;man 页
   https://lh3.github.io/minimap2/minimap2.html(`-L`/`--eqx`/preset 语义
   出处)
2. Li H. *Minimap2: pairwise alignment for nucleotide sequences*.
   Bioinformatics 2018, 34(18):3094–3100.
   DOI: 10.1093/bioinformatics/bty191
3. samtools —— https://www.htslib.org/
