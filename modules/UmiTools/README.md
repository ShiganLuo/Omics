# UmiTools — UMI 读段去重模块(UMI-tools dedup)

UMI-tools `dedup` 封装:按 UMI(唯一分子标识)把 PCR 重复折叠为每
分子一条读段,输出去重 BAM 与索引(临时文件,供下游计数/变异分析)。
同目录 `extract/` 子模块负责 UMI 提取(读段名附加 UMI),本模块不覆盖。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | BAM 目录(上游 config 的 `outdir`);两条规则对应两种上游布局:`umi_tools_dedup_for_hisat2` 取 `{indir}/{sample_id}.bam`(平铺,hisat2 后接),`umi_tools_dedup_for_star` 取 `{indir}/{sample_id}/{sample_id}.bam`(样本子目录,STAR 后接) |
| `outdir` / `logdir` (config) | 输出/日志目录 |

UMI 前提:读段名已带 UMI(由上游 UmiTools/extract 或测序仪写入),
UMI-tools 从 read name 末段解析 UMI。

### 1.2 输出

```
{outdir}/{sample_id}.dedup.bam        # 去重 BAM(temp,下游消费后可回收)
{outdir}/{sample_id}.dedup.bam.bai    # samtools index(temp)
```

回执:`{outdir}/{sample_id}/umi_tools_dedup_for_{hisat2|star}.sh`;日志
`{logdir}/{sample_id}/umi_tools_dedup_run.txt`。

### 1.3 参数

config 位于 `Params.umi_tools` 与 `Procedure` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.umi_tools.method` | `--method` | `unique` | UMI 分组算法(UMI-tools 取值:unique / percentile / adjacency / directional / cluster;官方推荐 directional) |
| `Procedure.umi_tools` | 引擎可执行 | `umi_tools` | — |
| `Procedure.samtools` | `samtools index` | `samtools` | 建索引 |
| (threads) | — | 2 | `rule_res` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
umi_tools dedup --method=unique \
    -I sample.bam -S sample.dedup.bam
samtools index sample.dedup.bam
```

## 2. 算法原理与生物学意义

PCR 扩增会使同一原始分子产生多条完全相同的读段,常规坐标去重
(markdup)在低起始量 RNA-seq/扩增子数据中会误删真实高表达分子。
UMI(随机分子标签)在建库时给每个原始分子打唯一条码:比对后
「UMI + 比对位置」相同的读段组视为同一分子的 PCR 子代,`dedup` 每组
只保留一条代表读段。`--method` 决定 UMI 组的归并算法:

- `unique`:仅完全相同 UMI 归组(最保守);
- `percentile` / `adjacency` / `directional` / `cluster`:允许把测序
  错误产生的邻居 UMI 折叠回真实 UMI,灵敏度递进(官方推荐
  `directional`)。

生物学意义:去重后的计数直接对应原始分子数,是 UMI 建库定量
(RNA 表达、低频突变)的基石;`dedup.bam` 标 temp,只作下游中间物。
注意 `--method=unique` 对 UMI 测序错误敏感,高错误率平台可考虑
directional。

## 3. 参考文献

1. UMI-tools — https://github.com/CGATOxford/UMI-tools
   (文档 https://umi-tools.readthedocs.io/)
2. Smith TS, Heger A, Sudbery I. *UMI-tools: Modelling sequencing errors
   in Unique Molecular Identifiers to improve quantification accuracy*.
   Genome Research, 2017. DOI: 10.1101/gr.209601.116.
