# mimseq/clusters — isodecoder 去卷积子模块(mim-tRNAseq)

mimseq 家族的聚类/拆分步骤:把模糊比对的 tRNA 读段按序列相似度聚类、
拆分为唯一定位的 isodecoder 转录本。属 mimseq 家族(共享父目录
`mimseq.yaml`、`bin/` wrapper 与 `mimseq/` Python 库);家族结构、参数
总表与执行顺序见 [mimseq/README.md](../README.md)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `align_done`(规则 input) | `{outdir}/align.done` —— 上游 align 子模块(GSNAP 比对)的完成标记;真实数据由 mimseq 库按 `state/` pickle 与 `samples/{sample}/` 约定在共享 outdir 内传递 |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 分析根目录 / 日志 / 项目根(定位 `modules/mimseq/bin/clusters/run.py`) |

### 1.2 输出

```
{outdir}/clusters.done        # 完成标记(touch);分析中间产物由 mimseq 库写在 outdir 内
{outdir}/mimseq_clusters.sh   # 回执
```

日志:`{logdir}/mimseq_clusters.log`。`mimseq_clusters_result` 为聚合
规则。mimseq 输出目录须带尾部 `/`(家族约定,见父 README 注意事项)。

### 1.3 参数

config 位于 `Params.mimseq` 下(家族共享参数段):

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.mimseq.name` | `--name` | `tRNAseq` | 分析名称前缀 |
| `Params.mimseq.cluster_id` | `--cluster-id` | 0.97 | 聚类序列相似度阈值 |
| `Params.mimseq.cov_diff` | `--cov-diff` | 0.5 | 覆盖度差异阈值(拆分判据) |
| `Params.mimseq.cluster` | `--cluster` | true | false 时不传 `--cluster`(禁用聚类) |
| (threads) | `--threads` | 4 | `rule_res("mimseq_clusters")` 可覆盖 |
| (mem_mb) | — | 8000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python modules/mimseq/bin/clusters/run.py \
    --out out/ --name tRNAseq --threads 4 \
    --cluster-id 0.97 --cov-diff 0.5 --cluster
```

注:解释器为规则 conda 环境(`../mimseq.yaml`)的 `python`。

## 2. 算法原理与生物学意义

tRNA 基因家族高度同源,isodecoder 间仅差个别碱基,读段常模糊比对到
多个拷贝。clusters 步骤按 `cluster_id`(默认 0.97)把近乎相同的
tRNA 序列聚成簇,再依 `cov_diff` 覆盖度差异把模糊读段拆分回唯一定位
的转录本(去卷积),使后续修饰/表达定量落在有区分度的单元上。
生物学意义:mim-tRNAseq 用修饰诱导错配读出 tRNA 修饰与丰度;若不做
去卷积,同源拷贝的信号互相稀释,isodecoder 特异的修饰差异会被抹掉。
`--cluster=false`(或 `Params.mimseq.no_cluster`)保留原始粒度供对照。

## 3. 参考文献

1. mimseq(mim-tRNAseq 分析工具)—— https://github.com/drewjbehren/mimseq
   (家族内未附论文引文;参数与流程约定见 [mimseq/README.md](../README.md))
