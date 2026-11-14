# mimseq/deseq — mim-tRNAseq 差异表达子模块(DESeq2)

mim-tRNAseq 流程的差异表达步:基于上游 coverage 定量结果做 DESeq2
差异分析(isodecoder/基因水平 tRNA 丰度)。引擎为共享父目录包装脚本
`modules/mimseq/bin/deseq/run.py`(内部执行 vendored mimseq 包的
`deseq.R`),conda 环境共享父目录 `../mimseq.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `coverage.done`(config outdir) | 上游 coverage 步的标记文件 `{outdir}/coverage.done`(本模块唯一的 DAG 锚点;定量表在同目录由上游产出) |
| `genome.mito_trnas` | 可选:线粒体 tRNA 注释(转交 `--mito-trnas`) |

### 1.2 输出

```
{outdir}/deseq.done      # touch 标记(DESeq2 结果表由引擎脚本写在 outdir 下)
{outdir}/mimseq_deseq.sh # 命令回执
```

收尾规则 `mimseq_deseq_result` 消费 `deseq.done`。

### 1.3 参数

config 位于 `Params.mimseq` 下,与 `deseq/run.py` CLI 对应:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `control_cond` | `--control-cond` | "" | 对照条件名(实验设计的参照水平;模板 JSON 键名拼作 `control_condition`,以 .smk 读取的 `control_cond` 为准) |
| `p_adj` | `--p-adj` | 0.05 | 显著性校正后 p 值阈值 |
| `genome.mito_trnas` | `--mito-trnas` | "" | 线粒体 tRNA 注释(可选) |
| threads / mem_mb | — | 4 / 16000 | `rule_res("mimseq_deseq")` 配置 |

注:解释器为规则 conda 环境(`../mimseq.yaml`)的 `python`。

直接调用示例(脱离流程):

```bash
python modules/mimseq/bin/deseq/run.py \
  --out out/ --p-adj 0.05 --control-cond treated
```

## 2. 算法原理与生物学意义

DESeq2 用负二项分布建模计数、中位数比值归一化,size factor + 离散度
收缩后做 Wald/LRT 检验,BH 校正得 padj;`--control-cond` 指定参照条件,
`--p-adj` 判显著。mim-tRNAseq 的 tRNA 定量达到 isodecoder 分辨率,差异
分析可揭示条件间(处理/对照)tRNA 丰度重编程——与翻译需求、修饰变化
及疾病状态相关。结果读法:mimseq 输出的 DESeq2 表按 `padj < p_adj` 与
log2FC 筛选;`deseq.done` 仅是 DAG 标记,真结果在 `{outdir}/` 下由
引擎脚本落盘。本模块为流程封装。

## 3. 参考文献

1. mim-tRNAseq — https://github.com/nedialkova-lab/mim-tRNAseq
   (vendored 于 `modules/mimseq/mimseq/`,v1.3.6)。
2. Love MI, Huber W, Anders S. *Moderated estimation of fold change and
   dispersion for RNA-seq data with DESeq2*. Genome Biology 15:550, 2014.
   DOI: 10.1186/s13059-014-0550-8
