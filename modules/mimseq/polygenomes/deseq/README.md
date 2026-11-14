# mimseq/polygenomes/deseq — 多基因组 tRNA 差异表达子模块(DESeq2)

mimseq 家族的差异表达步骤(按基因组分目录的 polygenomes 变体):
对 coverage 阶段的 tRNA 计数做 DESeq2 分析,输出差异 tRNA 表。
共享父目录 `mimseq.yaml`、`bin/deseq/run.py`(与非 polygenomes 的
`mimseq/deseq/` 同一脚本)与 `mimseq/` Python 库;家族约定见
[mimseq/README.md](../../README.md)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `coverage_done`(规则 input) | `{outdir}/{genome}/coverage.done` —— 同基因组 coverage 子模块的完成标记;计数数据由 mimseq 库按 outdir 内约定(`state/coverageData.pkl` 等)传递 |
| wildcards `genome` | 基因组 ID(polygenomes 布局,输出按基因组分目录) |
| `genome.<genome_id>.mito_trnas` (config) | 线粒体 tRNA FASTA(可选) |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 分析根 / 日志 / 项目根 |

### 1.2 输出

```
{outdir}/{genome}/deseq.done        # 完成标记(touch);差异表由 mimseq 库写在同目录
{outdir}/{genome}/mimseq_deseq.sh   # 回执
```

日志:`{logdir}/{genome}/mimseq_deseq.log`。`mimseq_deseq_result` 为聚合
规则。输出目录须带尾部 `/`(家族约定)。

### 1.3 参数

config 位于 `Params.mimseq`(家族共享)与 `genome.<g>` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.mimseq.control_cond` | `--control-cond` | "" | DESeq2 对照条件名;非空才传(建议必设,否则对比方向不明) |
| `Params.mimseq.p_adj` | `--p-adj` | 0.05 | 校正 p 值(FDR)阈值 |
| `genome.<g>.mito_trnas` | `--mito-trnas` | "" | 线粒体 tRNA FASTA;非空才传(线粒体 tRNA 单独处理) |
| (threads) | — | 4 | `rule_res("mimseq_deseq")` 可覆盖 |
| (mem_mb) | — | 16000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python modules/mimseq/bin/deseq/run.py \
    --out out/genome1/ --p-adj 0.05 \
    --control-cond control --mito-trnas mito_trnas.fa
```

注:解释器为规则 conda 环境(`../../mimseq.yaml`)的 `python`。

## 2. 算法原理与生物学意义

DESeq2 用负二项模型 + 收缩估计对 tRNA 计数做差异表达:以
`control_cond` 为基线估计每个 tRNA(或聚类单元)的 log2FC 与 Wald
检验 p 值,BH 校正后按 `p_adj` 阈值筛显著差异。`--mito-trnas` 把
线粒体 tRNA 与细胞质 tRNA 分开建模(丰度与修饰背景差异大,混在一起
会稀释信号)。polygenomes 变体按 `{genome}` 分目录独立分析:同一读段
在不同参考(物种/组装)下的 tRNA 注释不同,差异结果须按基因组分别
解读,不可跨基因组合并。输出差异 tRNA 表(由 mimseq 库落盘)供
下游修饰-表达关联分析。

## 3. 参考文献

1. mimseq(mim-tRNAseq 分析工具)—— https://github.com/drewjbehren/mimseq
2. DESeq2 —— https://github.com/thelovelab/DESeq2
   (Love MI, Huber W, Anders S. *Moderated estimation of fold change and
   dispersion for RNA-seq data with DESeq2*. Genome Biology 15:550, 2014.
   DOI: 10.1186/s13059-014-0550-8.)
