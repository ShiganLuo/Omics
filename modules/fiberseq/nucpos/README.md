# fiberseq/nucpos — TSS 周围核小体定位/相位分析子模块

以 TSS(或 BED 锚点)为中心收集单纤维核小体中心,按「锚点上/下游第
i 个核小体」分槽,统计每槽平均偏移(相位)与离散度衰减,出表 + PDF。
属 fiberseq 家族(布局 C),共享父目录 `../fiberseq.yaml` 环境与
`bin/` 脚本;家族输入约定见 [fiberseq/README.md](../README.md)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | 元素目录(默认 `common/4_elements`,家族唯一权威布局):`{indir}/{sample_id}/{sample_id}.nuc.bed.gz`(核小体 BED12,bgzip) |
| `genome.gtf` (config) | 基因注释,取每 gene 的 TSS 作锚点(**必填**,缺失报错) |
| `samples` (config) | 样本 ID 列表 |
| `outdir` / `logdir` / `ROOT_DIR` (config) | 输出(默认 `fiberseq_analysis/nucpos`)/日志/项目根 |

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.nuc_offset_per_anchor.tsv   # 每锚点 × 槽位偏移
{outdir}/{sample_id}/{sample_id}.nuc_offset_summary.tsv      # 槽位汇总(均值/离散度)
{outdir}/{sample_id}/{sample_id}.nucpos.pdf                  # 相位图
```

回执:`{outdir}/{sample_id}/fiberseq_nucpos_run_{sample_id}.sh`(另有
`{logdir}/cmd/...` 时间戳历史,见父 README);日志
`{logdir}/{sample_id}/fiberseq_nucpos_run.log`。`fiberseq_nucpos_result`
聚合全部样本 summary。

### 1.3 参数

config 位于 `Params.fiberseq_nucpos` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `genome.gtf` | `-g` | — | TSS 锚点来源,必填 |
| (anchor) | `--anchor` | `tss` | 锚点语义(tss/bed,规则固定 tss) |
| `indir` | `-u` | `common/4_elements` | 核小体 BED12 |
| `Params.fiberseq_nucpos.upstream_nucs` | `-5` | 5 | 上游核小体槽数 |
| `Params.fiberseq_nucpos.downstream_nucs` | `-3` | 5 | 下游核小体槽数 |
| `Params.fiberseq_nucpos.min_fibers` | `-i` | 10 | 每锚点最少纤维数,不足剔除 |
| (threads) | — | 4 | `rule_res("fiberseq_nucpos")` 可覆盖 |
| (mem_mb) | — | 6000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python modules/fiberseq/bin/nucpos.py \
    -g genes.gtf --anchor tss \
    -u common/4_elements/s1/s1.nuc.bed.gz \
    -o fiberseq_analysis/nucpos/s1 --sample s1 \
    -5 5 -3 5 -i 10 -l log/s1/fiberseq_nucpos_run.log
```

## 2. 算法原理与生物学意义

`nucpos.py` 沿 Cell Rep Methods 2024 的核小体定位分析:

1. 从 GTF 提取每 gene 的 TSS(链方向校正)作锚点;
2. 收集窗口内每条纤维的核小体中心,把每个核小体指派到相对槽位
   (锚点上/下游第 i 个最近核小体,`-5`/`-3` 控制范围);
3. 逐槽统计平均偏移(相对理论间隔的相位偏移)与离散度:转录活跃
   基因 TSS 呈典型「+1 核小体定位 + 上游相位梯」,相位随距离衰减;
4. `min_fibers` 过滤低覆盖锚点,保证偏移估计可靠。

结果怎么读:`per_anchor` 表逐锚点给各槽偏移(可做基因级筛选);
`summary` 表逐槽给均值/离散度(画相位曲线);PDF 可视化 +1/-1 相位
结构。生物学意义:核小体定位是转录可及性的物理基础——TSS 周围
相位清晰度反映染色质状态(活跃基因 +1 紧凑、上游有序;沉默基因
相位模糊),样本间对比可观察衰老/处理对染色质组织的影响。

## 3. 参考文献

1. Peter CJ et al. *Single chromatin fiber profiling and nucleosome
   position mapping in the human brain*. Cell Reports Methods 4(12):100911,
   2024. DOI: 10.1016/j.crmeth.2024.100911(方法出处)。
2. 家族其余方法来源见 [fiberseq/README.md](../README.md) Sources:
   DOI: 10.1101/gr.279095.124、DOI: 10.1038/s41588-025-02324-w。
