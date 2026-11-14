# targets — 捕获区段 BED 推导模块(自研 derive_targets.py)

从比对 BAM 的实际覆盖度反推杂交捕获/外显子组的靶区 BED,一份区间表
同时供 GATK `-L` 与 CNVkit `--targets` 消费。厂商 panel BED 常缺失或
停留在旧组装版本,而从当前参考上的真实读段重建可保证 build 一致。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `indir` (config) | BAM 目录(上游 config 的 `outdir`),形态 A 布局:`{indir}/{sid}/{sid}{bam_suffix}.bam`(全部样本) |
| `bam_suffix` (config) | BAM 文件名后缀段(如 `.sorted_markdup` 前的中缀);缺省空 |
| `samples` (config) | 样本 ID 列表 |
| `ROOT_DIR` (config) | 项目根目录,用于定位 `modules/targets/bin/derive_targets.py` |
| `outdir` / `logdir` (config) | 输出/日志目录 |

环境借用 samtools 模块:`conda: "../samtools/samtools.yaml"`、容器
`sif("../samtools/samtools.yaml", ...)`(`samtools depth` 为引擎依赖)。

### 1.2 输出

```
{outdir}/targets.bed              # 排序合并后的捕获区段 BED(chrom/start/end)
{outdir}/derive_targets.sh        # 回执
```

日志:`{logdir}/targets/derive_targets.log`。

### 1.3 参数

config 位于 `Params.derive` 下:

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.derive.min_depth` | `--min-depth` | 10 | 判定碱基被"捕获"的最小深度 |
| `Params.derive.min_bq` | `--min-bq` | 20 | `samtools depth -q` 最小碱基质量 |
| `Params.derive.min_mq` | `--min-mq` | 20 | `samtools depth -Q` 最小比对质量 |
| `Params.derive.min_gap` | `--min-gap` | 50 | 捕获区内小于该间隔的空隙合并,避免外显子内小间隙碎裂区间 |
| `Params.derive.min_region` | `--min-region` | 50 | 丢弃短于该长度的最终区间 |
| (`--fraction`) | `--fraction` | 1.0 | 脚本支持(须覆盖样本比例,1.0 = 全样本交集),.smk 未暴露,用回执脚本手改可调 |
| (threads) | — | 4 | `rule_res("derive_targets")` 可覆盖 |
| (mem_mb) | — | 16000 | 资源声明 |

直接调用示例(引擎命令,与规则生成的回执一致):

```bash
python modules/targets/bin/derive_targets.py \
    -b s1.sorted_markdup.bam s2.sorted_markdup.bam \
    -o out/targets.bed \
    --min-depth 10 --min-bq 20 --min-mq 20 \
    --min-gap 50 --min-region 50
```

## 2. 算法原理与生物学意义

脚本流程(`derive_targets.py`):

1. **逐碱基深度**:`samtools depth -a -q <min_bq> -Q <min_mq>
   -G UNMAP,SECONDARY,QCFAIL,DUP` 流式统计每样本每碱基深度
   (排除未比对/次要比对/QC 失败/重复读段);
2. **单样本捕获区**:深度 ≥ `min_depth` 的连续碱基段成区(游程切分);
3. **多样本共识**:事件扫描线统计每区间被多少样本覆盖,保留覆盖数
   ≥ round(fraction × 样本数) 的区段(fraction=1.0 即全体交集);
4. **合并与过滤**:间隔 ≤ `min_gap` 的相邻区段合并(捕获效率波动、
   外显子内小间隙不致碎裂),再丢弃 < `min_region` 的短区间;
5. **输出**:排序 BED(`chrom<TAB>start<TAB>end`),原子写入(tmp+rename);
   结果为空时退出码 1 并提示降低 `--min-depth`。

生物学意义:靶区 BED 是外显子/杂交捕获分析的坐标基准——GATK 变异
调用 `-L` 限定活性区域提速并压制 off-target 假阳性,CNVkit `--targets`
分箱计数。用真实覆盖反推的区间天然匹配当前参考与本批文库,规避厂商
BED 的版本错配;但注意它反映的是"这批样本实际捕获到的区域",极端
覆盖不均的样本会拉低共识区,多样本合并(而非单样本)是关键。

## 3. 参考文献

1. derive_targets.py 为本仓库自研脚本(见脚本 docstring),无外部文献;
   深度统计依赖 samtools —— https://github.com/samtools/samtools
