# mimseq 模块

tRNA 修饰诱导错配测序（mim-tRNAseq）分析模块，基于 [mimseq](https://github.com/drewjbehren/mimseq) 工具。

## 目录结构

```
modules/mimseq/
  mimseq.smk              # 主入口：include 6 个子模块 + prepare_sample_data + mimseq_all/result
  mimseq.yaml             # 共享 conda 环境
  tRNAtools/
    tRNAtools.smk         # tRNA 注释、聚类、SNP 索引
  align/
    align.smk             # GSNAP 比对（SNP-tolerant）
  clusters/
    clusters.smk          # isodecoder 去卷积
  mods/
    mods.smk              # 错配 / 修饰定量
  coverage/
    coverage.smk          # 覆盖度分析 + CCA
  deseq/
    deseq.smk             # DESeq2 差异表达
  bin/                    # 各子模块的 Python wrapper
    tRNAtools/run.py
    align/run.py
    clusters/run.py
    mods/run.py
    coverage/run.py
    deseq/run.py
    serialize.py          # pickle 状态序列化 / 反序列化工具
    utils.py              # 共享工具函数（如 extract_condition）
  mimseq/                 # mimseq Python 库源码（上游）
    tRNAtools.py
    tRNAmap.py
    splitClusters.py
    mmQuant.py
    getCoverage.py
    ssAlign.py
    crosstalks.py
    modifications/        # 修饰参考数据
    data/                 # 物种 tRNA 参考数据
```

## 执行顺序（DAG 依赖链）

```
prepare_sample_data
    ↓
tRNAtools (state/)  →  align (align.done)  →  clusters (clusters.done)
    ↓
mods (mods.done)  →  coverage (coverage.done)  →  deseq (deseq.done)
    ↓
mimseq.done (result)
```

## 子模块说明

### 1. prepare_sample_data
从 meta.tsv 生成 `sample_data.tsv`（mimseq 所需的 sample-condition 映射文件）。

### 2. tRNAtools
- 解析 tRNA 序列，生成聚类和 SNP 索引
- 支持内置物种（`species` 参数：Hsap, Mmus, Scer 等）和自定义 tRNA 参考
- 输出：`state/` 目录（pickle 序列化的状态文件）

### 3. align
- 使用 GSNAP 进行 SNP-tolerant 比对
- 支持 remap（两轮比对发现新修饰位点）
- 输出：`samples/{sample}/` 子目录下的 BAM 文件，`align.done`

### 4. clusters
- isodecoder 去卷积：将模糊比对拆分为唯一定位的 tRNA 转录本
- 参数：`cluster_id`（默认 0.97）、`cov_diff`（默认 0.5）
- 输出：`clusters.done`

### 5. mods
- 错配 / 修饰定量分析
- 参数：`min_cov`（最小覆盖度阈值）、`misinc_thresh`（错配阈值，默认 0.1）
- 支持 remap、crosstalks 分析
- 输出：`mods.done`

### 6. coverage
- 计算每个样本的 tRNA 覆盖度
- 输出文件：
  - `coverage_byaa.txt` — 按氨基酸聚合的覆盖度（TSV：aa, bin, condition, bam, pos, cov, cov_norm）
  - `coverage_bygene.txt` — 按基因聚合的覆盖度
  - `coverage.log` — 运行日志
  - `state/coverageData.pkl` — 覆盖度数据 pickle
- 支持 CCA 分析（`--cca`）和双 CCA（`--double-cca`）
- 输出：`coverage.done`

### 7. deseq
- DESeq2 差异 tRNA 表达分析
- 需要指定 `control_cond` 作为对照条件
- 参数：`p_adj`（校正 p 值阈值，默认 0.05）
- 输出：`deseq.done`

## 配置参数

### Params.mimseq

| 参数 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `species` | str | "" | 内置物种（Hsap/Mmus/Scer 等），为空则使用自定义 tRNA 文件 |
| `name` | str | "tRNAseq" | 分析名称前缀 |
| `control_cond` | str | "" | DESeq2 对照条件名 |
| `cluster_id` | float | 0.97 | 聚类相似度阈值 |
| `threads` | int | 8 | 线程数 |
| `min_cov` | float | 0.0005 | 最小覆盖度阈值 |
| `max_mismatches` | float | 0.075 | 最大错配率 |
| `max_multi` | int | 4 | 最大多重比对数 |
| `misinc_thresh` | float | 0.1 | 错配阈值 |
| `remap_mismatches` | float | 0.05 | remap 错配率 |
| `p_adj` | float | 0.05 | DESeq2 校正 p 值阈值 |
| `cov_diff` | float | 0.5 | 覆盖度差异阈值 |
| `no_cluster` | bool | false | 禁用聚类 |
| `no_cca` | bool | false | 禁用 CCA 分析 |
| `double_cca` | bool | false | 启用双 CCA |
| `remap` | bool | true | 启用两轮比对 |
| `snp_tolerance` | bool | true | 启用 SNP 容忍 |
| `keep_temp` | bool | false | 保留临时文件 |
| `crosstalks` | bool | false | 启用串扰分析 |
| `pretRNAs` | bool | false | 包含 pre-tRNA |
| `posttrans_mod_off` | bool | false | 关闭转录后修饰 |
| `local_modomics` | bool | false | 使用本地 Modomics 数据 |

### genome

| 参数 | 说明 |
|---|---|
| `trnas` | 自定义 tRNA FASTA（species 为空时必填） |
| `trnaout` | 自定义 tRNAscan-SE 输出（species 为空时必填） |
| `mito_trnas` | 线粒体 tRNA FASTA（可选） |
| `plastid_trnas` | 质体 tRNA FASTA（可选） |

## 注意事项

1. **conda 环境路径**：当前 smk 文件中硬编码了 Python 路径，需根据部署环境调整
2. **mimseq 输出目录**：路径必须有尾部 `/`，否则内部路径拼接会出错
3. **pickle 状态**：tRNAtools 通过 pickle 传递状态到下游模块，`bin/serialize.py` 提供 load/save 工具
4. **coverage_byaa.txt 用途**：可直接用于绘制氨基酸水平表达热力图（如 PPT 中的 slide_cn-4）
