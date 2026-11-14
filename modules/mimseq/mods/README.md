# mimseq/mods — mim-tRNAseq 修饰/错配掺入定量子模块

mim-tRNAseq 流程的 mods 步:对 tRNA 比对结果做修饰诱导错配
(misincorporation)定量。C 形态子模块(单规则 `mimseq_mods` + result 锚点,
`touch` 标记文件驱动);共享父目录脚本 `modules/mimseq/bin/mods/run.py`
(经 `config["ROOT_DIR"]` 拼接)与 `Params.mimseq` 参数段;流程全貌见父目录
`modules/mimseq/README.md`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{outdir}/clusters.done` | 上游 mimseq clusters 阶段的标记文件(唯一的 DAG 输入锚点) |
| `Params.mimseq.name` 等 | 运行标识与定量参数(见下表) |

### 1.2 输出

```
{outdir}/mods.done        # touch 标记文件(rule output,驱动 DAG)
{outdir}/                 # 修饰/错配掺入定量产物由 run.py 写入(文件名随 mimseq 包)
{logdir}/mimseq_mods.log
```

### 1.3 参数

本规则实际消费的 `Params.mimseq` 字段:

| config 字段 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.mimseq.name` | `--name` | `tRNAseq` | 运行名(产物命名前缀) |
| `Params.mimseq.min_cov` | `--min-cov` | 0.0005 | 位点最小覆盖占比,低于则不判修饰 |
| `Params.mimseq.misinc_thresh` | `--misinc-thresh` | 0.1 | 错配掺入率阈值 |
| `Params.mimseq.no_cca` | `--cca` | false | 为 true 时追加 `--cca` 标志(按 .smk 原样转发) |
| `Params.mimseq.remap` | `--remap` | true | 启用重比对 |
| `Params.mimseq.crosstalks` | `--crosstalks` | false | 启用修饰间串扰分析 |
| `ROOT_DIR` | — | `.` | 定位父目录 bin/mods/run.py |

`Params.mimseq` 段其余字段(species、cluster_id、max_mismatches、p_adj 等)
由 mimseq 其他阶段(align/clusters/coverage/deseq)模块消费,本规则不转发。
资源默认:`threads = 4`(经 `rule_res("mimseq_mods")` 可覆盖)、
`mem_mb = 32000`。引擎解释器为规则 conda 环境(`../mimseq.yaml`)的
`python`。

直接调用示例(共享引擎脚本 `modules/mimseq/bin/mods/run.py`):

```bash
python $ROOT_DIR/modules/mimseq/bin/mods/run.py \
  --out output/ --name tRNAseq --threads 4 \
  --min-cov 0.0005 --misinc-thresh 0.1 --remap
```

## 2. 算法原理与生物学意义

mim-tRNAseq:tRNA 修饰(如 m1A、m1G、m2,2G)在逆转录时引起错配或终止,
错配掺入率正比于修饰占据率。本模块复用 mimseq 包的 mmQuant(逐位点错配
定量)、ssAlign(结构比对)与 crosstalks(相邻修饰串扰)组件,按
`min_cov`/`misinc_thresh` 判定位点修饰,可选 remap(把错配 reads 重比对
以提高 tRNA 家族归属准确性)。生物学意义:tRNA 修饰动态调控翻译效率与
应激响应,位点级修饰谱是 mim-tRNAseq 的核心产出。结果阅读:定量表按 tRNA
位点给出错配率/修饰判定;`mods.done` 标记模块完成。

## 3. 参考文献

1. mimseq(mim-tRNAseq 分析工具)—— https://github.com/drewjbehren/mimseq
   (与父目录 `modules/mimseq/README.md` 引用一致)
