# mimseq/polygenomes/coverage — 多基因组 tRNA 覆盖度分析子模块

C 形态子模块(单规则 `mimseq_coverage` + result 锚点):mim-tRNAseq 流程
的覆盖度步,按 genome 调共享父目录包装脚本
`modules/mimseq/bin/coverage/run.py` 计算 tRNA 覆盖度并出图。流程全貌与
共享 `Params.mimseq` 参数表见父目录 `modules/mimseq/README.md`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{outdir}/{genome}/mods.done` | 上游 mods 子模块(polygenomes 布局)完成标记 |

### 1.2 输出

```
{outdir}/{genome}/coverage.done          # touch 标记(output 锚点)
{outdir}/{genome}/coverage_byaa.txt      # 按氨基酸聚合覆盖度(run.py 产物)
{outdir}/{genome}/coverage_bygene.txt    # 按基因聚合覆盖度
{outdir}/{genome}/coverage.log           # 运行日志(run.py 内部)
{outdir}/{genome}/state/coverageData.pkl # 覆盖度状态 pickle
{outdir}/{genome}/mimseq_coverage.sh     # 回执
logdir/{genome}/mimseq_coverage.log      # 规则日志
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.mimseq.control_cond` | `--control-cond` | "" | 对照条件名(非空才传) |
| `Params.mimseq.misinc_thresh` | `--misinc-thresh` | 0.1 | 错配(错参)阈值 |
| `Params.mimseq.no_cca` | `--cca` | false | false 时**加** `--cca`(做 CCA 分析) |
| `Params.mimseq.double_cca` | `--double-cca` | false | 启用双 CCA |
| `genome.{genome}.mito_trnas` | `--mito-trnas` | "" | 线粒体 tRNA FASTA(每 genome,非空才传) |
| `ROOT_DIR` | — | `.` | 仓库根(拼 `modules/mimseq/bin/coverage/run.py`) |
| `outdir` / `logdir` | — | `output` / `log` | 目录基路径(输出按 `{genome}/` 分层) |

线程/内存经 `rule_res("mimseq_coverage")` 配置:threads 默认 4,mem_mb
默认 16000。运行的 Python 解释器与 `PATH` 在 .smk 中硬编码(站点路径)。

引擎命令(回执脚本原文):

```bash
python modules/mimseq/bin/coverage/run.py \
    --out {outdir}/{genome}/ --misinc-thresh 0.1 [ --control-cond <cond> ] \
    [ --mito-trnas <fa> ] [ --cca ] [ --double-cca ]
```

## 2. 算法原理与生物学意义

纯工具封装:coverage 步汇总各样本 BAM 在每个 tRNA 转录本各位置的覆盖,
按氨基酸/基因两个层次聚合输出(供表达热力图等),并可做 CCA(成熟
tRNA 3' 端 CCA 尾)状态分析;`misinc_thresh` 以下的错配忽略以免修饰
错配虚增覆盖差异。生物学意义:tRNA 覆盖度反映同工受体/isodecoder 的
表达量与 3' 端成熟度;多 genome 布局下每个参考的 tRNA 注释独立计算,
互不串扰。

## 3. 参考文献

1. mimseq(mim-tRNAseq 分析工具)— https://github.com/drewjbehren/mimseq
