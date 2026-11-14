# mimseq/coverage — tRNA 修饰与覆盖度分析子模块(mimseq coverage)

mim-tRNAseq 流程的覆盖度/修饰步骤:基于上游 mods.done 产物计算 tRNA 覆盖度
与错义掺入(misincorporation)轮廓并出图。C 形态子模块,共享父目录脚本
`modules/mimseq/bin/coverage/run.py`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `outdir/mods.done` | 上游 mimseq 修饰检测步骤的完成标记(DAG 锚点);分析中间产物均在 `outdir/` 内由脚本自行读取 |
| `genome.mito_trnas` | 线粒体 tRNA 注释(可选,透传 `--mito-trnas`) |

### 1.2 输出

```
{outdir}/coverage.done              # 完成标记(touch);图表产物由 run.py 写入 {outdir}/
```

日志:`{logdir}/mimseq_coverage.log`。

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `outdir` / `logdir` | `--out` | `output` / `log` | 输出/日志目录 |
| `Params.mimseq.misinc_thresh` | `--misinc-thresh` | `0.1` | 错义掺入阈值 |
| `Params.mimseq.control_cond` | `--control-cond` | `""` | 对照条件名,非空才传(json 示例键为 `control_condition`,本 .smk 读 `control_cond`,以 .smk 为准) |
| `genome.mito_trnas` | `--mito-trnas` | `""` | 非空才传 |
| `Params.mimseq.no_cca` | `--cca` | `false` | 默认传 `--cca`;`no_cca=true` 时不传 |
| `Params.mimseq.double_cca` | `--double-cca` | `false` | true 时传 |
| `ROOT_DIR` | — | `.` | 拼 `modules/mimseq/bin/coverage/run.py` |
| 线程/内存 | — | `threads=4`,`mem_mb=16000` | `rule_res("mimseq_coverage")` 可覆盖 |

环境:run 块内固定激活 conda env
`/data/pub/zhousha/env/mutation_0.1/eda061b3f191779ad16ff11ee6fe53b4_`
并以其 python 解释器运行(该绝对路径为环境现状,迁移时需替换)。

直接调用示例(等价引擎命令):

```bash
python modules/mimseq/bin/coverage/run.py \
  --out output/ --misinc-thresh 0.1 [--control-cond treated] \
  [--mito-trnas mito_tRNAs.bed] --cca [--double-cca]
```

## 2. 算法原理与生物学意义

mim-tRNAseq(酶促捕获 tRNA 测序)用改良的 YAMAT/mimseq 策略定量 tRNA 表达
并检测碱基修饰导致的错义掺入/截断信号。coverage 步骤按条件汇总每 tRNA/
同受体(isoacceptor)的覆盖度与 misincorporation 率:`misinc_thresh` 以上
的位点视为修饰候选,`control_cond` 用于条件间差异比较;CCA/双 CCA 选项
处理 tRNA 3' 端成熟态(转录后加 CCA)。tRNA 修饰动态(如 m5C、m1A)与
翻译调控、应激与衰老密切相关。结果怎么读:输出目录内覆盖度图看表达丰度
与覆盖均一性,修饰位点表看 misincorporation 是否超阈值。

## 3. 参考文献

1. mimseq / mim-tRNAseq — https://github.com/goodarzilab/mimseq
   (Behrens A et al., *New insights into tRNA biology using mim-tRNAseq*,
   2021;如引用请以原文为准)
