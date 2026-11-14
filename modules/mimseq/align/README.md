# mimseq/align — GSNAP 比对子模块(mim-tRNAseq)

C 形态子模块(单规则 `mimseq_align` + result 锚点):mim-tRNAseq 流程的
比对步,经共享父目录包装脚本 `modules/mimseq/bin/align/run.py` 调
GSNAP 做 SNP-tolerant 比对(可选两轮 remap)。流程全貌、共享
`Params.mimseq` 参数表与 bin/run.py 约定见父目录
`modules/mimseq/README.md`;conda 环境共享父目录 `mimseq.yaml`。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 说明 |
|---|---|
| `{outdir}/sample_data.tsv` | 样本—条件映射表(父模块 prepare_sample_data 产出) |
| `{outdir}/state/` | tRNAtools 产出的 pickle 状态目录(序列化状态) |

### 1.2 输出

```
{outdir}/align.done                      # touch 标记(output 锚点)
{outdir}/samples/{sample}/...            # 比对 BAM 等(run.py 内部产物,见父 README)
{outdir}/mimseq_align.sh                 # 回执
logdir/mimseq_align.log                  # 运行日志
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Params.mimseq.name` | `--name` | `tRNAseq` | 分析名称前缀 |
| `Params.mimseq.max_mismatches` | `--mismatches` | 0.075 | 最大错配率 |
| `Params.mimseq.remap` | `--remap` | true | 两轮比对(发现新修饰位点) |
| `Params.mimseq.keep_temp` | `--keep-temp` | false | 保留临时文件 |
| `ROOT_DIR` | — | `.` | 仓库根(拼 `modules/mimseq/bin/align/run.py`) |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录基路径 |

线程/内存经 `rule_res("mimseq_align")` 配置:threads 默认 10(传
`--threads`),mem_mb 默认 32000。运行的 Python 解释器为规则 conda 环境
(`../mimseq.yaml`)的 `python`。

引擎命令(回执脚本原文):

```bash
python modules/mimseq/bin/align/run.py \
    --sample-data {outdir}/sample_data.tsv --name tRNAseq --out {outdir}/ \
    --threads 10 --mismatches 0.075 --remap
```

注意:`--out` 必须带尾部 `/`(父 README 注意事项 2)。

## 2. 算法原理与生物学意义

纯工具封装:mimseq(https://github.com/drewjbehren/mimseq)面向
mim-tRNAseq——tRNA 修饰诱导错配测序。align 步用 GSNAP 把 reads 比对到
tRNA 参考(含 SNP 容错,容忍已知基因组 SNP 与修饰错配),`--remap` 第
一轮发现高错配位点后做第二轮比对以校正新修饰位点。生物学意义:tRNA
高度相似(isodecoder 仅差数碱基)且修饰造成系统性错配,SNP-tolerant
多轮比对是准确定量 tRNA 表达与修饰的前提(比对后由 clusters/mods/
coverage 子模块继续解析)。

## 3. 参考文献

1. mimseq(mim-tRNAseq 分析工具)— https://github.com/drewjbehren/mimseq
