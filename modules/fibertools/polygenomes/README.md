# fibertools/polygenomes — fibertools-rs 元素调用流水线子模块

C 形态子模块(四规则链:`ft_predict_m6a` → `ft_add_nucleosomes` →
`ft_fire` → `ft_extract`):用 fibertools-rs(`ft`)从 PacBio HiFi
kinetics BAM 依次调出 m6A、核小体、FIRE,并提取元素 BED。conda 环境共享
父目录 `../fibertools.yaml`(见 `modules/fibertools/`)。

## 1. 模块输入输出与参数

### 1.1 输入

| 输入 | 规则 | 说明 |
|---|---|---|
| `{indir}/{sample_id}/{sample_id}.bam` | ft_predict_m6a | PacBio CCS BAM(带 polymerase kinetics IPD/PL 标签) |
| 上一步产物 | 其余规则 | 逐规则串联,见输出 |

({genome} wildcard 只出现在 log 路径;BAM 布局不随 genome 分层。)

### 1.2 输出

```
{outdir}/{sample_id}/{sample_id}.fiberseq.bam       # 带 m6A 调用(MM/ML 标签)
{outdir}/{sample_id}/{sample_id}.fiberseq.nuc.bam   # 加核小体/MSP 调用
{outdir}/{sample_id}/{sample_id}.fiberseq.fire.bam  # 加 FIRE 调用
{outdir}/{sample_id}/{sample_id}.m6a.bed.gz         # ft extract:m6A BED12
{outdir}/{sample_id}/{sample_id}.nuc.bed.gz         # 核小体 BED12
{outdir}/{sample_id}/{sample_id}.msp.bed.gz         # MSP BED12
{outdir}/{sample_id}/{sample_id}.fire.bed.gz        # FIRE BED12
{outdir}/{sample_id}/ft_predict_m6a.sh 等            # 每规则回执
logdir/{genome}/{sample_id}/ft_predict_m6a.log 等    # 每规则日志
```

### 1.3 参数

| 参数 | CLI | 默认 | 说明 |
|---|---|---|---|
| `Procedure.fibertools` | 可执行文件 | null(回退 `ft`) | fibertools-rs 可执行路径 |
| `Params.fibertools.ont` | `ft fire --ont` | false | ONT 数据时加 `--ont`(仅 ft_fire) |
| `samples` | — | [] | 样本列表 |
| `indir` / `outdir` / `logdir` | — | `input` / `output` / `log` | 目录基路径 |

线程/内存经 `rule_res(...)` 配置:predict-m6a 16/32000,add-nucleosomes
16/8000,fire 8/8000,extract 8/8000。

引擎命令(每样本,回执脚本原文):

```bash
ft predict-m6a -t 16 {in.bam} {s}.fiberseq.bam
ft add-nucleosomes -t 16 {s}.fiberseq.bam {s}.fiberseq.nuc.bam
ft fire {s}.fiberseq.nuc.bam {s}.fiberseq.fire.bam [ --ont ]
ft extract -t 8 --m6a {s}.m6a.bed.gz --nuc {s}.nuc.bed.gz --msp {s}.msp.bed.gz --fire {s}.fire.bed.gz {s}.fiberseq.fire.bam
```

## 2. 算法原理与生物学意义

纯工具封装(fibertools-rs):`predict-m6a` 从 HiFi 聚合酶动力学(IPD/PL)
预测每条分子上的 m6A 位点;`add-nucleosomes` 按 m6A 足迹模式判读核小体
占位与之间裸露的 MSP(甲基化敏感区段);`fire` 用机器学习把 MSP 分类为
FIRE 调控元件(写入 Ma/Aq 标签,precision/q 约定见
`modules/fiberbrowser/README.md`);`extract` 把四类元素落成逐纤维 BED12。
生物学意义:该 BAM/BED 链是 Fiber-seq 单分子染色质分析的上游标准产物,
下游 fiberseq/ 子模块(quant、motif、coactuation 等)都消费
`{m6a,nuc,msp,fire}.bed.gz`。

## 3. 参考文献

1. fibertools-rs(`ft`)— https://github.com/fiberseq/fibertools-rs ;
   文档 https://docs.rs ;项目主页 https://fiberseq.github.io/
2. FIRE analysis pipeline — https://github.com/fiberseq/FIRE
